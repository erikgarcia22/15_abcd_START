# practice data----
#

library(bigsnpr)

# Read from bed/bim/fam----
  snp_readBed("tmp-data/public-data3.bed")

  obj.bigSNP <- snp_attach("tmp-data/public-data3.rds")

# see how the object looks in compacted view

  str(obj.bigSNP, max.level = 2, strict.width = 'cut')


# Get alias for useful slots----

  G <- obj.bigSNP$genotypes
  CHR <- obj.bigSNP$map$chromosome
  POS <- obj.bigSNP$map$physical.pos
  y <- obj.bigSNP$fam$affection

NCORES <- nb_cores()

# Read external summary statistics----

  sumstats <- bigreadr::fread2("tmp-data/public-data3-sumstats.txt")

  str(sumstats)

# Split the genotype data into 350 for validation and remaining as test set----

  set.seed(1)
  ind.val <- sample(nrow(G), 350)
  ind.test <- setdiff(rows_along(G), ind.val)

# Matching variants between genotype data and summary stats
# To match variants contained in genotype data and summary 
# statistics, the variables "chr" (chromosome number), "pos" 
# (genetic position), "a0" (reference allele) and "a1" (derived allele)
# should be available in the summary statistics and in the genotype data. 
# These 4 variables are used to match variants between the two data frames. 
# From the summary statistics, you need to get "beta", "beta_se" 
# (standard errors), and "n_eff" (effective sample size per variant for GWAS 
# with logistic regression, and just total sample size for continuous traits).
  
  # sumstats$n_eff <- 4 /(1 / sumstats$n_case + 1 / sumstats$n_control)
  # sumstats$n_case <- sumstats$n_control <- NULL
  
    sumstats$n_eff <- sumstats$N
    map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
    df_beta <- snp_match(sumstats, map)
  
# Problem with matching bec different genome builds
# can either convert between builds with snp_modifyBuild()
# (or directly use the converted positions in info), or match by rsIDs instead.
    
    df_beta <- snp_match(sumstats, map, join_by_pos = FALSE)  # use rsid instead of pos

# If no or few variants are actually flipped, you might want to disable 
# the strand flipping option (strand_flip = FALSE) 
    
# Computing LDpred2 scores genome-wide----
    
# Correlation, compute correlations between variants.----
    # POS2 <- snp_asGeneticPos(CHR, POS, dir = "tmp-data", ncores = NCORES)
    # To avoid downloading "large" files, this has been precomputed
    POS2 <- obj.bigSNP$map$genetic.dist    

# We create the on-disk sparse genome-wide correlation matrix on-the-fly:
    tmp <- tempfile(tmpdir = "tmp-data")
    
    for (chr in 1:22) {
      
      # print(chr)
      
      ## indices in 'df_beta'
      ind.chr <- which(df_beta$chr == chr)
      ## indices in 'G'
      ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
      
      corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3 / 1000,
                       infos.pos = POS2[ind.chr2], ncores = NCORES)
      
      if (chr == 1) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp, compact = TRUE)
      } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
      }
    }  

#check file size
    file.size(corr$sbk) / 1024^3  # file size in GB  
    

# LDpred2-inf: infinitesimal model----
    (ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                    sample_size = n_eff, blocks = NULL)))
    
    h2_est <- ldsc[["h2"]]
    beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
    pred_inf <- big_prodVec(G, beta_inf, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
    cor(pred_inf, y[ind.test])

# LDpred2(-grid): grid of models----
    
    (h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4))
    (p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2))
    (params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE)))
    
    # takes less than 2 min with 4 cores
    beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)
    pred_grid <- big_prodMat(G, beta_grid, ind.col = df_beta[["_NUM_ID_"]])
    params$score <- apply(pred_grid[ind.val, ], 2, function(x) {
      if (all(is.na(x))) return(NA)
      summary(lm(y[ind.val] ~ x))$coef["x", 3]
      # summary(glm(y[ind.val] ~ x, family = "binomial"))$coef["x", 3]
    })
    
    library(ggplot2)
    ggplot(params, aes(x = p, y = score, color = as.factor(h2))) +
      theme_bigstatsr() +
      geom_point() +
      geom_line() +
      scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
      facet_wrap(~ sparse, labeller = label_both) +
      labs(y = "GLM Z-Score", color = "h2") +
      theme(legend.position = "top", panel.spacing = unit(1, "lines"))
    
    library(dplyr)
    params %>%
      mutate(sparsity = colMeans(beta_grid == 0), id = row_number()) %>%
      arrange(desc(score)) %>%
      mutate_at(c("score", "sparsity"), round, digits = 3) %>%
      slice(1:10)
    
    best_beta_grid <- params %>%
      mutate(id = row_number()) %>%
      # filter(sparse) %>% 
      arrange(desc(score)) %>%
      slice(1) %>%
      pull(id) %>%
      beta_grid[, .]
    
    pred <- big_prodVec(G, best_beta_grid, ind.row = ind.test,
                        ind.col = df_beta[["_NUM_ID_"]])
    cor(pred, y[ind.test])
    
# LDpred2-auto: automatic model----
    # takes less than 2 min with 4 cores
    multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                                   vec_p_init = seq_log(1e-4, 0.5, length.out = 30),
                                   allow_jump_sign = FALSE, shrink_corr = 0.95,
                                   ncores = NCORES)
    str(multi_auto, max.level = 1)
    
    str(multi_auto[[1]], max.level = 1)
 
    library(ggplot2)
    auto <- multi_auto[[1]]
    plot_grid(
      qplot(y = auto$path_p_est) +
        theme_bigstatsr() +
        geom_hline(yintercept = auto$p_est, col = "blue") +
        scale_y_log10() +
        labs(y = "p"),
      qplot(y = auto$path_h2_est) +
        theme_bigstatsr() +
        geom_hline(yintercept = auto$h2_est, col = "blue") +
        labs(y = "h2"),
      ncol = 1, align = "hv"
    )   
    
    beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
    pred_auto <- big_prodMat(G, beta_auto, ind.row = ind.test, 
                             ind.col = df_beta[["_NUM_ID_"]])
    
    sc <- apply(pred_auto, 2, sd)
    keep <- abs(sc - median(sc)) < 3 * mad(sc)
    final_beta_auto <- rowMeans(beta_auto[, keep])
    final_pred_auto <- rowMeans(pred_auto[, keep])
    
    cor(final_pred_auto, y[ind.test])
    