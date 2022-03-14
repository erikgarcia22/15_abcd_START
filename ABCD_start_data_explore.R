#
# Erik Garcia ABCD START project----

# Written by Erik J. Garcia 

# Introduction this script will begin to explore the ABCD dataset

# Load Packages

library(tidyverse)

# Get data from external hard drive into R
  
  drug_use_disorders_data <- read_delim(
    file = "E:/ABCDStudyNDAr4/drug_use_disorders01.txt"
  ) 

  youth_substance_use_attitudes_data <- read_delim(
    file = "E:/ABCDStudyNDAr4/abcd_ysua01.txt"
  )

# Need to remove the Data information row
  
  drug_use_disorders_data <- data[-1, ]

  write_csv(data, file = "drug_use_disorders01.csv")
  