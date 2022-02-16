#
# Erik Garcia ABCD START project----

# Written by Erik J. Garcia 

# Introduction this script will begin to explore the ABCD dataset

# Load Packages

library(tidyverse)

# Get data
  
  data <- read_delim(
    file = "E:/ABCDStudyNDAr4/drug_use_disorders01.txt"
  ) 

# Need to remove the Data information row
  
  data <- data[-1, ]

  write_csv(data, file = "drug_use_disorders01.csv")
  