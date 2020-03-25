library(tidyverse)

######################
## Format counts - all
######################

dat <- read_tsv("data_raw/counts/Shah.featurecounts.paired.tsv", skip=1)

#Create simplified names list
name.list <- sub("genomics/results/bam_filter_paired/",
                 "", colnames(dat)[7:length(colnames(dat))])
name.list <- sub("_filter_paired.bam",
                 "", name.list)

#Format counts
dat.format <- dat %>% 
  select(-c(Chr:Length)) %>% 
  #Rename columns
  setNames(c("geneName", name.list))

write_csv(dat.format, "data_clean/Shah.counts.clean.csv")
