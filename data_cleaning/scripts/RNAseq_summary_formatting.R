library(tidyverse)

##### Raw and adapter trimmmed seqs ##### 
trim.files <- list.files("data_raw/fastq_trim/", pattern="*settings", 
                         all.files=FALSE, full.names=TRUE) %>% 
  gsub("//", "/", .)

#Create simplified names list
name.list <- sub("data_raw/fastq_trim/",
                 "", trim.files)
name.list <- sub(".settings",
                 "", name.list)
raw.summ <- data.frame()

for (file in trim.files){
  sampID <- gsub("data_raw/fastq_trim/", "", file) %>% 
    gsub(".settings", "", .)
  
  raw.summ.temp <- read_tsv(file, col_names = FALSE) %>% 
    filter(startsWith(X1, 'Total number of read pairs') |
             startsWith(X1, 'Number of retained reads')) %>% 
    separate(X1, into=c("metric","value"), sep=": ") %>% 
    mutate(sampID = sampID) %>% 
    pivot_wider(names_from = "metric", values_from = "value") %>% 
    rename(raw="Total number of read pairs",
           trim="Number of retained reads") %>% 
    mutate_at(vars(raw, trim), as.numeric) %>% 
    mutate(raw=raw*2)
  
  raw.summ <- bind_rows(raw.summ, raw.summ.temp)
}

##### Align summary ##### 
aligned <- read_tsv("results/results_cleaning/summary.alignment.tsv", col_names = FALSE) %>% 
  #Get sample name from filename
  separate(X1, into=c("a","b","c","d","e","f","g"), sep="/") %>% 
  select(a,g) %>% 
  separate(g, into="sampID", sep="[_]Align") %>% 
  #Give name to all data from that sample
  fill(sampID) %>% 
  #Separate data to column
  separate(a, into=c("h","i"), sep=" \\+ ") %>% 
  drop_na(i) %>% 
  
  #Remove failed values
  separate(i, into=c("fail","i"), sep="[0-9] ") %>% 
  select(-fail) %>% 
  
  #Recode data types (f)
  separate(i, into=c("i"), sep="[(]") %>% 
  mutate(i = fct_recode(factor(i),
                        to.be.aligned="in total ",
                        secondary.align="secondary",
                        chimeric.align="supplementary",
                        PCR.dups="duplicates",
                        align="mapped ",
                        paired="paired in sequencing",
                        R1.paired="read1", R2.paired="read2",
                        align.paired="properly paired ",
                        both.align.paired= "with itself and mate mapped",
                        one.align.paired="singletons " ,
                        both.align.paired.diffCHR="with mate mapped to a different chr",
                        both.align.paired.diffCHR.mapq="with mate mapped to a different chr ")) %>% 
  pivot_wider(names_from = "i", values_from = "h") %>% 
  mutate_at(vars(-sampID), as.numeric)

##### BAM summary ##### 
bam.colnames <- read_delim("results/results_cleaning/bam.metrics.tsv", 
                           delim = " ", col_names=FALSE)[7,1] %>% 
  separate(X1, into=as.character(c(1:30)), sep="\t") %>% 
  unlist(use.names = FALSE)

bam.summ <- read_delim("results/results_cleaning/bam.metrics.tsv", 
                       delim = " ", col_names=FALSE, comment = "#") %>%
  #Get sample name from file name
  separate(X1, into=c("dir1","dir2","dir3","dir4","dir5","dir6","file"), sep="/") %>% 
  select(dir1, file) %>% 
  separate(file, into="sampID", sep="[_]Align") %>% 
  #Fill in sampID
  fill(sampID) %>% 
  #Sep data columns
  separate(dir1, into = bam.colnames, sep="\t") %>% 
  #Remove histogram columns
  drop_na(CODING_BASES) %>% 
  #Keep only data rows
  filter(!startsWith(as.character(PF_BASES), 'PF')) %>% 
  #Make numeric
  mutate_at(vars(-sampID), as.numeric) %>% 
  #Calculate perc align
  mutate(PCT_PF_ALIGNED = PF_ALIGNED_BASES/PF_BASES) %>% 
  select(sampID, PCT_PF_ALIGNED, everything())

##### Filtered alignment ##### 
align.filter.summ <- read_tsv("results/results_cleaning/summary.align.filter.paired.tsv",
                              col_names = FALSE) %>% 
  #Get sample name from filename
  separate(X1, into=c("a","b","c","d"), sep="/") %>% 
  select(a,d) %>% 
  separate(d, into="sampID", sep="[_]filter") %>% 
  #Give name to all data from that sample
  fill(sampID) %>% 
  #Separate data to column
  separate(a, into=c("h","i"), sep=" \\+ ") %>% 
  drop_na(i) %>% 
  
  #Remove failed values
  separate(i, into=c("fail","i"), sep="[0-9] ") %>% 
  select(-fail) %>% 
  
  #Recode data types (f)
  separate(i, into=c("i"), sep="[(]") %>% 
  mutate(i = fct_recode(factor(i),
                        to.be.aligned="in total ",
                        secondary.align="secondary",
                        chimeric.align="supplementary",
                        PCR.dups="duplicates",
                        align="mapped ",
                        paired="paired in sequencing",
                        R1.paired="read1", R2.paired="read2",
                        align.paired="properly paired ",
                        both.align.paired= "with itself and mate mapped",
                        one.align.paired="singletons " ,
                        both.align.paired.diffCHR="with mate mapped to a different chr",
                        both.align.paired.diffCHR.mapq="with mate mapped to a different chr ")) %>% 
  pivot_wider(names_from = "i", values_from = "h") %>% 
  mutate_at(vars(-sampID), as.numeric) %>% 
  rename_if(is.numeric, ~paste(., "filter", sep="_"))

##### Gene counts ##### 
gene <- read_tsv("data_raw/counts/Shah.featurecounts.paired.tsv.summary") %>% 
  #Rename columns
  setNames(c("Status", name.list)) %>% 
  #Keep rows with non-zero values
  filter_if(is.numeric, any_vars(.>0)) %>% 
  #Transpose
  pivot_longer(-Status, names_to = "sampID", values_to = "value") %>% 
  pivot_wider(names_from = Status, values_from = value) 

##### MERGE ##### 
summ.all <- full_join(raw.summ, aligned, by="sampID") %>% 
  full_join(bam.summ, by="sampID") %>% 
  full_join(align.filter.summ, by="sampID") %>% 
  full_join(gene, by="sampID") %>% 

  #remove columns that are all blank or 0
  mutate_at(vars(-sampID), as.numeric) %>% 
  select_if(~sum(!is.na(.)) > 0) %>% 
  select_if(~sum(.!= 0) > 0)

dir.create("data_clean", showWarnings = FALSE)
write_csv(summ.all, "data_clean/Shah.data.cleaning.metrics.csv")

#### plots #####

seq.summ <- summ.all %>% 
  select(sampID, raw, trim, 
         both.align.paired,
         both.align.paired_filter,
         Assigned) %>% 
  
  rename(pair.align= both.align.paired,
         pair.align.filter= both.align.paired_filter,
         assign= Assigned) %>%
  
  pivot_longer(-sampID, names_to = "group", values_to = "sequences") %>% 
  mutate(group = factor(group, levels=c("raw", "trim", "pair.align",
                                        "pair.align.filter", "assign")))

seq.summ.plot <- seq.summ %>%   
  ggplot(aes(x=sampID, y=sequences, fill=group)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_discrete(name="", labels=c("Raw",
                                        "Min quality, adapter trimmed",
                                        "Aligned, paired",
                                        "High-quality aligned, paired",
                                        "Assigned to gene")) +
  labs(x="",y="Total sequences")

dir.create("figs", showWarnings = FALSE)
dir.create("figs/cleaning", showWarnings = FALSE)
ggsave("figs/cleaning/total.seqs.cleaning.png", seq.summ.plot,
       height=5, width=10)


seq.summ.plot2 <- seq.summ %>%
  group_by(sampID) %>% 
  mutate(pct = sequences/sequences[group == "raw"]*100) %>% 
  
  ggplot(aes(x=sampID, y=pct, fill=group)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_discrete(name="", labels=c("Raw",
                                        "Min quality, adapter trimmed",
                                        "Aligned, paired",
                                        "High-quality aligned, paired",
                                        "Assigned to gene")) +
  labs(x="",y="Percent of raw sequences")

ggsave("figs/cleaning/pct.seqs.cleaning.png", seq.summ.plot2,
       height=5, width=10)


align.cv <- summ.all %>% 
  select(sampID, PCT_PF_ALIGNED, MEDIAN_CV_COVERAGE) %>%
  
  ggplot(aes(x=MEDIAN_CV_COVERAGE, y=PCT_PF_ALIGNED)) +
  geom_point(aes(color=sampID), size=2) +
  theme_classic() +
  lims(x=c(0,1), y=c(0,1)) +
  labs(x="Median CV coverage", y="Percent alignment") +
  coord_fixed()

ggsave("figs/cleaning/bam.metrics.png", align.cv,
       height=5, width=6)

align.cv2 <- summ.all %>% 
  select(sampID, PCT_PF_ALIGNED, MEDIAN_CV_COVERAGE) %>%
  
  ggplot(aes(x=MEDIAN_CV_COVERAGE, y=PCT_PF_ALIGNED)) +
  geom_point(aes(color=sampID), size=2) +
  theme_classic() +
  labs(x="Median CV coverage", y="Percent alignment")

ggsave("figs/cleaning/bam.metrics_scaled.png", align.cv2,
       height=5, width=6)
