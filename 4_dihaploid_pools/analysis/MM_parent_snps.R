#' ---
#' title: "MM Parent SNPs"
#' author: Kirk Amundson
#' date: 2020_1007
#' output: html_notebook
#' ---
#' 
#' Aim: Define parent-specific SNP loci and inducer/non-inducer specific alleles
#' at these loci for low-pass SNP analysis of MM dihaploid cohorts.
#' 
#' Low quality sites were filtered out in the preceding step based on attributes
#' of that site across all samples. Here, I implement sample-specific filters
#' to generate a flat tsv of parent-informative SNP loci to use in binned genotyping
#' for chromosome parental origin tests.
#' 
#' ## Packages:
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)

#' 
#' ## Functions:
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
# separate sample-specific VCF attributes
sep <- function(...) {
  dots <- list(...)
  separate_(..., into = paste(dots[[2]], attributes[[1]], sep = "_"), convert = T, sep = ":")
}

## ------------------------------------------------------------------------------------------------------------------------------------------------------------
# filter to retain only those loci with called homozygous genotypes for alternate alleles in two specified samples
filter_homozygous_vars <- function(tetraploid, hi, snplist, dp_threshold) {
  nhi_gt <- enexpr(tetraploid)
  hi_gt <- enexpr(hi)
  
  nhi_dp <- str_replace(nhi_gt, "GT", "DP")
  hi_dp  <- str_replace(hi_gt, "GT", "DP")
  
  hom <- snplist %>%
    filter(!!nhi_gt %in% c("0/0/0/0", "1/1/1/1")) %>%
    filter(!!nhi_dp >= dp_threshold) %>% 
    filter(!!hi_gt %in% c("0/0", "1/1")) %>%
    filter(!!hi_dp >= dp_threshold) %>% 
    filter(!(!!nhi_gt == "0/0/0/0" & !!hi_gt == "0/0")) %>%
    filter(!(!!nhi_gt == "1/1/1/1" & !!hi_gt == "1/1"))

  plt <- hom %>%
    ggplot(., aes(x = POS)) +
    geom_histogram(binwidth = 1e6) +
    facet_wrap(~CHROM, strip.position = "r", nrow = 7) +
    scale_y_log10() +
    theme_bw()

  print(plt)
  print(paste(nrow(hom), "SNP retained", sep = "_"))
  return(hom)
}

#' 
#' ## Read in data:
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
# either download files to local or mount to server via, e.g., sshfs
files <- dir(pattern = "-filtered-",
             path = "../data/calls/",
             full.names = T)
files

#' 
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
snps <- map_dfr(files, function(x) read_tsv(x, col_names = T, na = "NA"))

#' 
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
names(table(snps$FORMAT)) # should only have one entry. does.
attributes <- str_split(names(table(snps$FORMAT[1])), ":", simplify = F)
attributes[[1]]
sample_vars <- colnames(snps)[-c(seq(1,49), ncol(snps))]

#' 
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
snps2 <- snps %>%
  Reduce(f = sep, x = sample_vars)

#' 
#' WA.077 x IVP101 (n=50)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(WA077_GT, IVP101_GT, snps2, 10) %>% 
  mutate(WA077 = ifelse(WA077_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, IVP101, WA077) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., "WA077-IVP101-hom-SNP.tsv", col_names = T)

#' 
#' WA.077 x IVP35 (n=134)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(WA077_GT, IVP35_GT, snps2, 10) %>% 
  mutate(WA077 = ifelse(WA077_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, IVP35, WA077) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., "WA077-IVP35-hom-SNP.tsv", col_names = T)

#' 
#' WA.077 x PL4 (n=107)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(WA077_GT, PL4_GT, snps2, 10) %>% 
  mutate(WA077 = ifelse(WA077_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, PL4, WA077) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., "WA077-PL4-hom-SNP.tsv", col_names = T)

#' 
#' LR00.014 x IVP101 (n=30)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(LR00014_GT, IVP101_GT, snps2, 10) %>% 
  mutate(LR00014 = ifelse(LR00014_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT)) %>% 
  dplyr::select(CHROM, POS, REF, IVP101, LR00014) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., "LR00014-IVP101-hom-SNP.tsv", col_names = T)

#' 
#' LR00.014 x IVP35 (n=77)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(LR00014_GT, IVP35_GT, snps2, 10) %>% 
  mutate(LR00014 = ifelse(LR00014_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, IVP35, LR00014) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., "LR00014-IVP35-hom-SNP.tsv", col_names = T)

#' 
#' LR00.014 x PL4 (n=67)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(LR00014_GT, PL4_GT, snps2, 10) %>% 
  mutate(LR00014 = ifelse(LR00014_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(PL4 = ifelse(IVP35_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, IVP35, LR00014) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., "LR00014-PL4-hom-SNP.tsv", col_names = T)

#' 
#' LR00.026 x IVP101 (n=4)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(LR00026_GT, IVP101_GT, snps2, 10) %>% 
  mutate(LR00026 = ifelse(LR00026_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, IVP101, LR00026) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., "LR00026-IVP101-hom-SNP.tsv", col_names = T) 

#' 
#' LR00.026 x IVP35 (n=36)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(LR00026_GT, IVP35_GT, snps2, 10) %>% 
  mutate(LR00026 = ifelse(LR00026_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, IVP35, LR00026) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., "LR00026-IVP35-hom-SNP.tsv", col_names = T)

#' 
#' LR00.026 x PL4 (n=35)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(LR00026_GT, PL4_GT, snps2, 10) %>% 
  mutate(LR00026 = ifelse(LR00026_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, PL4, LR00026) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., "LR00026-PL4-hom-SNP.tsv", col_names = T)

#' 
#' Atlantic x IVP35 (n=5)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(Atlantic_GT, IVP35_GT, snps2, 10) %>% 
  mutate(Atlantic = ifelse(Atlantic_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, IVP35, Atlantic) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., "Atlantic-IVP35-hom-SNP.tsv", col_names = T) 

#' 
#' Atlantic x PL4 (n=10)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(Atlantic_GT, PL4_GT, snps2, 10) %>% 
  mutate(Atlantic = ifelse(Atlantic_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, PL4, Atlantic) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., "Atlantic-PL4-hom-SNP.tsv", col_names = T)

#' 
#' Desiree x IVP101 (n=2)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(Desiree_GT, IVP101_GT, snps2, 10) %>% 
  mutate(Desiree = ifelse(Desiree_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, IVP101, Desiree) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>%
  write_tsv(., "Desiree-IVP101-hom-SNP.tsv", col_names = T)

#' 
#' Desiree x IVP35 (n=2)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(Desiree_GT, IVP35_GT, snps2, 10) %>% 
  mutate(Desiree = ifelse(Desiree_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, IVP35, Desiree) %>% 
  filter(CHROM %in% sprintf("chr%0.2d", 1:12)) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., "Desiree-IVP35-hom-SNP.tsv", col_names = T)

#' 
#' Desiree x PL4 (n=6)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(Desiree_GT, PL4_GT, snps2, 10) %>% 
  mutate(Desiree = ifelse(Desiree_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, PL4, Desiree) %>% 
  filter(CHROM %in% sprintf("chr%0.2d", 1:12)) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., "Desiree-PL4-hom-SNP.tsv", col_names = T)

#' 
#' 93.003 x IVP101 (n=12)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
clean_93003_dihaploids_IVP101 <- parent_snps(snps2, clean_93003_dihaploids_GT, clean_93003_dihaploids_DP, IVP101_GT, IVP101_DP) %>%
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT),
         clean_93003_dihaploids = ifelse(IVP101_GT == "0/0", ALT, REF)) %>%
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP101, clean_93003_dihaploids)
write_tsv(clean_93003_dihaploids_IVP101, "clean_93003_dihaploids-IVP101-SNP.tsv", col_names = T)

#' 
#' 93.003 x IVP35 (n=21)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
clean_93003_dihaploids_IVP35 <- parent_snps(snps2, clean_93003_dihaploids_GT, clean_93003_dihaploids_DP, IVP35_GT, IVP35_DP) %>%
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT),
         clean_93003_dihaploids = ifelse(IVP35_GT == "0/0", ALT, REF)) %>%
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP35, clean_93003_dihaploids)
write_tsv(clean_93003_dihaploids_IVP35, "clean_93003_dihaploids-IVP35-SNP.tsv", col_names = T)

#' 
#' 93.003 x PL4 (n=49)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
clean_93003_dihaploids_PL4 <- parent_snps(snps2, clean_93003_dihaploids_GT, clean_93003_dihaploids_DP, PL4_GT, PL4_DP) %>%
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT),
         clean_93003_dihaploids = ifelse(PL4_GT == "0/0", ALT, REF)) %>%
  select(Chrom = CHROM, Pos = POS, Ref = REF, PL4, clean_93003_dihaploids)
write_tsv(clean_93003_dihaploids_PL4, "clean_93003_dihaploids-PL4-SNP.tsv", col_names = T)

#' 
#' C91.640 x IVP101 (n=0)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
clean_C91640_dihaploids_IVP101 <- parent_snps(snps2, clean_C91640_dihaploids_GT, clean_C91640_dihaploids_DP, IVP101_GT, IVP101_DP) %>%
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT),
         clean_C91640_dihaploids = ifelse(IVP101_GT == "0/0", ALT, REF)) %>%
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP101, clean_C91640_dihaploids)
write_tsv(clean_C91640_dihaploids_IVP101, "clean_C91640_dihaploids-IVP101-SNP.tsv", col_names = T)

#' 
#' C91.640 x IVP35 (n=1)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
clean_C91640_dihaploids_IVP35 <- parent_snps(snps2, clean_C91640_dihaploids_GT, clean_C91640_dihaploids_DP, IVP35_GT, IVP35_DP) %>%
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT),
         clean_C91640_dihaploids = ifelse(IVP35_GT == "0/0", ALT, REF)) %>%
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP35, clean_C91640_dihaploids)
write_tsv(clean_C91640_dihaploids_IVP35, "clean_C91640_dihaploids-IVP35-SNP.tsv", col_names = T)

#' 
#' C91.640 x PL4 (n=86)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
clean_C91640_dihaploids_PL4 <- parent_snps(snps2, clean_C91640_dihaploids_GT, clean_C91640_dihaploids_DP, PL4_GT, PL4_DP) %>%
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT),
         clean_C91640_dihaploids = ifelse(PL4_GT == "0/0", ALT, REF)) %>%
  select(Chrom = CHROM, Pos = POS, Ref = REF, PL4, clean_C91640_dihaploids)
write_tsv(clean_C91640_dihaploids_PL4, "clean_C91640_dihaploids-PL4-SNP.tsv", col_names = T)

#' 
#' C93.154 x IVP101 (n=24)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
clean_C93154_dihaploids_IVP101 <- parent_snps(snps2, clean_C93154_dihaploids_GT, clean_C93154_dihaploids_DP, IVP101_GT, IVP101_DP) %>%
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT),
         clean_C93154_dihaploids = ifelse(IVP101_GT == "0/0", ALT, REF)) %>%
  select(Chrom = CHROM, Pos = POS, Ref = REF, clean_C93154_dihaploids, IVP101)
write_tsv(clean_C93154_dihaploids_IVP101, "clean_C93154_dihaploids-IVP101-SNP.tsv", col_names = T)

#' 
#' C93.154 x IVP35 (n=88)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
clean_C93154_dihaploids_IVP35 <- parent_snps(snps2, clean_C93154_dihaploids_GT, clean_C93154_dihaploids_DP, IVP35_GT, IVP35_DP) %>%
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT),
         clean_C93154_dihaploids = ifelse(IVP35_GT == "0/0", ALT, REF)) %>%
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP35, clean_C93154_dihaploids)
write_tsv(clean_C93154_dihaploids_IVP35, "clean_C93154_dihaploids-IVP35-SNP.tsv", col_names = T)

#' 
#' C93.154 x PL4 (n=161)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
clean_C93154_dihaploids_PL4 <- parent_snps(snps2, clean_C93154_dihaploids_GT, clean_C93154_dihaploids_DP, PL4_GT, PL4_DP) %>%
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT),
         clean_C93154_dihaploids = ifelse(PL4_GT == "0/0", ALT, REF)) %>%
  select(Chrom = CHROM, Pos = POS, Ref = REF, PL4, clean_C93154_dihaploids)
write_tsv(clean_C93154_dihaploids_PL4, "clean_C93154_dihaploids-PL4-SNP.tsv", col_names = T)

#' 
#' LR00.022 x IVP101 (n=2)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
clean_LR00022_dihaploids_IVP101 <- parent_snps(snps2, clean_LR00022_dihaploids_GT, clean_LR00022_dihaploids_DP, IVP101_GT, IVP101_DP) %>%
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT),
         clean_LR00022_dihaploids = ifelse(IVP101_GT == "0/0", ALT, REF)) %>%
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP101, clean_LR00022_dihaploids)
write_tsv(clean_LR00022_dihaploids_IVP101, "clean_LR00022_dihaploids-IVP101-SNP.tsv", col_names = T)

#' 
#' LR00.022 x IVP35 (n=2)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
clean_LR00022_dihaploids_IVP35 <- parent_snps(snps2, clean_LR00022_dihaploids_GT, clean_LR00022_dihaploids_DP, IVP35_GT, IVP35_DP) %>%
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT),
         clean_LR00022_dihaploids = ifelse(IVP35_GT == "0/0", ALT, REF)) %>%
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP35, clean_LR00022_dihaploids)
write_tsv(clean_LR00022_dihaploids_IVP35, "clean_LR00022_dihaploids-IVP35-SNP.tsv", col_names = T)

#' 
#' LR00.022 x PL4 (n=59)
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
clean_LR00022_dihaploids_PL4 <- parent_snps(snps2, clean_LR00022_dihaploids_GT, clean_LR00022_dihaploids_DP, PL4_GT, PL4_DP) %>%
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT),
         clean_LR00022_dihaploids = ifelse(PL4_GT == "0/0", ALT, REF)) %>%
  select(Chrom = CHROM, Pos = POS, Ref = REF, PL4, clean_LR00022_dihaploids)
write_tsv(clean_LR00022_dihaploids_PL4, "clean_LR00022_dihaploids-PL4-SNP.tsv", col_names = T)

#' 
## ------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::purl("MM_parent_snps.Rmd", documentation = 2)

