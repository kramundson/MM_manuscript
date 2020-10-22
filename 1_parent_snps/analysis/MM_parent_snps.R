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
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)

#' 
#' ## Functions:
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# separate sample-specific VCF attributes
sep <- function(...) {
  dots <- list(...)
  separate_(..., into = paste(dots[[2]], attributes[[1]], sep = "_"), convert = T, sep = ":")
}

## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# either download files to local or mount to server via, e.g., sshfs
files <- dir(path = "~/Desktop/mount/share/comailab/kramundson/MM_manuscript/1_parent_snps/data/calls/-filtered-",
             full.names = T)

# temp working version
# files <- dir(pattern = "-filtered-",
#              path = "../data/calls/",
#              full.names = T)
# files

#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
snps <- map_dfr(files, function(x) read_tsv(x, col_names = T, na = "NA"))

#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
names(table(snps$FORMAT)) # should only have one entry. does.
attributes <- str_split(names(table(snps$FORMAT[1])), ":", simplify = F)
attributes[[1]]
sample_vars <- colnames(snps)[-c(seq(1,49), ncol(snps))]

#' 
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
snps2 <- snps %>%
  Reduce(f = sep, x = sample_vars)

#' 
#' WA.077 x IVP101 (n=50)
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(MM302_GT, IVP101_GT, snps2, 10) %>% 
  mutate(MM302 = ifelse(MM302_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, IVP101, MM302) %>% 
  filter(CHROM %in% sprintf("chr%0.2d", 1:12)) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., paste0(Sys.Date(), "-IVP101-MM302-hom-SNP.tsv"), col_names = T)

#' 
#' WA.077 x IVP35 (n=134)
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(MM302_GT, IVP35_GT, snps2, 10) %>% 
  mutate(MM302 = ifelse(MM302_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, IVP35, MM302) %>% 
  filter(CHROM %in% sprintf("chr%0.2d", 1:12)) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., paste0(Sys.Date(), "-IVP35-MM302-hom-SNP.tsv"), col_names = T)

#' 
#' WA.077 x PL4 (n=107)
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(MM302_GT, PL4_GT, snps2, 10) %>% 
  mutate(MM302 = ifelse(MM302_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, PL4, MM302) %>% 
  filter(CHROM %in% sprintf("chr%0.2d", 1:12)) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., paste0(Sys.Date(), "-PL4-MM302-hom-SNP.tsv"), col_names = T)

#' 
#' LR00.014 x IVP101 (n=30)
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(MM293_GT, IVP101_GT, snps2, 10) %>% 
  filter(CHROM %in% sprintf("chr%0.2d", 1:12)) %>%
  mutate(MM293 = ifelse(MM293_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT)) %>% 
  dplyr::select(CHROM, POS, REF, IVP101, MM293) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(paste0(Sys.Date(), "-IVP101-MM293-hom-SNP.tsv"), col_names = T)

#' 
#' LR00.014 x IVP35 (n=77)
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(MM293_GT, IVP35_GT, snps2, 10) %>% 
  mutate(MM293 = ifelse(MM293_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, IVP35, MM293) %>% 
  filter(CHROM %in% sprintf("chr%0.2d", 1:12)) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., paste0(Sys.Date(), "-IVP35-MM293-hom-SNP.tsv"), col_names = T)

#' 
#' LR00.014 x PL4 (n=67)
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(MM293_GT, PL4_GT, snps2, 10) %>% 
  mutate(MM293 = ifelse(MM293_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(PL4 = ifelse(IVP35_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, IVP35, MM293) %>% 
  filter(CHROM %in% sprintf("chr%0.2d", 1:12)) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., paste0(Sys.Date(), "-PL4-MM293-hom-SNP.tsv"), col_names = T)

#' 
#' LR00.026 x IVP101 (n=4)
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(MM294_GT, IVP101_GT, snps2, 10) %>% 
  mutate(MM294 = ifelse(MM294_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, IVP101, MM294) %>% 
  filter(CHROM %in% sprintf("chr%0.2d", 1:12)) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., paste0(Sys.Date(), "-IVP101-MM294-hom-SNP.tsv"), col_names = T) 

#' 
#' LR00.026 x IVP35 (n=36)
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(MM294_GT, IVP35_GT, snps2, 10) %>% 
  mutate(MM294 = ifelse(MM294_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, IVP35, MM294) %>% 
  filter(CHROM %in% sprintf("chr%0.2d", 1:12)) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., paste0(Sys.Date(), "-IVP35-MM294-hom-SNP.tsv"), col_names = T)

#' 
#' LR00.026 x PL4 (n=35)
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(MM294_GT, PL4_GT, snps2, 10) %>% 
  mutate(MM294 = ifelse(MM294_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, PL4, MM294) %>% 
  filter(CHROM %in% sprintf("chr%0.2d", 1:12)) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., paste0(Sys.Date(), "-PL4-MM294-hom-SNP.tsv"), col_names = T)

#' 
#' Atlantic x IVP35 (n=5)
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(Atlantic_GT, IVP35_GT, snps2, 10) %>% 
  mutate(Atlantic = ifelse(Atlantic_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, IVP35, Atlantic) %>% 
  filter(CHROM %in% sprintf("chr%0.2d", 1:12)) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., paste0(Sys.Date(), "-IVP35-Atlantic-hom-SNP.tsv"), col_names = T) 

#' 
#' Atlantic x PL4 (n=10)
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(Atlantic_GT, PL4_GT, snps2, 10) %>% 
  mutate(Atlantic = ifelse(Atlantic_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, PL4, Atlantic) %>% 
  filter(CHROM %in% sprintf("chr%0.2d", 1:12)) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., paste0(Sys.Date(), "-PL4-Atlantic-hom-SNP.tsv"), col_names = T)

#' 
#' Desiree x IVP101 (n=2)
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(Desiree_GT, IVP101_GT, snps2, 10) %>% 
  mutate(Desiree = ifelse(Desiree_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, IVP101, Desiree) %>% 
  filter(CHROM %in% sprintf("chr%0.2d", 1:12)) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>%
  write_tsv(., paste0(Sys.Date(), "-IVP101-Desiree-hom-SNP.tsv"), col_names = T)

#' 
#' Desiree x IVP35 (n=2)
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(Desiree_GT, IVP35_GT, snps2, 10) %>% 
  mutate(Desiree = ifelse(Desiree_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, IVP35, Desiree) %>% 
  filter(CHROM %in% sprintf("chr%0.2d", 1:12)) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., paste0(Sys.Date(), "-IVP35-Desiree-hom-SNP.tsv"), col_names = T)

#' 
#' Desiree x PL4 (n=6)
## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
filter_homozygous_vars(Desiree_GT, PL4_GT, snps2, 10) %>% 
  mutate(Desiree = ifelse(Desiree_GT == "0/0/0/0", REF, ALT)) %>% 
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT)) %>% 
  select(CHROM, POS, REF, PL4, Desiree) %>% 
  filter(CHROM %in% sprintf("chr%0.2d", 1:12)) %>% 
  rename(Chrom = CHROM, Pos = POS, Ref = REF) %>% 
  write_tsv(., paste0(Sys.Date(), "-PL4-Desiree-hom-SNP.tsv"), col_names = T)

