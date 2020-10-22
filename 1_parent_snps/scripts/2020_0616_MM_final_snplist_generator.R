library(tidyverse)
library(stringr)

sep <- function(...) {
  dots <- list(...)
  separate_(..., into = paste(dots[[2]], attributes[[1]], sep = "_"), convert = T, sep = ":")
}

parent_snps <- function(df, tetgt, tetdp, higt, hidp, tetmin=10, himin=10) {
  out <- df %>% 
    filter({{tetdp}} >= tetmin,
           {{tetgt}} %in% c("0/0/0/0", "1/1/1/1"),
           {{hidp}} >= himin,
           {{higt}} %in% c("0/0", "1/1"),
           !({{higt}} == "0/0" & {{tetgt}} == "0/0/0/0"),
           !({{higt}} == "1/1" & {{tetgt}} == "1/1/1/1"))
  return(out)
}

files <- dir(pattern = "-MM-filtered-")
print(files)

snps <- files %>% 
  map(read_tsv, col_names = T, na = "NA") %>% 
  bind_rows()

names(table(snps$FORMAT)) # should only have one entry. does.
attributes <- str_split(names(table(snps$FORMAT[1])), ":", simplify = F)
attributes[[1]]
sample_vars <- colnames(snps)[-c(seq(1,49), ncol(snps))] # last one is MQ.diff, which we also do not want. Will want to have a more robust approach than indexing in the future.
print(sample_vars)
site_vars <- colnames(snps)[c(1:49,ncol(snps))]
print(site_vars)

snps2 <- snps %>% 
  Reduce(f = sep, x = sample_vars)

# clean C01020 dihaploids
# c01020_IVP101 <- parent_snps(snps2, clean_C01020_dihaploids_GT, clean_C01020_dihaploids_DP, IVP101_GT, IVP101_DP) %>% 
#   mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT),
#          clean_C01020_dihaploids = ifelse(IVP101_GT == "0/0", ALT, REF)) %>% 
#   select(Chrom = CHROM, Pos = POS, Ref = REF, clean_C01020_dihaploids, IVP101)
# write_tsv(c01020_IVP101, paste0(Sys.Date(),"-clean_C01020_dihaploids-IVP101-SNP.tsv"), col_names = T)
# 
# c01020_PL4 <- parent_snps(snps2, clean_C01020_dihaploids_GT, clean_C01020_dihaploids_DP, PL4_GT, PL4_DP) %>% 
#   mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT),
#          clean_C01020_dihaploids = ifelse(PL4_GT == "0/0", ALT, REF)) %>% 
#   select(Chrom = CHROM, Pos = POS, Ref = REF, clean_C01020_dihaploids, PL4)
# write_tsv(c01020_PL4, paste0(Sys.Date(),"-clean_C01020_dihaploids-PL4-SNP.tsv"), col_names = T)

# from snps2, make lists of parental snps for dosage analyses
# clean C01020 dihaploids
c01020_IVP101 <- parent_snps(snps2, clean_C01020_dihaploids_GT, clean_C01020_dihaploids_DP, IVP101_GT, IVP101_DP) %>% 
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT),
         clean_C01020_dihaploids = ifelse(IVP101_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP101, clean_C01020_dihaploids)
write_tsv(c01020_IVP101, paste0(Sys.Date(),"-clean_C01020_dihaploids-IVP101-SNP.tsv"), col_names = T)

c01020_PL4 <- parent_snps(snps2, clean_C01020_dihaploids_GT, clean_C01020_dihaploids_DP, PL4_GT, PL4_DP) %>% 
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT),
         clean_C01020_dihaploids = ifelse(PL4_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, PL4, clean_C01020_dihaploids)
write_tsv(c01020_PL4, paste0(Sys.Date(),"-clean_C01020_dihaploids-PL4-SNP.tsv"), col_names = T)

# LR00014 dihaploids (MM293)
MM293_IVP101 <- parent_snps(snps2, MM293_GT, MM293_DP, IVP101_GT, IVP101_DP) %>% 
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT),
         MM293 = ifelse(IVP101_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP101, MM293)
write_tsv(MM293_IVP101, paste0(Sys.Date(),"-MM293-IVP101-SNP.tsv"), col_names = T)

MM293_IVP35 <- parent_snps(snps2, MM293_GT, MM293_DP, IVP35_GT, IVP35_DP) %>% 
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT),
         MM293 = ifelse(IVP35_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP35, MM293)
write_tsv(MM293_IVP35, paste0(Sys.Date(),"-MM293-IVP35-SNP.tsv"), col_names = T)

MM293_PL4 <- parent_snps(snps2, MM293_GT, MM293_DP, PL4_GT, PL4_DP) %>% 
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT),
         MM293 = ifelse(PL4_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, PL4, MM293)
write_tsv(MM293_PL4, paste0(Sys.Date(),"-MM293-PL4-SNP.tsv"), col_names = T)

# LR00026 dihaploids (MM294)
MM294_IVP101 <- parent_snps(snps2, MM294_GT, MM294_DP, IVP101_GT, IVP101_DP) %>% 
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT),
         MM294 = ifelse(IVP101_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP101, MM294)
write_tsv(MM294_IVP101, paste0(Sys.Date(),"-MM294-IVP101-SNP.tsv"), col_names = T)

MM294_IVP35 <- parent_snps(snps2, MM294_GT, MM294_DP, IVP35_GT, IVP35_DP) %>% 
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT),
         MM294 = ifelse(IVP35_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP35, MM294)
write_tsv(MM294_IVP35, paste0(Sys.Date(),"-MM294-IVP35-SNP.tsv"), col_names = T)

MM294_PL4 <- parent_snps(snps2, MM294_GT, MM294_DP, PL4_GT, PL4_DP) %>% 
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT),
         MM294 = ifelse(PL4_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, PL4, MM294)
write_tsv(MM294_PL4, paste0(Sys.Date(),"-MM294-PL4-SNP.tsv"), col_names = T)

# WA077 (MM302)
MM302_IVP101 <- parent_snps(snps2, MM302_GT, MM302_DP, IVP101_GT, IVP101_DP) %>% 
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT),
         MM302 = ifelse(IVP101_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP101, MM302)
write_tsv(MM302_IVP101, paste0(Sys.Date(),"-MM302-IVP101-SNP.tsv"), col_names = T)

MM302_IVP35 <- parent_snps(snps2, MM302_GT, MM302_DP, IVP35_GT, IVP35_DP) %>% 
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT),
         MM302 = ifelse(IVP35_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP35, MM302)
write_tsv(MM302_IVP35, paste0(Sys.Date(),"-MM302-IVP35-SNP.tsv"), col_names = T)

MM302_PL4 <- parent_snps(snps2, MM302_GT, MM302_DP, PL4_GT, PL4_DP) %>% 
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT),
         MM302 = ifelse(PL4_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, PL4, MM302)
write_tsv(MM302_PL4, paste0(Sys.Date(),"-MM302-PL4-SNP.tsv"), col_names = T)

# LR00022 (dihaploid pool)
clean_LR00022_dihaploids_IVP101 <- parent_snps(snps2, clean_LR00022_dihaploids_GT, clean_LR00022_dihaploids_DP, IVP101_GT, IVP101_DP) %>% 
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT),
         clean_LR00022_dihaploids = ifelse(IVP101_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP101, clean_LR00022_dihaploids)
write_tsv(clean_LR00022_dihaploids_IVP101, paste0(Sys.Date(),"-clean_LR00022_dihaploids-IVP101-SNP.tsv"), col_names = T)

clean_LR00022_dihaploids_IVP35 <- parent_snps(snps2, clean_LR00022_dihaploids_GT, clean_LR00022_dihaploids_DP, IVP35_GT, IVP35_DP) %>% 
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT),
         clean_LR00022_dihaploids = ifelse(IVP35_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP35, clean_LR00022_dihaploids)
write_tsv(clean_LR00022_dihaploids_IVP35, paste0(Sys.Date(),"-clean_LR00022_dihaploids-IVP35-SNP.tsv"), col_names = T)

clean_LR00022_dihaploids_PL4 <- parent_snps(snps2, clean_LR00022_dihaploids_GT, clean_LR00022_dihaploids_DP, PL4_GT, PL4_DP) %>% 
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT),
         clean_LR00022_dihaploids = ifelse(PL4_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, PL4, clean_LR00022_dihaploids)
write_tsv(clean_LR00022_dihaploids_PL4, paste0(Sys.Date(),"-clean_LR00022_dihaploids-PL4-SNP.tsv"), col_names = T)

# C93154 dihaploid pool
clean_C93154_dihaploids_IVP101 <- parent_snps(snps2, clean_C93154_dihaploids_GT, clean_C93154_dihaploids_DP, IVP101_GT, IVP101_DP) %>% 
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT),
         clean_C93154_dihaploids = ifelse(IVP101_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, clean_C93154_dihaploids, IVP101)
write_tsv(clean_C93154_dihaploids_IVP101, paste0(Sys.Date(),"-clean_C93154_dihaploids-IVP101-SNP.tsv"), col_names = T)

clean_C93154_dihaploids_IVP35 <- parent_snps(snps2, clean_C93154_dihaploids_GT, clean_C93154_dihaploids_DP, IVP35_GT, IVP35_DP) %>% 
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT),
         clean_C93154_dihaploids = ifelse(IVP35_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP35, clean_C93154_dihaploids)
write_tsv(clean_C93154_dihaploids_IVP35, paste0(Sys.Date(),"-clean_C93154_dihaploids-IVP35-SNP.tsv"), col_names = T)

clean_C93154_dihaploids_PL4 <- parent_snps(snps2, clean_C93154_dihaploids_GT, clean_C93154_dihaploids_DP, PL4_GT, PL4_DP) %>% 
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT),
         clean_C93154_dihaploids = ifelse(PL4_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, PL4, clean_C93154_dihaploids)
write_tsv(clean_C93154_dihaploids_PL4, paste0(Sys.Date(),"-clean_C93154_dihaploids-PL4-SNP.tsv"), col_names = T)

# C91640 dihaploid pool
clean_C91640_dihaploids_IVP101 <- parent_snps(snps2, clean_C91640_dihaploids_GT, clean_C91640_dihaploids_DP, IVP101_GT, IVP101_DP) %>% 
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT),
         clean_C91640_dihaploids = ifelse(IVP101_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP101, clean_C91640_dihaploids)
write_tsv(clean_C91640_dihaploids_IVP101, paste0(Sys.Date(),"-clean_C91640_dihaploids-IVP101-SNP.tsv"), col_names = T)

clean_C91640_dihaploids_IVP35 <- parent_snps(snps2, clean_C91640_dihaploids_GT, clean_C91640_dihaploids_DP, IVP35_GT, IVP35_DP) %>% 
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT),
         clean_C91640_dihaploids = ifelse(IVP35_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP35, clean_C91640_dihaploids)
write_tsv(clean_C91640_dihaploids_IVP35, paste0(Sys.Date(),"-clean_C91640_dihaploids-IVP35-SNP.tsv"), col_names = T)

clean_C91640_dihaploids_PL4 <- parent_snps(snps2, clean_C91640_dihaploids_GT, clean_C91640_dihaploids_DP, PL4_GT, PL4_DP) %>% 
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT),
         clean_C91640_dihaploids = ifelse(PL4_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, PL4, clean_C91640_dihaploids)
write_tsv(clean_C91640_dihaploids_PL4, paste0(Sys.Date(),"-clean_C91640_dihaploids-PL4-SNP.tsv"), col_names = T)

# CIP Desiree
CIP_Desiree_IVP101 <- parent_snps(snps2, CIP_Desiree_GT, CIP_Desiree_DP, IVP101_GT, IVP101_DP) %>% 
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT),
         CIP_Desiree = ifelse(IVP101_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP101, CIP_Desiree)
write_tsv(CIP_Desiree_IVP101, paste0(Sys.Date(),"-CIP_Desiree-IVP101-SNP.tsv"), col_names = T)

CIP_Desiree_IVP35 <- parent_snps(snps2, CIP_Desiree_GT, CIP_Desiree_DP, IVP35_GT, IVP35_DP) %>% 
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT),
         CIP_Desiree = ifelse(IVP35_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP35, CIP_Desiree)
write_tsv(CIP_Desiree_IVP35, paste0(Sys.Date(),"-CIP_Desiree-IVP35-SNP.tsv"), col_names = T)

CIP_Desiree_PL4 <- parent_snps(snps2, CIP_Desiree_GT, CIP_Desiree_DP, PL4_GT, PL4_DP) %>% 
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT),
         CIP_Desiree = ifelse(PL4_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, PL4, CIP_Desiree)
write_tsv(CIP_Desiree_PL4, paste0(Sys.Date(),"-CIP_Desiree-PL4-SNP.tsv"), col_names = T)

# Atlantic
Atlantic_IVP101 <- parent_snps(snps2, Atlantic_GT, Atlantic_DP, IVP101_GT, IVP101_DP) %>% 
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT),
         Atlantic = ifelse(IVP101_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP101, Atlantic)
write_tsv(Atlantic_IVP101, paste0(Sys.Date(),"-Atlantic-IVP101-SNP.tsv"), col_names = T)

Atlantic_IVP35 <- parent_snps(snps2, Atlantic_GT, Atlantic_DP, IVP35_GT, IVP35_DP) %>% 
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT),
         Atlantic = ifelse(IVP35_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP35, Atlantic)
write_tsv(Atlantic_IVP35, paste0(Sys.Date(),"-Atlantic-IVP35-SNP.tsv"), col_names = T)

Atlantic_PL4 <- parent_snps(snps2, Atlantic_GT, Atlantic_DP, PL4_GT, PL4_DP) %>% 
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT),
         Atlantic = ifelse(PL4_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, PL4, Atlantic)
write_tsv(Atlantic_PL4, paste0(Sys.Date(),"-Atlantic-PL4-SNP.tsv"), col_names = T)

# clean_93003_dihaploids
clean_93003_dihaploids_IVP101 <- parent_snps(snps2, clean_93003_dihaploids_GT, clean_93003_dihaploids_DP, IVP101_GT, IVP101_DP) %>% 
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT),
         clean_93003_dihaploids = ifelse(IVP101_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP101, clean_93003_dihaploids)
write_tsv(clean_93003_dihaploids_IVP101, paste0(Sys.Date(),"-clean_93003_dihaploids-IVP101-SNP.tsv"), col_names = T)

clean_93003_dihaploids_IVP35 <- parent_snps(snps2, clean_93003_dihaploids_GT, clean_93003_dihaploids_DP, IVP35_GT, IVP35_DP) %>% 
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT),
         clean_93003_dihaploids = ifelse(IVP35_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP35, clean_93003_dihaploids)
write_tsv(clean_93003_dihaploids_IVP35, paste0(Sys.Date(),"-clean_93003_dihaploids-IVP35-SNP.tsv"), col_names = T)

clean_93003_dihaploids_PL4 <- parent_snps(snps2, clean_93003_dihaploids_GT, clean_93003_dihaploids_DP, PL4_GT, PL4_DP) %>% 
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT),
         clean_93003_dihaploids = ifelse(PL4_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, PL4, clean_93003_dihaploids)
write_tsv(clean_93003_dihaploids_PL4, paste0(Sys.Date(),"-clean_93003_dihaploids-PL4-SNP.tsv"), col_names = T)

# clean_LR00022_dihaploids
clean_LR00022_dihaploids_IVP101 <- parent_snps(snps2, clean_LR00022_dihaploids_GT, clean_LR00022_dihaploids_DP, IVP101_GT, IVP101_DP) %>% 
  mutate(IVP101 = ifelse(IVP101_GT == "0/0", REF, ALT),
         clean_LR00022_dihaploids = ifelse(IVP101_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP101, clean_LR00022_dihaploids)
write_tsv(clean_LR00022_dihaploids_IVP101, paste0(Sys.Date(),"-clean_LR00022_dihaploids-IVP101-SNP.tsv"), col_names = T)

clean_LR00022_dihaploids_IVP35 <- parent_snps(snps2, clean_LR00022_dihaploids_GT, clean_LR00022_dihaploids_DP, IVP35_GT, IVP35_DP) %>% 
  mutate(IVP35 = ifelse(IVP35_GT == "0/0", REF, ALT),
         clean_LR00022_dihaploids = ifelse(IVP35_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, IVP35, clean_LR00022_dihaploids)
write_tsv(clean_LR00022_dihaploids_IVP35, paste0(Sys.Date(),"-clean_LR00022_dihaploids-IVP35-SNP.tsv"), col_names = T)

clean_LR00022_dihaploids_PL4 <- parent_snps(snps2, clean_LR00022_dihaploids_GT, clean_LR00022_dihaploids_DP, PL4_GT, PL4_DP) %>% 
  mutate(PL4 = ifelse(PL4_GT == "0/0", REF, ALT),
         clean_LR00022_dihaploids = ifelse(PL4_GT == "0/0", ALT, REF)) %>% 
  select(Chrom = CHROM, Pos = POS, Ref = REF, PL4, clean_LR00022_dihaploids)
write_tsv(clean_LR00022_dihaploids_PL4, paste0(Sys.Date(),"-clean_LR00022_dihaploids-PL4-SNP.tsv"), col_names = T)