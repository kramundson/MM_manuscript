#!/share/comailab/kramundson/miniconda3/bin/Rscript
# 2019_0225
# Kirk Amundson
# 2019_0225_filter_cheat_MM.R

# USAGE: Rscript 2019_0225_filter_cheat_MM.R <input vcf>

library(tidyverse)
library(stringr)

# generic_qual_filter <- function(file, mom, dad) {
generic_qual_filter <- function(file) {
  # parse header, requires zgrep installed on system
  print("Parsing VCF header")
  vcf_header <- system(paste("grep '#C'", file), intern = T) %>%
    str_replace("^#C", "C") %>%
    str_replace_all("[0-9]x_", "") %>%
    str_split(pattern = "\t")
  
  # read in file
  print("Reading in VCF")
  vcf <- read_tsv(file, col_names = vcf_header[[1]], comment = "#", na = c("NA", ".", "./.", "././.", "./././."))
  
  # parse INFO column; this is the per-locus attributes
  # this assumes that the INFO column format of the first row is valid for all rows.
  # For the MM cheat dataset, I validated this at the command line:
  # zgrep -v "^#" all-calls.vcf.gz | cut -f 8 | sed -e 's/=[A-Za-z0-9.,-]\+//g' | uniq -c
  # returns 26135465 AB;ABP;AC;AF;AN;AO;CIGAR;DP;DPB;DPRA;EPP;EPPR;GTI;LEN;MEANALT;MQM;MQMR;NS;NUMALT;ODDS;PAIRED;PAIREDR;PAO;PQA;PQR;PRO;QA;QR;RO;RPL;RPP;RPPR;RPR;RUN;SAF;SAP;SAR;SRF;SRP;SRR;TYPE
  info <- str_split(vcf$INFO[1], ";")[[1]] %>% 
    str_replace("=.+", "")
  print(info)
  
  # parse FORMAT column, this is the per-sample attributes
  print("Parsing FORMAT")
  attributes <- str_split(names(table(vcf$FORMAT)), ":", simplify = F)
  print(attributes)[[1]]
  
  # apply filters, opening up per-sample columns as necessary
  print("Applying site quality filters")
  vcf_filter <- vcf %>% 
    mutate(INFO = str_replace_all(INFO, "[A-Za-z]*=", "")) %>%
    separate(INFO, into = info, sep = ";", convert = T) %>% 
    filter(FORMAT != "GT:GQ:DP:AD:RO:QR:AO:QA:GL") %>% 
    filter(NUMALT == 1) %>% # added this filter step on 4 May 2020. Output was entirely NUMALT==1 anyway, but better to enforce it
    filter(CIGAR == "1X") %>% 
    filter(QUAL >= 20) %>% 
    mutate(MQM = as.numeric(MQM)) %>% # for example, MQM is read in as a character by default
    filter(MQM >= 50) %>% 
    filter(MQMR >= 50) %>% # MQMR usually gets imputed correctly, either as double or numeric
    mutate(MQ.diff = abs(MQM - MQMR)) %>%
    filter(MQ.diff < 10) %>% 
    mutate(RPPR = as.numeric(RPPR)) %>% # obscure, would like to remove
    filter(RPPR <= 20) %>% # obscure, would like to remove
    mutate(RPP = as.numeric(RPP)) %>% 
    filter(RPP <= 20) %>% 
    mutate(EPP = as.numeric(EPP)) %>% 
    filter(EPP <= 20) %>% 
    mutate(EPPR = as.numeric(EPPR)) %>% 
    filter(EPPR <= 20) %>% 
    mutate(SAP = as.numeric(SAP)) %>% 
    filter(SAP <= 20) %>% 
    mutate(SRP = as.numeric(SRP)) %>% 
    filter(SRP <= 20)
  
  # write out filtered chromosome variant as tsv, then when all are done, read back in, open up sample-specific attributes, and filter on those.
  # want to generate a list of parental SNPs for each family comparison
  writeout <- gsub("calls", "filtered-calls", file)
  print(paste("Writing out filtered calls to:", writeout))
  write.table(arrange(vcf_filter, POS), gzfile(writeout),
              quote=F, sep = '\t', eol = '\n',
              na = "NA", row.names = F, col.names = T)
}

args = commandArgs(trailingOnly = TRUE)
file <- args[1]
generic_qual_filter(file)
