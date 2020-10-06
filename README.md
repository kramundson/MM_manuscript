Pipelines and analyses of MM manuscript

The project is composed of six main analyses, each of which is written as its own
Snakemake workflow. The analyses are meant to be run in order, i.e., the output of earlier
analyses are required input for later analyses.

A description of the workflows, in the order they should be run:

1. ```1_dosage```: Low pass sequencing of 1,001 dihaploids and 134 hybrids to identify
   whole-chromosome and segmental aneuploids.
2. ```2_parent_snps```: Variant calling to identify SNPs used for inferring parental
   origin of whole-chromosome and segmental aneuploidy.
3. ```3_chromosome_ori```: Scripts for evaluating parental origin of whole-chromosome and
   segmental aneuploidy.
4. ```4_cheat_test```: Proof of concept in light of results from ```3_chromosome_ori```
   Here, I evaluated the use of pooling dihaploids as a stand-in for the tetraploid parent.
5. ```5_offchrom```: Call variants and look for segmental secondary introgression events
   in HI addition dihaploids.
6. ```6_haplotypes```: 




1. Clone this repo

```
git clone https://github.com/kramundson/MM_manuscript
cd MM_manuscript
```

2. Build conda environment. Need to manually install legacy GATK.

```
conda env create -n MM -f environment.yaml

# todo add instructions for gatk legacy installation
```

3. Build reference genome (DM1-3 v4.04) used for all pipelines.

```
cd ref
snakemake
```

4. Short read QC and alignment for 1,001 dihaploids, 134 hybrids and 14 tetraploid selfs

```
cd ../1_dosage
snakemake
```

5. Short read QC, alignment and variant calling for tetraploid parents and haploid inducers

```
cd ../2_parent_snps
snakemake
```

6. Trisomy parental origin analysis

```
cd ../3_chromosome_ori
snakemake
```

7. Evaluate pooled low coverage sequencing for HI addition lines

```
cd ../4_cheat_test
snakemake
```

8. Short read QC, alignment and variant calling for HI addition lines. Add in pooled 
   euploid dihaploids as a substitute for the 4x parent for variant calling.

```
cd ../5_offchrom
snakemake
```

8. Haplotype extraction from addition lines and analyses in hybrids.

```
cd ../6_haplotypes
snakemake
```