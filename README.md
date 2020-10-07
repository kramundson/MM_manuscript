Pipelines and analyses of MM manuscript

The project is composed of six main analyses, each of which is written as its own
Snakemake workflow. The analyses are meant to be run in order, i.e., the output of earlier
analyses are required input for later analyses.

A description of the workflows, in the order they should be run:

1. ```1_parent_snps```: Variant calling to identify SNPs used for inferring parental
   origin of whole-chromosome and segmental aneuploidy.
2. ```2_dosage```: Analysis of low pass sequencing for 1,001 dihaploids and 134 hybrids to
   identify whole-chromosome and segmental aneuploids, then determine parental origin.
3. ```3_dihaploid_sub_test```: Proof of concept in light of results from ```2_dosage```
   Here, I evaluated the use of pooling dihaploids as a stand-in for the tetraploid parent.
4. ```4_offchrom```: Call variants and look for segmental secondary introgression events
   in HI addition dihaploids.
5. ```5_haplotypes```: For tetraploid hybrids, genotype the disome inherited by the HI.

To run these analyses and reproduce the figures and tables of the manuscript:

1. Clone this repo

```
git clone https://github.com/kramundson/MM_manuscript
cd MM_manuscript
```

2. Install dependencies from conda environment. Need also to install GATK3 from source.

```
# todo debug conda env conflicts
conda env create -n MM -f environment.yaml
# todo add instructions for gatk legacy installation
```

3. Build reference genome (DM1-3 v4.04) used for all pipelines.

```
cd ref
snakemake
```

4. Short read QC, alignment and variant calling for tetraploid parents and haploid inducers

```
cd ../1_parent_snps
snakemake
```

5. Short read QC and alignment for 1,001 dihaploids, 134 hybrids and 14 tetraploid selfs

```
cd ../2_dosage
snakemake
```

6. Evaluate pooled low coverage sequencing for HI addition lines

```
cd ../3_dihaploid_sub_test
snakemake
```

7. Short read QC, alignment and variant calling for HI addition lines. Add in pooled 
   euploid dihaploids as a substitute for the 4x parent for variant calling.

```
cd ../4_offchrom
snakemake
```

8. Haplotype extraction from addition lines and analyses in hybrids.

```
cd ../5_haplotypes
snakemake
```