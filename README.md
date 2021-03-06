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
4. ```4_dihaploid_pools```: Variant calling including pools of euploid dihaploids as a
   stand-in for tetraploid parents that we didn't sequence.
5. ```5_dosage_with_pools```: Chromosome parental origin for dihaploids and hybrids for
   which maternal SNPs were determined from dihaploid pools.
6. ```6_offchrom```: Call variants and look for segmental secondary introgression events
   in HI addition dihaploids.
7. ```7_haplotypes```: For tetraploid hybrids, genotype the disome inherited by the HI.

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
snakemake --cores <cores> --use-conda
```

4. Short read QC, alignment and variant calling for tetraploid parents and haploid inducers

```
cd ../1_parent_snps
snakemake --cores <cores> --use-conda
```

5. Short read QC and alignment for 1,001 dihaploids, 134 hybrids and 14 tetraploid selfs

```
cd ../2_dosage
snakemake --cores <cores> --use-conda
```

6. Evaluate pooled low coverage sequencing for HI addition lines

```
cd ../3_dihaploid_sub_test
snakemake --cores <cores> --use-conda
```

7. Variant calling including dihaploid pools

```
cd ../4_dihaploid_pools
snakemake --cores <cores> --use-conda
```

8. Chromosome dosage including dihaploid pools

```
cd ../5_dosage_with_pools
snakemake --cores <cores> --use-conda
```

9. Short read QC, alignment and variant calling for HI addition lines. Add in pooled 
   euploid dihaploids as a substitute for the 4x parent for variant calling.

```
cd ../6_offchrom
snakemake --cores <cores> --use-conda
```

10. Haplotype extraction from addition lines and analyses in hybrids.

```
cd ../7_haplotypes
snakemake --cores <cores> --use-conda
```