From results of ```1_parent_snps``` and ```2_dosage``` determine whether using
pooled euploid dihaploids from WA.077 yields similar results as using WA.077 itself.

Specify dihaploids extracted from WA.077 that did not have a trisomy according to 
```2_dosage``` as "clean" dihaploids. Then, process the alignments of these dihaploids
for SNP calling and genotyping. Output is a SNP list that can be used for binned
genotyping as described in ```2_dosage```.

Compare output from both SNP lists. Specifically, I will check the genome-wide and
per-chromosome HI SNP % and check whether the chr08 breakpoints in MM247 are concordant.

