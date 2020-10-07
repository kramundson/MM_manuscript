From results of ```1_dosage``` and ```3_chromosome_ori``` determine whether using
pooled euploid dihaploids from WA.077 yields similar results as using WA.077 itself.

There are a couple factors at play:
If dihaploids have smaller introgressions, loci will not be included in the analysis. The
number of markers will decrease.

Read type. WA.077 was sequenced 2x 150, all dihaploids were sequenced single-end 100.

From ```1_dosage```, keep the "clean" dihaploids of WA.077, i.e., those that don't have
a trisomy. Merge reads of these dihaploids and substitute WA.077 for this pool during
the variant calling step.

