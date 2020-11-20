# RidgelinesCWRdist

## Aim

The R code in this repository infers K2P genetic distances between CWRs and crop species (using aligned DNA matrices in FASTA format) and displays these interspecific distances using ridgeline plots.

## Steps

1. Infer Kimura's 2-parameters distances.
2. Convert pairwise distance matrix and add species IDs.
3. Tidy matrix to infer interspecific distances between CWRs and crop.
4. Order species (or CWRs) for plot (based on genetic distances, here by comparing minimum distances per species).
5. Draw ridgeline plots.
6. Export plot in pdf format.