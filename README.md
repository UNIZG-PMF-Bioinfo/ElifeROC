# ElifeROC
ROC analysis for Samples+-BIX

## 1_PreparingIS_and_EpigenTables.R
First, cut the genome into 1kb bins and calculate average score for epigenetic features on this regions. Add the values to integration sites +- bix.
https://github.com/UNIZG-PMF-Bioinfo/ElifeROC/blob/master/1_PreparingIS_and_EpigenTables.R

## 2_Generation_Of_Radnom_Matched_sites.R
Generate random atched control integration sites, 10 for each read integration site.
Add values of epigenetic features averaged for 1kb bins to each site.
