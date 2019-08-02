# ElifeROC
ROC analysis for Samples+-BIX

## 1_PreparingIS_and_EpigenTables.R
First, cut the genome into 1kb bins. Remove bins in blacklisted regions. Add distance to nearest protein coding gene. Calculate average score for epigenetic features on this regions. 

Add the values to integration sites +- bix.
https://github.com/UNIZG-PMF-Bioinfo/ElifeROC/blob/master/1_PreparingIS_and_EpigenTables.R

## 2_Generation_Of_Radnom_Matched_sites.R
Generate random atched control integration sites, 10 for each read integration site.

# Comparison of Epigenetic marks scores on CTRL IS VS BIX IS:

![](https://github.com/UNIZG-PMF-Bioinfo/ElifeROC/blob/master/p%20values%20CTRL%20VS%20BIX%20IS%20Epigens.PNG)
