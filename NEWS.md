# mspms 0.99

## mspms 0.99.0

* Initial Bioconductor submission.

## mspms 0.99.1

* Major edits to the structure of the package based on review.
Now utilizes QFeatures objects and MSCoreutils to perform data processing steps. 

## mspms 0.99.2

* Removed logo from vignette to pass Bioconductor build. Fixed bugs to 
nterm_cleavages() and cterm_cleavages() that were causing problems when 
n_residue != 4. Simplified the logic behind these functions and added extensive 
testing to confirm that they work as expected with a range of peptide sequences.
Minor change to the plot_volcano() function allowing the user to specify a 
padj threshold.