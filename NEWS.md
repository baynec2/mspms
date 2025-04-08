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

## mspms 0.99.3

* Removed some redundant tests of the log2fc_t_test() function that took a long
time to perform and were causing timeouts on the Bioconductor build.


## mspms 0.99.4

* Removed the cathepsin B-D data from the internal datasets and from extdata/
. This keeps the example run time down, as well as greatly reduces the
overall size of the package. 

## mspms 0.99.5

* Removed interactive heatmap from vignette to get rid of warning that vignette 
size was too large. Minor changes to vignette aesthetics.  

## mspms 0.99.6

* Improved vignette according to review feedback. Added peptide length 
calculation to row data. 

## mspms 0.99.7

* Added ability to remove dendrogram plotting from plot_heatmap(). This was 
done because the dendrogram plotting feature in heatmaply is recursive and 
uses a ton of memory which was problematic for the shiny instance.

## mspms 0.99.8

* Added support for limma based statistics via limma_stats() function

