
<!-- README.md is generated from README.Rmd. Please edit that file -->

# msp-ms

<!-- badges: start -->
<!-- badges: end -->

The goal of mspms is provide a concise code-base for the normalization
and data processing required to analyze data from the [Multiplex
Substrate Profiling by Mass Spectrometry (MSP-MS)
method](https://pubmed.ncbi.nlm.nih.gov/36948708/).

Additionally, we provide a [graphical user interface powered by shiny
apps](https://gonzalezlab.shinyapps.io/mspms_shiny/) that allows for a
user to utilize the method without requiring any R coding knowledge.

## Installation

You can install the released version of mspms from github

``` r
devtools::install_github("baynec2/mspms")
```

## Workflow

So how does the msps data normalization process work?

1.  Takes two input files from PEAKS and combines them
2.  Normalizes values and then does a reverse log2 transformation
3.  Looks for outliers across replicates. Removes them.
4.  Imputes data for missing values (likely to be very low intesity).
5.  Figures out the locations of the detected clevages within the
    library of peptide sequences used. Is it cleaved at the N or C
    terminus, or both?
6.  Calculates the fold change and the p/q value across experimental
    conditions using T-tests.

### Combining Peaks file outputs

The files coming from peaks should have….. **Describe how to do this
here**.

Importantly, this. data contains… what????

``` r
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(mspms)

### Loading the files ###
lfq_filename = "legacy/protein-peptides-lfq.csv"
#file "protein-peptides.csv" exported from PEAKS identification
id_filename = "legacy/protein-peptides-id.csv"

# Prepare the data for normalyzer analysis
prepared_data = prepare_for_normalyzer(lfq_filename,id_filename)
#> Warning: One or more parsing issues, call `problems()` on your data frame for details,
#> e.g.:
#>   dat <- vroom(...)
#>   problems(dat)
#> Rows: 1099 Columns: 46
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr  (6): Protein Accession, Peptide, Used, Candidate, Sample Profile (Ratio...
#> dbl (40): Protein Group, Protein ID, Quality, Significance, Avg. ppm, Avg. A...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> Rows: 1381 Columns: 66
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr  (4): Protein Accession, Peptide, Unique, Source File
#> dbl (62): Protein Group, Protein ID, -10lgP, Mass, Length, ppm, m/z, z, RT, ...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

### Normalizing values.

msp-ms uses the …. package to do normalization under the hood.

First we need to extract a design matrix.

Then we can normalyze the data.

``` r
# Extracting design matrix
design_matrix = extract_design_matrix(prepared_data)

# normalyze the data 

source("R/normalyze.R")

normalyzed_data = normalyze(prepared_data,design_matrix)
#> You are running version 1.19.7 of NormalyzerDE
#> [Step 1/5] Load data and verify input
#> Input data checked. All fields are valid.
#> Sample check: More than one sample group found
#> Sample replication check: All samples have replicates
#> RT annotation column found (23)
#> [Step 1/5] Input verified, job directory prepared at:./2024-04-12_msp-ms_normalyze_output
#> [Step 2/5] Performing normalizations
#> [Step 2/5] Done!
#> [Step 3/5] Generating evaluation measures...
#> [Step 3/5] Done!
#> [Step 4/5] Writing matrices to file
#> [Step 4/5] Matrices successfully written
#> [Step 5/5] Generating plots...
#> [Step 5/5] Plots successfully generated
#> All done! Results are stored in: ./2024-04-12_msp-ms_normalyze_output, processing time was 0.4 minutes
#> Rows: 835 Columns: 53
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: "\t"
#> chr  (6): Protein Accession, Peptide, Used, Candidate, Sample Profile (Ratio...
#> dbl (47): Protein Group, Protein ID, Quality, Significance, Avg. ppm, Avg. A...
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

### Handling Outliers

Here we use a dixon test from the outliers package to detect outliers
from each of our replicates.

``` r
outliers = handle_outliers(normalyzed_data)
```

### Imputation of data

We have a lot of missing, or 0 values. For these, we need to impute them
so we can do downstream statistics Data is imputated by ….

``` r
imputed = impute(outliers)
```

### Joining with Library

Next we need to join everything with the sequences of the peptide
library

``` r
joined_with_library = join_with_library(imputed)
```

### Calcuating clevages

Next, we need to determine where the peptide sequences are cleaved.

We check both the N and C terminus.

Sequences are presented as the 4 amino acids on both sides of a
cleavage. X indicates that there was nothing on that side in the library
because it was cleaved close to the edge.

``` r
final_data = add_cleavages(joined_with_library)
```

### Conducting stats

Right now this is left to the user to do on their own. Will update after
talkng to Bri.
