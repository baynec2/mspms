
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

### Loading design matrix

Before we normalyze data, we need to know what samples are in what
groups. We can do that by defining a design matrix.

``` r
design_matrix = readr::read_csv("annotation.csv")
#> Rows: 24 Columns: 2
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (2): sample, group
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

head(design_matrix)
#> # A tibble: 6 × 2
#>   sample      group   
#>   <chr>       <chr>   
#> 1 DMSO_T000_1 DMSO_T0 
#> 2 DMSO_T000_2 DMSO_T0 
#> 3 DMSO_T000_3 DMSO_T0 
#> 4 DMSO_T000_4 DMSO_T0 
#> 5 DMSO_T060_1 DMSO_T60
#> 6 DMSO_T060_2 DMSO_T60
```

### Normalyzing data

msp-ms uses the …. package to do normalization under the hood.

Now we can normalyze the data.

``` r
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
from each of our replicates. We need to know what samples are part of
which groups, so we need to specify the design matrix here too. Make
sure that the column header names are “sample” and “group” just like
before.

``` r
design_matrix = readr::read_csv("annotation.csv")
#> Rows: 24 Columns: 2
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (2): sample, group
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
outliers = handle_outliers(normalyzed_data,design_matrix)
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

head(final_data)
#>    library_reference_id library_real_sequence         Peptide    nterm
#> 1 TDP1|TDP1|generation1        LVATVYEFGHIDHM  LVATVYEFGHIDHL     <NA>
#> 2 TDP1|TDP1|generation1        LVATVYEFGHIDHM L_VATVYEFGHIDHL XXXLVATV
#> 3 TDP1|TDP1|generation1        LVATVYEFGHIDHM L_VATVYEFGHIDHL XXXLVATV
#> 4 TDP1|TDP1|generation1        LVATVYEFGHIDHM    T_VYEFGHIDHL LVATVYEF
#> 5 TDP1|TDP1|generation1        LVATVYEFGHIDHM   A_TVYEFGHIDHL XLVATVYE
#> 6 TDP1|TDP1|generation1        LVATVYEFGHIDHM  LVATVYEFGHID_H     <NA>
#>   nterm_cleavage_pos    cterm cterm_cleavage_pos library_match_sequence
#> 1                 NA     <NA>                 NA         LVATVYEFGHIDHL
#> 2                  1     <NA>                 NA         LVATVYEFGHIDHL
#> 3                  1     <NA>                 NA         LVATVYEFGHIDHL
#> 4                  4     <NA>                 NA         LVATVYEFGHIDHL
#> 5                  3     <NA>                 NA         LVATVYEFGHIDHL
#> 6                 NA GHIDHMXX                 12         LVATVYEFGHIDHL
#>   Protein Group Protein ID Peptide_no_cleavage Used Candidate Quality
#> 1            12        383      LVATVYEFGHIDHL    Y         Y    1.56
#> 2            12        383      LVATVYEFGHIDHL    N         Y    1.16
#> 3            12        383      LVATVYEFGHIDHL    N         Y    0.49
#> 4            12        383         TVYEFGHIDHL    N         Y    0.61
#> 5            12        383        ATVYEFGHIDHL    N         Y    0.88
#> 6            12        383       LVATVYEFGHIDH    N         Y    0.45
#>   Significance Avg. ppm Avg. Area
#> 1         0.82      1.6  1.68e+08
#> 2         1.10      1.3  2.88e+06
#> 3         0.78      1.5  1.41e+06
#> 4        60.00      0.7  6.71e+05
#> 5        60.00      1.6  4.54e+04
#> 6         1.86      2.2  1.85e+05
#>                                                                                                    Sample Profile (Ratio)
#> 1 1.00:0.92:0.75:0.60:0.75:0.88:0.71:0.41:0.69:0.57:1.05:0.65:0.88:1.00:0.88:0.77:0.81:0.90:0.67:0.80:0.82:1.30:0.72:0.81
#> 2 1.00:1.16:0.67:0.51:0.73:0.58:0.68:0.40:0.69:0.56:0.90:0.62:0.60:0.86:0.95:0.76:0.89:0.73:0.58:0.69:0.81:1.09:0.77:0.76
#> 3 1.00:0.94:0.80:0.64:0.82:0.94:0.74:0.48:1.13:0.85:1.07:0.63:0.88:0.84:0.90:0.57:1.02:0.75:0.66:0.75:0.84:1.65:1.00:1.01
#> 4                                                   0:0:0:0:0:0:1.00:0.83:16.64:10.26:14.99:15.20:0:0:0:0:0:0:0:0:0:0:0:0
#> 5                                           0:1.00:0:1.21:1.01:0.81:1.23:0.76:0:1.27:0:0:0:0:0:0:0:0:1.04:1.61:0:0:1.49:0
#> 6    1.00:1.07:0.92:0.05:1.09:1.24:0.51:0.31:1.18:1.63:2.50:0.86:1.66:0.34:1.05:0.50:1.21:0.46:1.83:0.60:0.44:0:1.35:0.69
#>    Group 1  Group 2  Group 3  Group 4  Group 5  Group 6
#> 1 1.70e+08 1.43e+08 1.54e+08 1.84e+08 1.66e+08 1.91e+08
#> 2 3.18e+06 2.27e+06 2.65e+06 3.02e+06 2.75e+06 3.27e+06
#> 3 1.37e+06 1.21e+06 1.50e+06 1.29e+06 1.29e+06 1.82e+06
#> 4       NA 1.25e+05 3.90e+06       NA       NA       NA
#> 5 5.27e+04 9.08e+04 3.02e+04       NA 6.32e+04 3.56e+04
#> 6 1.50e+05 1.56e+05 3.05e+05 1.76e+05 2.03e+05 1.22e+05
#>           Group Profile (Ratio) Max Ratio #Vector Start End    RT -10lgP
#> 1 1.00:0.84:0.91:1.08:0.97:1.12      1.33       6     1  14 86.93 155.13
#> 2 1.00:0.71:0.83:0.95:0.87:1.03      1.44       2     2  14 82.05 104.01
#> 3 1.00:0.88:1.09:0.94:0.94:1.33      1.51       4     2  14 82.05 104.01
#> 4            0:1.00:31.22:0:0:0    256.00       2     5  14 77.84  64.99
#> 5    1.00:1.72:0.57:0:1.20:0.67    256.00       1     4  14 83.12 100.25
#> 6 1.00:1.03:2.03:1.17:1.35:0.81      2.50       3     1  12 79.10 126.76
#>       Mass Length  ppm      m/z z DMSO_T060_1 DMSO_T060_2 DMSO_T060_3
#> 1 1612.825     14 -7.5 807.4137 2 249501444.0 241745551.6 272260347.0
#> 2 1499.741     13 -6.2 500.9178 3   4430250.0   2890435.9   4796968.0
#> 3 1499.741     13 -6.2 500.9178 3   2143153.4   1997028.5   2222533.4
#> 4 1228.588     10 -6.1 615.2974 2      4869.0      8188.0    505626.4
#> 5 1371.646     11 -5.3 686.8266 2    153219.5    101690.8    218549.1
#> 6 1362.682     12 -5.9 682.3442 2    343864.2    323203.3    187063.2
#>   DMSO_T060_4 DMSO_T240_1 DMSO_T240_2 DMSO_T240_3 DMSO_T240_4 MZB_T240_2
#> 1 186722266.4 153354807.7 123062500.0 188028052.3 159479618.8  157146302
#> 2   3324889.4   2800855.2   2204869.8   2962085.8   2779166.9    2412283
#> 3   1724098.3   1948884.0   1425474.0   1485335.8   1207823.6    1548268
#> 4    497632.5   4845585.9   2871458.3   3520159.9   4866473.7       8252
#> 5    159638.7      5976.0    124088.0      4557.0     12087.0       6457
#> 6    134096.5    249201.6    330217.7    424994.9    199349.5      12246
#>    MZB_T240_3   MZB_T240_1  MZB_T240_4 DMSO_T000_1  DMSO_T000_2 DMSO_T000_3
#> 1 143555311.2 149433035.71 162315194.0 139446122.0 174110648.15 170226354.7
#> 2   2776036.5   2716165.18   2768906.2   2542056.1   4038273.15   2782546.2
#> 3   1540129.8   1195464.29   1565864.2   1080874.2   1394708.33   1418553.0
#> 4     12414.0      8002.00      8620.0      6213.0      5282.00      4360.0
#> 5    134999.0     11678.00      7989.0     11261.0     86873.01      5561.0
#> 6    253836.2     75683.44    129852.2    132106.9    193253.70    198597.4
#>    DMSO_T000_4  MZB_T000_1   MZB_T000_2  MZB_T000_3   MZB_T000_4  MZB_T060_1
#> 1 139358447.49 187670312.5 155657432.43 141737019.2 157803152.65 139685121.3
#> 2   2169046.80   2338187.5   2462081.50   2811504.8   2842417.04   2810233.2
#> 3   1168812.79   1466494.8   1017760.14   1130798.1    901732.30   1372055.0
#> 4     10179.00      4969.0     10172.00      9968.0      5995.00      8444.0
#> 5    130367.58     10300.0      4244.00      8566.0      6852.00      9854.0
#> 6     10294.54    336370.8     51037.68    161874.5     96838.21    198369.4
#>     MZB_T060_2  MZB_T060_3  MZB_T060_4         peptide
#> 1 184267035.40 206092620.5 209471050.4  LVATVYEFGHIDHL
#> 2   2724799.78   3276724.4   3311398.6 L_VATVYEFGHIDHL
#> 3   1195775.44   1586468.4   1530267.6 L_VATVYEFGHIDHL
#> 4      9063.00      8685.0      9454.0    T_VYEFGHIDHL
#> 5      9129.00    147378.5    191910.6   A_TVYEFGHIDHL
#> 6     89095.07    535247.7    149263.8  LVATVYEFGHID_H
```

### Conducting stats

Right now this is left to the user to do on their own. Will update after
talkng to Bri.
