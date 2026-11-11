
# Overview

The Clinical Proteomic Tumor Analysis Consortium (CPTAC) profiled the
proteomes of several hundred TCGA tumors. `CPTACdata` downloads those
data from the Proteomic Data Commons (PDC) and appends them to a
`MultiAssayExperiment` from `curatedTCGAData`.

``` r
library(CPTACdata)
library(curatedTCGAData)
library(MultiAssayExperiment)
```

# Which TCGA projects have CPTAC data?

``` r
cptacProjects
#>   tcga_project pdc_study_id                                          pdc_name
#> 1    TCGA-BRCA    PDC000173                       TCGA Breast Cancer Proteome
#> 2    TCGA-BRCA    PDC000174                TCGA Breast Cancer Phosphoproteome
#> 3    TCGA-COAD    PDC000111                        TCGA Colon Cancer Proteome
#> 4    TCGA-READ    PDC000111                        TCGA Colon Cancer Proteome
#> 5      TCGA-OV    PDC000112                    TCGA Ovarian JHU Glycoproteome
#> 6      TCGA-OV    PDC000113                         TCGA Ovarian JHU Proteome
#> 7      TCGA-OV    PDC000114                        TCGA Ovarian PNNL Proteome
#> 8      TCGA-OV    PDC000115 TCGA Ovarian PNNL Phosphoproteome Velos Qexactive
#>                   experiment analytical_fraction  data_type
#> 1             CPTAC_Proteome            Proteome log2_ratio
#> 2      CPTAC_Phosphoproteome     Phosphoproteome log2_ratio
#> 3             CPTAC_Proteome            Proteome       <NA>
#> 4             CPTAC_Proteome            Proteome       <NA>
#> 5    CPTAC_JHU_Glycoproteome       Glycoproteome       <NA>
#> 6         CPTAC_JHU_Proteome            Proteome log2_ratio
#> 7        CPTAC_PNNL_Proteome            Proteome log2_ratio
#> 8 CPTAC_PNNL_Phosphoproteome     Phosphoproteome log2_ratio
#>        pdc_disease_type verified
#> 1                  <NA>     TRUE
#> 2                  <NA>     TRUE
#> 3  Colon Adenocarcinoma     TRUE
#> 4 Rectum Adenocarcinoma     TRUE
#> 5                  <NA>     TRUE
#> 6                  <NA>     TRUE
#> 7                  <NA>     TRUE
#> 8                  <NA>     TRUE
```

# Appending to a MultiAssayExperiment

Start with any `curatedTCGAData` object for one TCGA project. Only
patients in its `colData` can be matched, so include enough assays to
cover the patients of interest.

``` r
mae <- curatedTCGAData("OV", "RPPAArray", version = "2.1.1", dry.run = FALSE)
```

A dry run reports how many patients each PDC study shares with the
object, downloading only sample metadata:

``` r
appendCPTAC(mae, dry_run = TRUE)
#>            experiment pdc_study_id quant_available mae_patients pdc_patients
#> 1  CPTAC_JHU_Proteome    PDC000113            TRUE          420          122
#> 2 CPTAC_PNNL_Proteome    PDC000114            TRUE          420           84
#>   matched_patients
#> 1              119
#> 2               81
```

Then add the data. The project is detected from `colData`, and the
**proteome is added by default**. Call `appendCPTAC()` again for other
assay types:

``` r
mae <- appendCPTAC(mae)
#> Adding CPTAC_JHU_Proteome from PDC000113 (119 patients)...
#> Adding CPTAC_PNNL_Proteome from PDC000114 (81 patients)...
mae <- appendCPTAC(mae, fraction = "Phosphoproteome")
#> Adding CPTAC_PNNL_Phosphoproteome from PDC000115 (66 patients)...
mae
#> A MultiAssayExperiment object of 4 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 4:
#>  [1] OV_RPPAArray-20160128: SummarizedExperiment with 208 rows and 427 columns
#>  [2] CPTAC_JHU_Proteome: SummarizedExperiment with 8597 rows and 119 columns
#>  [3] CPTAC_PNNL_Proteome: SummarizedExperiment with 7480 rows and 81 columns
#>  [4] CPTAC_PNNL_Phosphoproteome: SummarizedExperiment with 21478 rows and 66 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
```

The new experiments have one column per TCGA sample, and the `sampleMap`
links them to patients:

``` r
tail(sampleMap(mae))
#> DataFrame with 6 rows and 3 columns
#>                             assay      primary          colname
#>                          <factor>  <character>      <character>
#> [688,] CPTAC_PNNL_Phosphoproteome TCGA-61-1907 TCGA-61-1907-01A
#> [689,] CPTAC_PNNL_Phosphoproteome TCGA-61-1919 TCGA-61-1919-01A
#> [690,] CPTAC_PNNL_Phosphoproteome TCGA-61-1995 TCGA-61-1995-01A
#> [691,] CPTAC_PNNL_Phosphoproteome TCGA-61-2008 TCGA-61-2008-01A
#> [692,] CPTAC_PNNL_Phosphoproteome TCGA-61-2095 TCGA-61-2095-01A
#> [693,] CPTAC_PNNL_Phosphoproteome TCGA-61-2612 TCGA-61-2612-01A
mae[["CPTAC_PNNL_Proteome"]]
#> class: SummarizedExperiment 
#> dim: 7480 81 
#> metadata(4): pdc_study_id pdc_name data_type source
#> assays(1): log2_ratio
#> rownames(7480): A1BG A2M ... ZYX ZZEF1
#> rowData names(0):
#> colnames(81): TCGA-09-1664-01A TCGA-13-1484-01A ... TCGA-61-2096-01A
#>   TCGA-61-2612-01A
#> colData names(4): case_submitter_id aliquot_submitter_id sample_type
#>   pdc_column
```

The JHU and PNNL ovarian proteomes are separate experiments that share
some patients:

``` r
upsetSamples(mae)
```

<img src="/home/mramos/gh/CPTACdata/README_files/figure-gfm/unnamed-chunk-14-1.png" alt="" width="100%" />

# Session info

``` r
sessionInfo()
#> R version 4.6.1 Patched (2026-07-12 r90246)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 26.04.1 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.1 
#> LAPACK: /home/mramos/src/svn/r-4-6/inst/lib/R/lib/libRlapack.so;  LAPACK version 3.12.1
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: America/New_York
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] curatedTCGAData_1.35.0      MultiAssayExperiment_1.39.2
#>  [3] SummarizedExperiment_1.43.0 Biobase_2.73.2             
#>  [5] GenomicRanges_1.65.4        Seqinfo_1.3.2              
#>  [7] IRanges_2.47.5              S4Vectors_0.51.10          
#>  [9] BiocGenerics_0.59.12        generics_0.1.4             
#> [11] MatrixGenerics_1.25.0       matrixStats_1.5.0          
#> [13] CPTACdata_0.99.4           
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyselect_1.2.1          farver_2.1.2             
#>  [3] dplyr_1.2.1               blob_1.3.0               
#>  [5] S7_0.2.2                  filelock_1.0.3           
#>  [7] Biostrings_2.81.9         bitops_1.1-0             
#>  [9] fastmap_1.2.0             RCurl_1.98-1.20          
#> [11] BiocFileCache_3.3.2       GenomicAlignments_1.49.2 
#> [13] XML_3.99-0.25             digest_0.6.39            
#> [15] lifecycle_1.0.5           KEGGREST_1.53.6          
#> [17] RSQLite_3.53.3            magrittr_2.0.5           
#> [19] compiler_4.6.1            rlang_1.3.0              
#> [21] tools_4.6.1               yaml_2.3.12              
#> [23] rtracklayer_1.73.0        knitr_1.52               
#> [25] labeling_0.4.3            S4Arrays_1.13.2          
#> [27] bit_4.6.0                 curl_8.1.0               
#> [29] DelayedArray_0.39.8       plyr_1.8.9               
#> [31] RColorBrewer_1.1-3        xml2_1.6.0               
#> [33] abind_1.4-8               BiocParallel_1.47.0      
#> [35] withr_3.0.3               grid_4.6.1               
#> [37] ExperimentHub_3.3.2       ggplot2_4.0.3            
#> [39] scales_1.4.0              dichromat_2.0-1          
#> [41] cli_3.6.6                 UpSetR_1.4.1             
#> [43] rmarkdown_2.32            crayon_1.5.3             
#> [45] otel_0.3.0                rstudioapi_0.19.0        
#> [47] httr_1.4.9                tzdb_0.5.0               
#> [49] rjson_0.2.23              BiocBaseUtils_1.15.1     
#> [51] GenomicDataCommons_1.37.0 DBI_1.3.0                
#> [53] cachem_1.1.0              stringr_1.6.0            
#> [55] rvest_1.0.5               parallel_4.6.1           
#> [57] AnnotationDbi_1.75.2      BiocManager_1.30.27      
#> [59] XVector_0.53.0            restfulr_0.0.17          
#> [61] vctrs_0.7.3               Matrix_1.7-6             
#> [63] jsonlite_2.0.0            hms_1.1.4                
#> [65] bit64_4.8.6               GenomicFeatures_1.65.0   
#> [67] glue_1.8.1                codetools_0.2-20         
#> [69] gtable_0.3.6              stringi_1.8.9            
#> [71] BiocVersion_3.24.0        GenomeInfoDb_1.49.1      
#> [73] BiocIO_1.23.3             UCSC.utils_1.9.0         
#> [75] tibble_3.3.1              pillar_1.11.1            
#> [77] rappdirs_0.3.4            htmltools_0.5.9          
#> [79] R6_2.6.1                  dbplyr_2.6.0             
#> [81] httr2_1.3.0               evaluate_1.0.5           
#> [83] lattice_0.23-1            readr_2.2.0              
#> [85] AnnotationHub_4.3.3       cigarillo_1.3.1          
#> [87] png_0.1-9                 Rsamtools_2.29.0         
#> [89] memoise_2.0.1             Rcpp_1.1.2               
#> [91] gridExtra_2.3.1           SparseArray_1.13.4       
#> [93] xfun_0.62                 TCGAutils_1.33.3         
#> [95] pkgconfig_2.0.3
```
