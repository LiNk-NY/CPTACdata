# Changes in CPTACdata 1.0.0

## New Features

- `appendCPTAC()` adds CPTAC proteomics from the Proteomic Data Commons (PDC)
  to a `MultiAssayExperiment` from `curatedTCGAData`. Only patients already in
  the object are added, and each new experiment comes with a `sampleMap` that
  links its samples to TCGA barcodes.
- `appendCPTAC()` supports the Proteome, Phosphoproteome and Glycoproteome
  fractions, and `dry_run = TRUE` reports how many patients each PDC study
  shares with the object without downloading quantification matrices.
- `appendCPTAC()` keeps every technical replicate in the PDC matrix as its
  own column. Use `replicates = "first"` to keep only the first replicate of
  each sample.
- The `cptacProjects` data set lists the TCGA projects with CPTAC data
  generated on the same tumors, the PDC study for each one, and whether every
  qualified TCGA patient in the study was found in the GDC.
- The vignette "Adding CPTAC proteomics to a TCGA MultiAssayExperiment"
  walks through a complete example with an ovarian cancer object.
