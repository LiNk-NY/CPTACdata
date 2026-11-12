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
- `getCPTAC()` downloads a CPTAC-3 or CPTAC-2 Prospective cohort from the
  PDC as a `MultiAssayExperiment`, with one experiment per PDC study (e.g.
  proteome and phosphoproteome) and the PDC clinical data of its patients as
  `colData`. `dry_run = TRUE` reports what would be downloaded.
- The `cptacCohorts` data set lists the CPTAC cohorts that `getCPTAC()` can
  download, with one row per PDC study.
- The vignette "Adding CPTAC proteomics to a TCGA MultiAssayExperiment"
  walks through a complete example with an ovarian cancer object.
- The vignette "Downloading CPTAC cohorts as MultiAssayExperiment objects"
  shows how to list the cohorts and download one with `getCPTAC()`.
