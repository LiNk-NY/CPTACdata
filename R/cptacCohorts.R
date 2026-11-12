#' CPTAC cohorts in the PDC
#'
#' The CPTAC patient cohorts in the Proteomic Data Commons (PDC) that
#' [getCPTAC()] can download, with one row per PDC study. A cohort, such as
#' the CPTAC-3 clear cell renal cell carcinoma discovery cohort, is profiled
#' in several PDC studies, one for each experiment (proteome,
#' phosphoproteome, acetylome, ...). [getCPTAC()] combines them in one
#' `MultiAssayExperiment`.
#'
#' The table covers the PDC projects "CPTAC3 Discovery and Confirmatory" and
#' "CPTAC-2" (the Prospective collections). Studies of reference samples
#' only ("CompRef") and method studies are left out. Cohorts run on TCGA
#' tumors are in [cptacProjects] and are added to TCGA data with
#' [appendCPTAC()].
#'
#' Some studies, such as most metabolome and lipidome studies, have no
#' quantification matrix available through the PDC API. They are listed with
#' `data_type = NA` and are skipped by [getCPTAC()].
#'
#' @format A `data.frame` with one row per PDC study and 8 columns:
#'   - `cohort`: Cohort name, the PDC study name without the experiment
#'     (e.g. `"CPTAC CCRCC Discovery Study"`), or a curated name for the
#'     CPTAC-2 Prospective studies (e.g. `"Prospective Breast"`).
#'   - `experiment`: Name of the experiment in the `MultiAssayExperiment`,
#'     from the rest of the study name (e.g. `"Proteome"`,
#'     `"DIA_Proteome"`).
#'   - `pdc_study_id`: PDC study ID.
#'   - `pdc_name`: PDC study name.
#'   - `pdc_project`: PDC project.
#'   - `analytical_fraction`: PDC analytical fraction, e.g. Proteome or
#'     Phosphoproteome.
#'   - `experiment_type`: PDC experiment type, e.g. `"TMT10"` or
#'     `"Label Free"`.
#'   - `data_type`: `quantDataMatrix` data type that is downloaded; `NA` if
#'     the PDC has no matrix for the study.
#'
#' @source Built from the PDC GraphQL API (<https://pdc.cancer.gov/graphql>)
#'   by `data-raw/cptacCohorts.R` in the package sources. The cohort and
#'   experiment names of the CPTAC-2 studies, and the excluded studies, are
#'   curated in that script.
#'
#' @examples
#' data("cptacCohorts", package = "CPTACdata")
#' unique(cptacCohorts[["cohort"]])
#' cptacCohorts[
#'     cptacCohorts[["cohort"]] == "CPTAC CCRCC Discovery Study",
#'     c("experiment", "pdc_study_id", "data_type")
#' ]
#'
#' @docType data
#' @keywords datasets
"cptacCohorts"
