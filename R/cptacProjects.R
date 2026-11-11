#' TCGA projects with matching CPTAC proteomics in the PDC
#'
#' A crosswalk between TCGA projects (e.g. `"TCGA-BRCA"`) and the
#' Proteomic Data Commons (PDC) studies holding CPTAC data generated on the
#' *same TCGA tumors* (the PDC "CPTAC-TCGA" project). It is used by
#' [appendCPTAC()] to select studies.
#'
#' Cancer types profiled by CPTAC-3 (e.g. LUAD, ccRCC) and the "Prospective"
#' CPTAC-2 collections are not included: those cohorts use new patients
#' rather than TCGA samples, so there is nothing to match against a TCGA
#' MultiAssayExperiment.
#'
#' One study can cover more than one TCGA project (PDC000111 holds colon and
#' rectum tumors); `pdc_disease_type` separates them. The JHU and PNNL
#' ovarian proteomes share 33 patients and are separate experiments.
#'
#' Some TCGA-derived studies (the colon proteome, PDC000111, and the ovarian
#' glycoproteome, PDC000112) have biospecimen records in the PDC but no
#' quantification matrix available through its API. They are listed with
#' `data_type = NA` and are skipped by [appendCPTAC()].
#'
#' @format A `data.frame` with one row per TCGA project and PDC study
#'   (8 rows) and 8 columns:
#'   - `tcga_project`: TCGA project ID.
#'   - `pdc_study_id`: PDC study ID.
#'   - `pdc_name`: PDC study name.
#'   - `experiment`: Name given to the experiment added to the
#'     `MultiAssayExperiment`.
#'   - `analytical_fraction`: Proteome, Phosphoproteome or Glycoproteome.
#'   - `data_type`: `quantDataMatrix` data type that is downloaded; `NA` if
#'     the PDC has no matrix for the study.
#'   - `pdc_disease_type`: PDC disease label used to split studies that span
#'     several TCGA projects; `NA` if not needed.
#'   - `verified`: `TRUE` if every qualified TCGA patient in the study was
#'     found in the GDC project when the data were built.
#'
#' @source Built from the PDC GraphQL API (<https://pdc.cancer.gov/graphql>)
#'   and the GDC API (<https://api.gdc.cancer.gov>) by
#'   `data-raw/cptacProjects.R` in the package sources. Experiment names and
#'   data types are curated in that script.
#'
#' @examples
#' cptacProjects
#'
#' @docType data
#' @keywords datasets
"cptacProjects"
