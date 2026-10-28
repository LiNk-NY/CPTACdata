#' TCGA projects with matching CPTAC proteomics in the PDC
#'
#' A crosswalk between TCGA projects (e.g. \code{"TCGA-BRCA"}) and the
#' Proteomic Data Commons (PDC) studies holding CPTAC data generated on the
#' \emph{same TCGA tumors} (the PDC "CPTAC-TCGA" project). It is used by
#' \code{\link{appendCPTAC}} to select studies.
#'
#' Cancer types profiled by CPTAC-3 (e.g. LUAD, ccRCC) and the "Prospective"
#' CPTAC-2 collections are not included: those cohorts use new patients
#' rather than TCGA samples, so there is nothing to match against a TCGA
#' MultiAssayExperiment.
#'
#' One study can cover more than one TCGA project (PDC000111 holds colon and
#' rectum tumors); \code{pdc_disease_type} separates them. The JHU and PNNL
#' ovarian proteomes share 33 patients and are separate experiments.
#'
#' Some TCGA-derived studies (the colon proteome, PDC000111, and the ovarian
#' glycoproteome, PDC000112) have biospecimen records in the PDC but no
#' quantification matrix available through its API. They are listed with
#' \code{data_type = NA} and are skipped by \code{\link{appendCPTAC}}.
#'
#' @format A \code{data.frame} with one row per TCGA project and PDC study
#'   (8 rows) and 8 columns:
#'   \describe{
#'     \item{tcga_project}{TCGA project ID.}
#'     \item{pdc_study_id}{PDC study ID.}
#'     \item{pdc_name}{PDC study name.}
#'     \item{experiment}{Name given to the experiment added to the
#'       \code{MultiAssayExperiment}.}
#'     \item{analytical_fraction}{Proteome, Phosphoproteome or
#'       Glycoproteome.}
#'     \item{data_type}{\code{quantDataMatrix} data type that is downloaded;
#'       \code{NA} if the PDC has no matrix for the study.}
#'     \item{pdc_disease_type}{PDC disease label used to split studies that
#'       span several TCGA projects; \code{NA} if not needed.}
#'     \item{verified}{\code{TRUE} if every qualified TCGA patient in the
#'       study was found in the GDC project when the data were built.}
#'   }
#'
#' @source Built from the PDC GraphQL API (\url{https://pdc.cancer.gov/graphql})
#'   and the GDC API (\url{https://api.gdc.cancer.gov}) by
#'   \code{data-raw/cptacProjects.R} in the package sources. Experiment names
#'   and data types are curated in that script.
#'
#' @examples
#' cptacProjects
#'
#' @docType data
#' @keywords datasets
"cptacProjects"
