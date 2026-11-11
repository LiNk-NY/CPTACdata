#' CPTACdata: CPTAC proteomics for TCGA MultiAssayExperiments
#'
#' Adds proteomic data from the Clinical Proteomic Tumor Analysis Consortium
#' (CPTAC), held in the NCI Proteomic Data Commons (PDC), to a
#' `MultiAssayExperiment` from `curatedTCGAData`. The main function is
#' [appendCPTAC()]. The studies it can use are listed in [cptacProjects].
#'
#' @section Data use agreement:
#' The PDC data are provided under the PDC Data Use Agreement, available at
#' <https://proteomic.datacommons.cancer.gov/pdc/data-use-guidelines>. Before
#' the package makes its first request to the PDC, it asks you to read these
#' terms and accept them. Your answer is kept only for the current R session.
#'
#' In a non-interactive session, such as a script or a batch job, the
#' package stops with an error until you accept the terms. To accept them,
#' read the agreement and run `options(CPTACdata.acceptDUA = TRUE)` before
#' calling the package.
#'
#' @section Licence and attribution:
#' The package code is licensed under the Artistic-2.0 licence (see the
#' `License` field of `DESCRIPTION`). The PDC data are not covered by that
#' licence. Check the PDC Data Use Agreement for the terms on reuse and for
#' any acknowledgement or citation you must include in publications.
#'
#' @keywords internal
"_PACKAGE"
