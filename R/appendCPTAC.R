#' Append CPTAC proteomics from the PDC to a TCGA MultiAssayExperiment
#'
#' Takes a \code{MultiAssayExperiment} from \code{curatedTCGAData}, finds the
#' CPTAC studies in the Proteomic Data Commons (PDC) that were run on the same
#' TCGA tumors, and adds them as new experiments with a correct
#' \code{sampleMap}.
#'
#' Only patients already in \code{colData(mae)} are added. Specimens are kept
#' if they carry a TCGA barcode, are "Qualified" and not pooled. Each new
#' experiment has one column per TCGA sample (e.g. \code{TCGA-24-1555-01A});
#' technical replicates are reduced to the first one. Missing values are
#' \code{NA}. The quantification is
#' the PDC Common Data Analysis Pipeline value given by the \code{data_type}
#' column of \code{\link{cptacProjects}} (log2 ratios to a common
#' reference).
#'
#' Biospecimen tables are cached for the R session, so a \code{dry_run}
#' followed by a real run downloads them only once.
#'
#' @param mae A \code{MultiAssayExperiment} for a single TCGA project, such
#'   as one from \code{curatedTCGAData}.
#' @param project TCGA project, e.g. \code{"OV"} or \code{"TCGA-OV"}. If
#'   \code{NULL}, it is read from the \code{admin.disease_code} column of
#'   \code{colData(mae)}.
#' @param fraction The assay type to add: one of \code{"Proteome"}
#'   (default), \code{"Phosphoproteome"} or \code{"Glycoproteome"}. Call
#'   \code{appendCPTAC} once per type to add several.
#' @param dry_run If \code{TRUE}, download only sample metadata and return a
#'   one-row-per-study summary of what would be added.
#'
#' @return The input \code{mae} with one new experiment per PDC study (named
#'   as in the \code{experiment} column of \code{\link{cptacProjects}}), or a
#'   summary \code{data.frame} if \code{dry_run = TRUE}.
#'
#' @seealso \code{\link{cptacProjects}}
#'
#' @examples
#' \dontrun{
#' mae <- curatedTCGAData::curatedTCGAData(
#'     "OV", "RNASeq2GeneNorm", version = "2.1.1", dry.run = FALSE
#' )
#' appendCPTAC(mae, dry_run = TRUE)
#' mae <- appendCPTAC(mae)
#' mae <- appendCPTAC(mae, fraction = "Phosphoproteome")
#' }
#' @export
appendCPTAC <- function(mae, project = NULL,
    fraction = c("Proteome", "Phosphoproteome", "Glycoproteome"),
    dry_run = FALSE) {
    stopifnot(methods::is(mae, "MultiAssayExperiment"),
              is.logical(dry_run), length(dry_run) == 1L)
    fraction <- match.arg(fraction)

    project <- if (is.null(project)) .detectTcgaProject(mae)
               else paste0("TCGA-", sub("^TCGA-", "", toupper(project)))

    studies <- CPTACdata::cptacProjects
    studies <- studies[studies$tcga_project %in% project &
                       studies$analytical_fraction %in% fraction, ,
                       drop = FALSE]
    if (nrow(studies) == 0L)
        stop(sprintf("No CPTAC %s data for %s. ", fraction, project),
             "See ?cptacProjects.", call. = FALSE)

    patients <- rownames(MultiAssayExperiment::colData(mae))
    # NULL for studies without a quantification matrix: nothing to fetch
    specimens <- lapply(seq_len(nrow(studies)), function(i)
        if (!is.na(studies$data_type[i]))
            .tcgaSpecimens(studies[i, , drop = FALSE]))
    matched <- lapply(specimens, function(studySpecimens)
        studySpecimens[studySpecimens$case_submitter_id %in% patients, ,
                       drop = FALSE])

    nPatients <- function(studySpecimens)
        if (is.null(studySpecimens)) NA_integer_
        else length(unique(studySpecimens$case_submitter_id))
    if (dry_run)
        return(data.frame(
            experiment       = studies$experiment,
            pdc_study_id     = studies$pdc_study_id,
            quant_available  = !is.na(studies$data_type),
            mae_patients     = length(patients),
            pdc_patients     = vapply(specimens, nPatients, integer(1L)),
            matched_patients = vapply(matched, nPatients, integer(1L)),
            stringsAsFactors = FALSE
        ))

    for (i in seq_len(nrow(studies))) {
        row <- studies[i, , drop = FALSE]
        experimentName <- row$experiment
        if (is.na(row$data_type)) {
            warning(sprintf("%s (%s) has no quantification matrix in the ",
                            row$pdc_name, row$pdc_study_id),
                    "PDC; skipped.", call. = FALSE)
            next
        }
        if (nrow(matched[[i]]) == 0L) {
            warning(sprintf("No patients in %s match the MAE; skipped.",
                            row$pdc_study_id), call. = FALSE)
            next
        }
        if (experimentName %in% names(mae)) {
            message(sprintf("'%s' already present; skipped.", experimentName))
            next
        }
        message(sprintf("Adding %s from %s (%d patients)...", experimentName,
                        row$pdc_study_id, nPatients(matched[[i]])))
        experiment <- .buildCptacSE(row, matched[[i]])
        sampleMap <- S4Vectors::DataFrame(
            assay   = factor(rep(experimentName, ncol(experiment))),
            primary = TCGAutils::TCGAbarcode(colnames(experiment)),
            colname = colnames(experiment)
        )
        mae <- c(mae, stats::setNames(list(experiment), experimentName),
                 sampleMap = sampleMap)
    }
    mae
}


# ---------------------------------------------------------------- helpers --

.detectTcgaProject <- function(mae) {
    clinical <- MultiAssayExperiment::colData(mae)
    diseaseColumn <- grep("disease_code$", colnames(clinical), value = TRUE)[1L]
    if (is.na(diseaseColumn))
        stop("Cannot detect the TCGA project from colData; ",
             "supply it, e.g. project = \"TCGA-OV\".", call. = FALSE)
    codes <- unique(stats::na.omit(
        toupper(as.character(clinical[[diseaseColumn]]))
    ))
    if (length(codes) != 1L)
        stop("MAE spans several projects (", paste(codes, collapse = ", "),
             "); supply project explicitly.", call. = FALSE)
    paste0("TCGA-", codes)
}

# Specimens in a PDC study that are TCGA-barcoded, qualified, not pooled and
# (for studies spanning several projects) from the right disease. Non-TCGA
# IDs (e.g. "Internal Reference") are dropped before TCGAbarcode(), which
# stops on them; the patient part of the sample barcode must then agree with
# case_submitter_id.
.tcgaSpecimens <- function(row) {
    biospecimens <- .pdcBiospecimenTable(row$pdc_study_id)
    keep <- startsWith(biospecimens$sample_submitter_id, "TCGA-") &
            biospecimens$aliquot_status %in% "Qualified" &
            !(tolower(biospecimens$pool) %in% "yes")
    if (!is.na(row$pdc_disease_type))
        keep <- keep & biospecimens$disease_type %in% row$pdc_disease_type
    biospecimens <- biospecimens[keep, , drop = FALSE]
    if (nrow(biospecimens))
        biospecimens <- biospecimens[
            TCGAutils::TCGAbarcode(biospecimens$sample_submitter_id) ==
                biospecimens$case_submitter_id, , drop = FALSE]
    biospecimens
}

# Map quant matrix column labels to rows of a biospecimen table via the
# sample barcode. Labels are either the full aliquot barcode
# (phosphoproteome studies) or the sample barcode without "TCGA-" (proteome
# studies, e.g. "BH-A18U-01A"), with ".1", ".2" appended to technical
# replicates.
.matchQuantColumns <- function(labels, biospecimens) {
    baseLabels <- sub("\\.[0-9]+$", "", labels)
    isFullBarcode <- startsWith(baseLabels, "TCGA-")
    sampleBarcodes <- paste0("TCGA-", baseLabels)
    if (any(isFullBarcode))
        sampleBarcodes[isFullBarcode] <- TCGAutils::TCGAbarcode(
            baseLabels[isFullBarcode], sample = TRUE
        )
    match(sampleBarcodes, biospecimens$sample_submitter_id)
}

.buildCptacSE <- function(row, specimens) {
    quantMatrix <- .pdcQuantMatrix(row$pdc_study_id, row$data_type)
    # PDC column order is arbitrary; sorting keeps replicate ".1" over ".2"
    quantMatrix <- quantMatrix[, order(colnames(quantMatrix)), drop = FALSE]
    specimenIndex <- .matchQuantColumns(colnames(quantMatrix), specimens)
    keep <- !is.na(specimenIndex)
    keep[keep] <- !duplicated(
        specimens$sample_submitter_id[specimenIndex[keep]]
    )
    sampleData <- specimens[specimenIndex[keep],
                            c("case_submitter_id", "aliquot_submitter_id",
                              "sample_type"), drop = FALSE]
    sampleData$pdc_column <- colnames(quantMatrix)[keep]

    nMissing <- length(setdiff(specimens$case_submitter_id,
                               sampleData$case_submitter_id))
    if (nMissing)
        message(sprintf("  %d matched patient(s) not in the %s matrix.",
                        nMissing, row$pdc_study_id))

    quantMatrix <- quantMatrix[, keep, drop = FALSE]
    colnames(quantMatrix) <- specimens$sample_submitter_id[specimenIndex[keep]]
    SummarizedExperiment::SummarizedExperiment(
        assays   = stats::setNames(list(quantMatrix), row$data_type),
        colData  = S4Vectors::DataFrame(sampleData,
                                        row.names = colnames(quantMatrix)),
        metadata = list(pdc_study_id = row$pdc_study_id,
                        pdc_name = row$pdc_name, data_type = row$data_type,
                        source = .PDC_URL)
    )
}
