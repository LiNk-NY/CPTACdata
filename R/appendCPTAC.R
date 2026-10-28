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

    xw <- CPTACdata::cptacProjects
    xw <- xw[xw$tcga_project %in% project &
             xw$analytical_fraction %in% fraction, , drop = FALSE]
    if (nrow(xw) == 0L)
        stop(sprintf("No CPTAC %s data for %s. ", fraction, project),
             "See ?cptacProjects.", call. = FALSE)

    patients <- rownames(MultiAssayExperiment::colData(mae))
    # NULL for studies without a quantification matrix: nothing to fetch
    specimens <- lapply(seq_len(nrow(xw)), function(i)
        if (!is.na(xw$data_type[i])) .tcgaSpecimens(xw[i, , drop = FALSE]))
    matched <- lapply(specimens,
                      function(s) s[s$case_submitter_id %in% patients, ,
                                    drop = FALSE])

    nPatients <- function(s)
        if (is.null(s)) NA_integer_ else length(unique(s$case_submitter_id))
    if (dry_run)
        return(data.frame(
            experiment       = xw$experiment,
            pdc_study_id     = xw$pdc_study_id,
            quant_available  = !is.na(xw$data_type),
            mae_patients     = length(patients),
            pdc_patients     = vapply(specimens, nPatients, integer(1L)),
            matched_patients = vapply(matched, nPatients, integer(1L)),
            stringsAsFactors = FALSE
        ))

    for (i in seq_len(nrow(xw))) {
        row <- xw[i, , drop = FALSE]
        nm  <- row$experiment
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
        if (nm %in% names(mae)) {
            message(sprintf("'%s' already present; skipped.", nm))
            next
        }
        message(sprintf("Adding %s from %s (%d patients)...", nm,
                        row$pdc_study_id, nPatients(matched[[i]])))
        se <- .buildCptacSE(row, matched[[i]])
        smap <- S4Vectors::DataFrame(
            assay   = factor(rep(nm, ncol(se))),
            primary = TCGAutils::TCGAbarcode(colnames(se)),
            colname = colnames(se)
        )
        mae <- c(mae, stats::setNames(list(se), nm), sampleMap = smap)
    }
    mae
}


# ---------------------------------------------------------------- helpers --

.detectTcgaProject <- function(mae) {
    cd  <- MultiAssayExperiment::colData(mae)
    col <- grep("disease_code$", colnames(cd), value = TRUE)[1L]
    if (is.na(col))
        stop("Cannot detect the TCGA project from colData; ",
             "supply it, e.g. project = \"TCGA-OV\".", call. = FALSE)
    codes <- unique(stats::na.omit(toupper(as.character(cd[[col]]))))
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
    bio <- .pdcBiospecimenTable(row$pdc_study_id)
    ok  <- startsWith(bio$sample_submitter_id, "TCGA-") &
           bio$aliquot_status %in% "Qualified" &
           !(tolower(bio$pool) %in% "yes")
    if (!is.na(row$pdc_disease_type))
        ok <- ok & bio$disease_type %in% row$pdc_disease_type
    bio <- bio[ok, , drop = FALSE]
    if (nrow(bio))
        bio <- bio[TCGAutils::TCGAbarcode(bio$sample_submitter_id) ==
                   bio$case_submitter_id, , drop = FALSE]
    bio
}

# Map quant matrix column labels to rows of a biospecimen table via the
# sample barcode. Labels are either the full aliquot barcode
# (phosphoproteome studies) or the sample barcode without "TCGA-" (proteome
# studies, e.g. "BH-A18U-01A"), with ".1", ".2" appended to technical
# replicates.
.matchQuantColumns <- function(labels, bio) {
    lab <- sub("\\.[0-9]+$", "", labels)
    full <- startsWith(lab, "TCGA-")
    sample <- paste0("TCGA-", lab)
    if (any(full))
        sample[full] <- TCGAutils::TCGAbarcode(lab[full], sample = TRUE)
    match(sample, bio$sample_submitter_id)
}

.buildCptacSE <- function(row, specimens) {
    mat <- .pdcQuantMatrix(row$pdc_study_id, row$data_type)
    # PDC column order is arbitrary; sorting keeps replicate ".1" over ".2"
    mat <- mat[, order(colnames(mat)), drop = FALSE]
    idx <- .matchQuantColumns(colnames(mat), specimens)
    keep <- !is.na(idx)
    keep[keep] <- !duplicated(specimens$sample_submitter_id[idx[keep]])
    cd <- specimens[idx[keep], c("case_submitter_id", "aliquot_submitter_id",
                                 "sample_type"), drop = FALSE]
    cd$pdc_column <- colnames(mat)[keep]

    nmiss <- length(setdiff(specimens$case_submitter_id, cd$case_submitter_id))
    if (nmiss)
        message(sprintf("  %d matched patient(s) not in the %s matrix.",
                        nmiss, row$pdc_study_id))

    mat <- mat[, keep, drop = FALSE]
    colnames(mat) <- specimens$sample_submitter_id[idx[keep]]
    SummarizedExperiment::SummarizedExperiment(
        assays   = stats::setNames(list(mat), row$data_type),
        colData  = S4Vectors::DataFrame(cd, row.names = colnames(mat)),
        metadata = list(pdc_study_id = row$pdc_study_id,
                        pdc_name = row$pdc_name, data_type = row$data_type,
                        source = .PDC_URL)
    )
}
