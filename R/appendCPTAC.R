#' Append CPTAC proteomics from the PDC to a TCGA MultiAssayExperiment
#'
#' Takes a `MultiAssayExperiment` from `curatedTCGAData`, finds the CPTAC
#' studies in the Proteomic Data Commons (PDC) that were run on the same TCGA
#' tumors, and adds them as new experiments with a correct `sampleMap`.
#'
#' Only patients already in `colData(mae)` are added. Specimens are kept if
#' they carry a TCGA barcode, are "Qualified" and not pooled. Each new
#' experiment has one column per column of the PDC matrix, including every
#' technical replicate unless `replicates = "first"`. Missing values are
#' `NA`. The quantification is the PDC Common Data Analysis Pipeline value
#' given by the `data_type` column of [cptacProjects] (log2 ratios to a common
#' reference).
#'
#' Biospecimen tables are cached for the R session, so a `dry_run` followed by
#' a real run downloads them only once.
#'
#' @param mae A `MultiAssayExperiment` for a single TCGA project, such as one
#'   from `curatedTCGAData`.
#' @param project TCGA project, e.g. `"OV"` or `"TCGA-OV"`. If `NULL`, it is
#'   read from the `admin.disease_code` column of `colData(mae)`.
#' @param fraction The assay type to add: one of `"Proteome"` (default),
#'   `"Phosphoproteome"` or `"Glycoproteome"`. Call `appendCPTAC` once per
#'   type to add several.
#' @param replicates How to handle technical replicates, which the PDC marks
#'   with a `.1`, `.2` suffix on the matrix column. `"all"` (default) keeps
#'   every replicate as its own column; `"first"` keeps only the first
#'   replicate of each sample.
#' @param dry_run If `TRUE`, download only sample metadata and return a
#'   one-row-per-study summary of what would be added.
#'
#' @return If `dry_run = FALSE`, the input `mae` with one new experiment per
#'   PDC study, named as in the `experiment` column of [cptacProjects]. Each
#'   experiment is a `SummarizedExperiment` with:
#'   - One assay, named after the `data_type` of the study (e.g.
#'     `"log2_ratio"`), holding the PDC quantification matrix for the matched
#'     samples. Rows are genes or protein sites, and missing values are `NA`.
#'   - Columns named by the PDC matrix label, with `TCGA-` added where the PDC
#'     leaves it out. Proteome labels are sample barcodes (e.g.
#'     `TCGA-BH-A18U-01A`) and phosphoproteome labels are aliquot barcodes
#'     (e.g. `TCGA-24-1555-01A-32-2224-29`). The PDC marks technical
#'     replicates with a `.1`, `.2` suffix, which is kept (e.g.
#'     `TCGA-C8-A131-01A.1` and `TCGA-C8-A131-01A.2`).
#'   - `colData` with `case_submitter_id` (TCGA patient barcode),
#'     `sample_submitter_id` (TCGA sample barcode, shared by the replicates
#'     of a sample), `aliquot_submitter_id`, `sample_type`, and `pdc_column`
#'     (the column label in the PDC matrix).
#'   - `metadata` with `pdc_study_id`, `pdc_name`, `data_type`, `source` (the
#'     PDC GraphQL URL), `pdc_study_version` and `pdc_study_uuid` (the PDC
#'     version of the study and its identifier, which pins that version),
#'     `retrieved` (the time the matrix was downloaded) and
#'     `CPTACdata_version` (the package version that downloaded it). Read them
#'     with `S4Vectors::metadata(mae[["CPTAC_Proteome"]])`.
#'
#'   The new `sampleMap` links every new column to its TCGA patient. Studies
#'   that are skipped, because they have no matrix, no matching patients or
#'   are already in `mae`, add no experiment and give a warning or message.
#'
#'   If `dry_run = TRUE`, a `data.frame` with one row per study in
#'   [cptacProjects] for the requested project and fraction, and these
#'   columns:
#'   - `experiment`: experiment name the study would be added as.
#'   - `pdc_study_id`: PDC study ID.
#'   - `quant_available`: `TRUE` if the PDC has a quantification matrix for
#'     the study.
#'   - `mae_patients`: number of patients in `colData(mae)`. Same value for
#'     every row.
#'   - `pdc_patients`: number of distinct TCGA patients in the PDC study that
#'     pass the specimen filters. `NA` if there is no matrix.
#'   - `matched_patients`: number of those patients that are also in
#'     `colData(mae)`. `NA` if there is no matrix.
#'
#' @seealso [cptacProjects]
#'
#' @examples
#' ## The example queries the PDC API, so it runs only after you have accepted
#' ## the PDC Data Use Agreement with options(CPTACdata.acceptDUA = TRUE).
#' ## A small object stands in for a curatedTCGAData object here.
#' if (isTRUE(getOption("CPTACdata.acceptDUA"))) {
#' patients <- c("TCGA-BH-A18U", "TCGA-BH-A0C0", "TCGA-A8-A08Z")
#' mae <- MultiAssayExperiment::MultiAssayExperiment(
#'     experiments = list(
#'         RNA = SummarizedExperiment::SummarizedExperiment(
#'             matrix(
#'                 0, nrow = 1L, ncol = 3L,
#'                 dimnames = list("gene", patients)
#'             )
#'         )
#'     ),
#'     colData = S4Vectors::DataFrame(
#'         admin.disease_code = rep("brca", 3L), row.names = patients
#'     ),
#'     sampleMap = S4Vectors::DataFrame(
#'         assay = factor(rep("RNA", 3L)),
#'         primary = patients,
#'         colname = patients
#'     )
#' )
#' appendCPTAC(mae, fraction = "Proteome", dry_run = TRUE)
#' mae <- appendCPTAC(mae, fraction = "Proteome")
#' }
#' @export
appendCPTAC <- function(
    mae,
    project = NULL,
    fraction = c("Proteome", "Phosphoproteome", "Glycoproteome"),
    replicates = c("all", "first"),
    dry_run = FALSE
) {
    stopifnot(
        methods::is(mae, "MultiAssayExperiment"),
        is.logical(dry_run),
        length(dry_run) == 1L
    )
    fraction <- match.arg(fraction)
    replicates <- match.arg(replicates)

    project <- if (is.null(project)) {
        .detectTcgaProject(mae)
    } else {
        paste0("TCGA-", sub("^TCGA-", "", toupper(project)))
    }

    studies <- .selectStudies(project, fraction)
    patients <- rownames(MultiAssayExperiment::colData(mae))
    # One single-row table per experiment, named by experiment so that every
    # later step looks studies up by name rather than by position
    studyRows <- split(
        studies,
        factor(studies[["experiment"]], levels = studies[["experiment"]])
    )
    # NULL for studies without a quantification matrix: nothing to fetch
    specimens <- lapply(studyRows, function(row) {
        if (!is.na(row[["data_type"]])) {
            .tcgaSpecimens(row)
        }
    })
    matched <- lapply(specimens, .matchPatients, patients = patients)

    if (dry_run) {
        return(.dryRunSummary(studies, specimens, matched, patients))
    }
    for (experimentName in names(studyRows)) {
        mae <- .addStudy(
            mae,
            studyRows[[experimentName]],
            matched[[experimentName]],
            replicates
        )
    }
    mae
}

# ----------------------------------------------------------- appendCPTAC --

# Rows of cptacProjects for one TCGA project and analytical fraction.
.selectStudies <- function(project, fraction) {
    studies <- .loadData("cptacProjects")
    proj_and_frac <- studies[["tcga_project"]] %in% project &
        studies[["analytical_fraction"]] %in% fraction
    studies <- studies[proj_and_frac, , drop = FALSE]
    if (nrow(studies) == 0L) {
        stop(
            sprintf("No CPTAC %s data for %s. ", fraction, project),
            "See ?cptacProjects.",
            call. = FALSE
        )
    }
    studies
}

# The specimens of one study whose patients are in the MAE.
.matchPatients <- function(studySpecimens, patients) {
    studySpecimens[
        studySpecimens[["case_submitter_id"]] %in% patients,
        ,
        drop = FALSE
    ]
}

# Number of distinct patients, or NA when the study has no matrix.
.nPatients <- function(studySpecimens) {
    if (is.null(studySpecimens)) {
        NA_integer_
    } else {
        length(unique(studySpecimens[["case_submitter_id"]]))
    }
}

.dryRunSummary <- function(studies, specimens, matched, patients) {
    data.frame(
        experiment = studies[["experiment"]],
        pdc_study_id = studies[["pdc_study_id"]],
        quant_available = !is.na(studies[["data_type"]]),
        mae_patients = length(patients),
        pdc_patients = unname(vapply(specimens, .nPatients, integer(1L))),
        matched_patients = unname(vapply(matched, .nPatients, integer(1L))),
        stringsAsFactors = FALSE
    )
}

# Add one study to the MAE, or skip it with a warning or message.
.addStudy <- function(mae, row, matchedSpecimens, replicates) {
    experimentName <- row[["experiment"]]
    if (is.na(row[["data_type"]])) {
        warning(
            sprintf(
                "%s (%s) has no quantification matrix in the ",
                row[["pdc_name"]],
                row[["pdc_study_id"]]
            ),
            "PDC; skipped.",
            call. = FALSE
        )
        return(mae)
    }
    if (nrow(matchedSpecimens) == 0L) {
        warning(
            sprintf(
                "No patients in %s match the MAE; skipped.",
                row[["pdc_study_id"]]
            ),
            call. = FALSE
        )
        return(mae)
    }
    if (experimentName %in% names(mae)) {
        message(sprintf("'%s' already present; skipped.", experimentName))
        return(mae)
    }
    message(sprintf(
        "Adding %s from %s (%d patients)...",
        experimentName,
        row[["pdc_study_id"]],
        .nPatients(matchedSpecimens)
    ))
    experiment <- .buildCptacSE(row, matchedSpecimens, replicates)
    sampleMap <- S4Vectors::DataFrame(
        assay = factor(rep(experimentName, ncol(experiment))),
        primary = experiment[["case_submitter_id"]],
        colname = colnames(experiment)
    )
    c(
        mae,
        stats::setNames(list(experiment), experimentName),
        sampleMap = sampleMap
    )
}


# ---------------------------------------------------------------- helpers --

.detectTcgaProject <- function(mae) {
    clinical <- MultiAssayExperiment::colData(mae)
    diseaseColumn <- grep("disease_code$", colnames(clinical), value = TRUE)[1L]
    if (is.na(diseaseColumn)) {
        stop(
            "Cannot detect the TCGA project from colData; ",
            "supply it, e.g. project = \"TCGA-OV\".",
            call. = FALSE
        )
    }
    codes <- unique(stats::na.omit(
        toupper(as.character(clinical[[diseaseColumn]]))
    ))
    if (length(codes) != 1L) {
        stop(
            "MAE spans several projects (",
            paste(codes, collapse = ", "),
            "); supply project explicitly.",
            call. = FALSE
        )
    }
    paste0("TCGA-", codes)
}

# Specimens in a PDC study that are TCGA-barcoded, qualified, not pooled and
# (for studies spanning several projects) from the right disease. Non-TCGA
# IDs (e.g. "Internal Reference") are dropped before TCGAbarcode(), which
# stops on them; the patient part of the sample barcode must then agree with
# case_submitter_id.
.tcgaSpecimens <- function(row) {
    biospecimens <- .pdcBiospecimenTable(row[["pdc_study_id"]])
    keep <- startsWith(biospecimens[["sample_submitter_id"]], "TCGA-") &
        biospecimens[["aliquot_status"]] %in% "Qualified" &
        !(tolower(biospecimens[["pool"]]) %in% "yes")
    if (!is.na(row[["pdc_disease_type"]])) {
        keep <- keep &
            biospecimens[["disease_type"]] %in% row[["pdc_disease_type"]]
    }
    biospecimens <- biospecimens[keep, , drop = FALSE]
    if (nrow(biospecimens)) {
        biospecimens <- biospecimens[
            TCGAutils::TCGAbarcode(biospecimens[["sample_submitter_id"]]) ==
                biospecimens[["case_submitter_id"]],
            ,
            drop = FALSE
        ]
    }
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
    if (any(isFullBarcode)) {
        sampleBarcodes[isFullBarcode] <- TCGAutils::TCGAbarcode(
            baseLabels[isFullBarcode],
            sample = TRUE
        )
    }
    match(sampleBarcodes, biospecimens[["sample_submitter_id"]])
}

# Column names are the PDC matrix labels, which are unique and keep the
# replicate suffix. Proteome labels lack the "TCGA-" prefix, so it is added.
.labelsToColnames <- function(labels) {
    ifelse(startsWith(labels, "TCGA-"), labels, paste0("TCGA-", labels))
}

.buildCptacSE <- function(row, specimens, replicates = "all") {
    studyVersion <- .pdcStudyVersion(row[["pdc_study_id"]])
    quantMatrix <- .pdcQuantMatrix(row[["pdc_study_id"]], row[["data_type"]])
    # PDC column order is arbitrary; sorting puts replicate ".1" before ".2"
    quantMatrix <- quantMatrix[, order(colnames(quantMatrix)), drop = FALSE]
    specimenIndex <- .matchQuantColumns(colnames(quantMatrix), specimens)
    keep <- !is.na(specimenIndex)
    if (identical(replicates, "first")) {
        keep[keep] <- !duplicated(
            specimens[["sample_submitter_id"]][specimenIndex[keep]]
        )
    }
    sampleData <- specimens[
        specimenIndex[keep],
        c(
            "case_submitter_id",
            "sample_submitter_id",
            "aliquot_submitter_id",
            "sample_type"
        ),
        drop = FALSE
    ]
    sampleData[["pdc_column"]] <- colnames(quantMatrix)[keep]
    .reportUnmatched(specimens, sampleData, row[["pdc_study_id"]])

    quantMatrix <- quantMatrix[, keep, drop = FALSE]
    colnames(quantMatrix) <- .labelsToColnames(sampleData[["pdc_column"]])
    SummarizedExperiment::SummarizedExperiment(
        assays = stats::setNames(list(quantMatrix), row[["data_type"]]),
        colData = S4Vectors::DataFrame(
            sampleData,
            row.names = colnames(quantMatrix)
        ),
        metadata = .experimentMetadata(row, studyVersion)
    )
}

# Report patients in the MAE that the PDC matrix does not contain.
.reportUnmatched <- function(specimens, sampleData, pdc_study_id) {
    nMissing <- length(setdiff(
        specimens[["case_submitter_id"]],
        sampleData[["case_submitter_id"]]
    ))
    if (nMissing) {
        message(sprintf(
            "  %d matched patient(s) not in the %s matrix.",
            nMissing,
            pdc_study_id
        ))
    }
    invisible(nMissing)
}

# Provenance stored with each experiment: where and when the matrix came
# from, and which study version and package build it used.
.experimentMetadata <- function(row, studyVersion) {
    c(
        list(
            pdc_study_id = row[["pdc_study_id"]],
            pdc_name = row[["pdc_name"]],
            data_type = row[["data_type"]],
            source = .PDC_URL
        ),
        studyVersion,
        list(
            retrieved = Sys.time(),
            CPTACdata_version = as.character(
                utils::packageVersion("CPTACdata")
            )
        )
    )
}
