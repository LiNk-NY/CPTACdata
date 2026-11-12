#' Get a CPTAC cohort from the PDC as a MultiAssayExperiment
#'
#' Downloads the studies of one CPTAC cohort listed in [cptacCohorts], such as
#' the CPTAC-3 clear cell renal cell carcinoma discovery cohort, and combines
#' them in a `MultiAssayExperiment`, with one experiment per study and the
#' PDC clinical data of the cohort's patients as `colData`.
#'
#' Specimens are kept if they are "Qualified", not pooled, and come from a
#' patient: reference channels such as the internal reference pool, QC
#' samples and the NCI-7 cell lines (PDC disease type "Other") are dropped.
#' Columns are named by the PDC matrix label, which for these cohorts is the
#' `aliquot_submitter_id`; technical replicates keep their `.1`, `.2`
#' suffix. The quantification is the PDC Common Data Analysis Pipeline value
#' given by the `data_type` column of [cptacCohorts].
#'
#' Biospecimen and clinical tables are cached for the R session, so a
#' `dry_run` followed by a real run downloads them only once.
#'
#' @param cohort The cohort to download, a value of the `cohort` column of
#'   [cptacCohorts], e.g. `"CPTAC CCRCC Discovery Study"`.
#' @param experiments The experiments to include, values of the `experiment`
#'   column of [cptacCohorts] for that cohort, e.g.
#'   `c("Proteome", "Phosphoproteome")`. `NULL` (default) includes every
#'   study of the cohort with a quantification matrix.
#' @param replicates How to handle technical replicates, which the PDC marks
#'   with a `.1`, `.2` suffix on the matrix column. `"all"` (default) keeps
#'   every replicate as its own column; `"first"` keeps only the first
#'   replicate of each sample.
#' @param dry_run If `TRUE`, download only sample metadata and return a
#'   one-row-per-study summary of what would be downloaded.
#'
#' @return If `dry_run = FALSE`, a `MultiAssayExperiment` with:
#'   - One experiment per study, named as in the `experiment` column of
#'     [cptacCohorts]. Each is a `SummarizedExperiment` with one assay, named
#'     after the `data_type` of the study, holding the PDC quantification
#'     matrix (missing values are `NA`). Its `colData` has
#'     `case_submitter_id`, `sample_submitter_id`, `aliquot_submitter_id`,
#'     `sample_type` and `pdc_column` (the column label in the PDC matrix),
#'     and its `metadata` has the same provenance fields as the experiments
#'     added by [appendCPTAC()].
#'   - `colData` with one row per patient (`case_submitter_id`), holding the
#'     PDC disease type, demographics (e.g. `gender`, `vital_status`,
#'     `days_to_birth`) and diagnosis (e.g. `primary_diagnosis`,
#'     `tumor_grade`, `ajcc_pathologic_stage`, `age_at_diagnosis`). Values
#'     are kept as the PDC reports them, as character. When a patient has
#'     several diagnosis records, their distinct values are joined with
#'     `"; "`.
#'   - A `sampleMap` linking every column to its patient.
#'
#'   If `dry_run = TRUE`, a `data.frame` with one row per selected study and
#'   these columns:
#'   - `experiment`: experiment name.
#'   - `pdc_study_id`: PDC study ID.
#'   - `quant_available`: `TRUE` if the PDC has a quantification matrix for
#'     the study.
#'   - `patients`: number of distinct patients in the study that pass the
#'     specimen filters.
#'   - `aliquots`: number of their aliquots.
#'
#' @seealso [cptacCohorts], [appendCPTAC()]
#'
#' @examples
#' ## Cohort and experiment names come from cptacCohorts, which ships with
#' ## the package and needs no network access
#' data("cptacCohorts", package = "CPTACdata")
#' unique(cptacCohorts[["cohort"]])
#'
#' ## Experiments of one cohort; those with data_type NA have no matrix
#' ccrccStudies <- cptacCohorts[
#'     cptacCohorts[["cohort"]] == "CPTAC CCRCC Discovery Study",
#'     c("experiment", "pdc_study_id", "data_type")
#' ]
#' ccrccStudies
#'
#' ## Downloading queries the PDC API, so it runs only after you have
#' ## accepted the PDC Data Use Agreement with
#' ## options(CPTACdata.acceptDUA = TRUE).
#' if (isTRUE(getOption("CPTACdata.acceptDUA"))) {
#' getCPTAC("CPTAC CCRCC Discovery Study", dry_run = TRUE)
#' ccrcc <- getCPTAC(
#'     "CPTAC CCRCC Discovery Study",
#'     experiments = c("Proteome", "Phosphoproteome")
#' )
#' ccrcc
#' }
#' @export
getCPTAC <- function(
    cohort,
    experiments = NULL,
    replicates = c("all", "first"),
    dry_run = FALSE
) {
    stopifnot(
        is.character(cohort),
        length(cohort) == 1L,
        is.null(experiments) || is.character(experiments),
        is.logical(dry_run),
        length(dry_run) == 1L
    )
    replicates <- match.arg(replicates)
    studies <- .selectCohortStudies(cohort, experiments)
    studyRows <- split(
        studies,
        factor(studies[["experiment"]], levels = studies[["experiment"]])
    )
    specimens <- lapply(studyRows, .cohortSpecimens)

    if (dry_run) {
        return(.cohortDryRun(studies, specimens))
    }
    available <- !is.na(studies[["data_type"]])
    experimentList <- .downloadCohortStudies(
        studyRows,
        specimens,
        available,
        replicates
    )
    .cohortMAE(experimentList, studyRows[available])
}

# One SummarizedExperiment per study with a quantification matrix; studies
# without one are skipped with a warning.
.downloadCohortStudies <- function(
    studyRows,
    specimens,
    available,
    replicates
) {
    for (row in studyRows[!available]) {
        warning(
            sprintf(
                "%s (%s) has no quantification matrix in the PDC; skipped.",
                row[["pdc_name"]],
                row[["pdc_study_id"]]
            ),
            call. = FALSE
        )
    }
    if (!any(available)) {
        stop(
            "None of the selected studies has a quantification matrix.",
            call. = FALSE
        )
    }
    experimentNames <- names(studyRows)[available]
    experimentList <- lapply(experimentNames, function(name) {
        message(sprintf(
            "Downloading %s (%s)...",
            name,
            studyRows[[name]][["pdc_study_id"]]
        ))
        .buildCptacSE(
            studyRows[[name]],
            specimens[[name]],
            replicates,
            matchColumns = .matchAliquotLabels,
            nameColumns = identity
        )
    })
    stats::setNames(experimentList, experimentNames)
}

# Rows of cptacCohorts for one cohort, optionally limited to some
# experiments.
.selectCohortStudies <- function(cohort, experiments) {
    cohorts <- .loadData("cptacCohorts")
    studies <- cohorts[cohorts[["cohort"]] %in% cohort, , drop = FALSE]
    if (!nrow(studies)) {
        stop(
            sprintf("'%s' is not a cohort in cptacCohorts. ", cohort),
            "See ?cptacCohorts.",
            call. = FALSE
        )
    }
    if (is.null(experiments)) {
        return(studies[!is.na(studies[["data_type"]]), , drop = FALSE])
    }
    unknown <- setdiff(experiments, studies[["experiment"]])
    if (length(unknown)) {
        stop(
            sprintf(
                "Not experiments of '%s': %s. Available: %s.",
                cohort,
                paste(unknown, collapse = ", "),
                paste(studies[["experiment"]], collapse = ", ")
            ),
            call. = FALSE
        )
    }
    studies[match(experiments, studies[["experiment"]]), , drop = FALSE]
}

# Specimens of a cohort study that are qualified, not pooled and from a
# patient. Reference channels (internal reference pool, QC samples, NCI-7
# cell lines) have the PDC disease type "Other".
.cohortSpecimens <- function(row) {
    biospecimens <- .pdcBiospecimenTable(row[["pdc_study_id"]])
    keep <- biospecimens[["aliquot_status"]] %in% "Qualified" &
        !(tolower(biospecimens[["pool"]]) %in% "yes") &
        !(biospecimens[["disease_type"]] %in% "Other")
    biospecimens[keep, , drop = FALSE]
}

# Map PDC matrix labels to rows of a biospecimen table: in the CPTAC
# cohorts the label is the aliquot_submitter_id, with ".1", ".2" appended to
# technical replicates.
.matchAliquotLabels <- function(labels, biospecimens) {
    match(
        sub("\\.[0-9]+$", "", labels),
        biospecimens[["aliquot_submitter_id"]]
    )
}

.cohortDryRun <- function(studies, specimens) {
    data.frame(
        experiment = studies[["experiment"]],
        pdc_study_id = studies[["pdc_study_id"]],
        quant_available = !is.na(studies[["data_type"]]),
        patients = unname(vapply(specimens, .nPatients, integer(1L))),
        aliquots = unname(vapply(specimens, nrow, integer(1L))),
        stringsAsFactors = FALSE
    )
}

# Combine the experiments of a cohort with the PDC clinical data of their
# patients.
.cohortMAE <- function(experimentList, studyRows) {
    sampleMap <- do.call(rbind, lapply(names(experimentList), function(name) {
        experiment <- experimentList[[name]]
        S4Vectors::DataFrame(
            assay = rep(name, ncol(experiment)),
            primary = experiment[["case_submitter_id"]],
            colname = colnames(experiment)
        )
    }))
    sampleMap[["assay"]] <- factor(
        sampleMap[["assay"]],
        levels = names(experimentList)
    )
    patients <- sort(unique(sampleMap[["primary"]]))
    clinical <- do.call(rbind, lapply(studyRows, function(row) {
        .pdcClinical(row[["pdc_study_id"]])
    }))
    clinical <- clinical[!duplicated(clinical[["case_submitter_id"]]), ,
        drop = FALSE]
    clinical <- clinical[
        match(patients, clinical[["case_submitter_id"]]),
        ,
        drop = FALSE
    ]
    # patients without a clinical record still get a row
    clinical[["case_submitter_id"]] <- patients
    MultiAssayExperiment::MultiAssayExperiment(
        experiments = MultiAssayExperiment::ExperimentList(experimentList),
        colData = S4Vectors::DataFrame(clinical, row.names = patients),
        sampleMap = sampleMap
    )
}
