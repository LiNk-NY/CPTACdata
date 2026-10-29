## Build the `cptacProjects` crosswalk and save it as exported package data
## (data/cptacProjects.rda).
##
## Run from the package root: Rscript data-raw/cptacProjects.R
## Needs internet access to pdc.cancer.gov and api.gdc.cancer.gov.
##
## Only the choices below are curated by hand; study names, fractions, TCGA
## projects and disease labels come from the PDC API. The script stops if
## the API and the curated table disagree (e.g. a new TCGA study in the PDC,
## or a quantification matrix appearing for a study marked NA).
library(tibble)
library(dplyr)
devtools::load_all(quiet = TRUE)

## Curated: experiment name used in the MultiAssayExperiment, and the
## quantDataMatrix data type to download (NA: the PDC has no matrix).
curated <- tribble(
    ~tcga_project , ~pdc_study_id , ~experiment                  , ~data_type  ,
    "TCGA-BRCA"   , "PDC000173"   , "CPTAC_Proteome"             , "log2_ratio",
    "TCGA-BRCA"   , "PDC000174"   , "CPTAC_Phosphoproteome"      , "log2_ratio",
    "TCGA-COAD"   , "PDC000111"   , "CPTAC_Proteome"             , NA          ,
    "TCGA-READ"   , "PDC000111"   , "CPTAC_Proteome"             , NA          ,
    "TCGA-OV"     , "PDC000112"   , "CPTAC_JHU_Glycoproteome"    , NA          ,
    "TCGA-OV"     , "PDC000113"   , "CPTAC_JHU_Proteome"         , "log2_ratio",
    "TCGA-OV"     , "PDC000114"   , "CPTAC_PNNL_Proteome"        , "log2_ratio",
    "TCGA-OV"     , "PDC000115"   , "CPTAC_PNNL_Phosphoproteome" , "log2_ratio"
)

## PDC studies run on TCGA tumors: those in the "CPTAC-TCGA" project whose
## biospecimens carry TCGA barcodes. Study names are not reliable for this:
## the "TCGA ... CompRef Samples" studies (PDC000289-294) hold no TCGA
## patients. The PDC disease_type ("Colon Adenocarcinoma;Rectum
## Adenocarcinoma") is mapped to TCGA projects with TCGAutils::diseaseCodes.
pdcTcgaStudies <- function(project_submitter_id = "CPTAC-TCGA") {
    progs <- .pdcGQL(
        '{ allPrograms(acceptDUA: true) {
        projects { project_submitter_id studies { pdc_study_id } } } }'
    )
    projects <- unlist(
        lapply(progs$allPrograms, `[[`, "projects"),
        recursive = FALSE
    )
    projects <- Filter(
        function(p) identical(p$project_submitter_id, project_submitter_id),
        projects
    )
    ids <- sort(unique(unlist(lapply(projects, function(p) {
        vapply(p$studies, `[[`, character(1L), "pdc_study_id")
    }))))

    data("diseaseCodes", package = "TCGAutils", envir = environment())
    bind_rows(lapply(ids, function(id) {
        bio <- .pdcBiospecimenTable(id)
        if (!any(startsWith(bio$case_submitter_id, "TCGA-"))) {
            return(NULL)
        }
        s <- .pdcGQL(sprintf(
            '{ study(pdc_study_id: "%s" acceptDUA: true) {
            study_name analytical_fraction disease_type } }',
            id
        ))$study[[1L]]
        diseases <- setdiff(strsplit(s$disease_type, ";")[[1L]], "Other")
        codes <- diseaseCodes$Study.Abbreviation[
            match(tolower(diseases), tolower(diseaseCodes$Study.Name))
        ]
        if (anyNA(codes)) {
            stop(
                id,
                ": no TCGA code for disease type(s) ",
                paste(diseases[is.na(codes)], collapse = ", ")
            )
        }
        tibble(
            tcga_project = paste0("TCGA-", codes),
            pdc_study_id = id,
            pdc_name = s$study_name,
            analytical_fraction = s$analytical_fraction,
            ## only needed to split studies spanning several projects
            pdc_disease_type = if (length(diseases) > 1L) {
                diseases
            } else {
                NA_character_
            }
        )
    }))
}

gdcCaseIds <- function(project) {
    filt <- sprintf(
        '{"op":"=","content":{"field":"project.project_id","value":"%s"}}',
        project
    )
    hits <- httr2::request("https://api.gdc.cancer.gov/cases") |>
        httr2::req_url_query(
            filters = filt,
            fields = "submitter_id",
            size = 5000L,
            format = "JSON"
        ) |>
        httr2::req_retry(max_tries = 3L) |>
        httr2::req_perform() |>
        httr2::resp_body_json()
    vapply(hits$data$hits, `[[`, character(1L), "submitter_id")
}

hasMatrix <- function(id) {
    tryCatch(
        {
            .pdcQuantMatrix(id, "log2_ratio")
            TRUE
        },
        error = function(e) FALSE
    )
}

## -------------------------------------------------------------- build --

cptacProjects <- full_join(
    pdcTcgaStudies(),
    curated,
    by = c("tcga_project", "pdc_study_id")
)

drift <- filter(cptacProjects, is.na(pdc_name) | is.na(experiment))
if (nrow(drift)) {
    print(drift)
    stop("PDC studies and the curated table disagree; update `curated`.")
}

## matrix availability must match the curated data_type
avail <- vapply(unique(cptacProjects$pdc_study_id), hasMatrix, logical(1L))
wrong <- avail[cptacProjects$pdc_study_id] != !is.na(cptacProjects$data_type)
if (any(wrong)) {
    stop(
        "quantDataMatrix availability changed for: ",
        paste(unique(cptacProjects$pdc_study_id[wrong]), collapse = ", "),
        "; update `data_type` in `curated`."
    )
}

## verified: every qualified TCGA patient in the study is in the GDC project
gdc <- lapply(
    stats::setNames(nm = unique(cptacProjects$tcga_project)),
    gdcCaseIds
)
cptacProjects$verified <- vapply(
    seq_len(nrow(cptacProjects)),
    function(i) {
        row <- cptacProjects[i, ]
        cases <- unique(.tcgaSpecimens(row)$case_submitter_id)
        length(cases) > 0L && all(cases %in% gdc[[row$tcga_project]])
    },
    logical(1L)
)

cptacProjects <- cptacProjects |>
    select(
        tcga_project,
        pdc_study_id,
        pdc_name,
        experiment,
        analytical_fraction,
        data_type,
        pdc_disease_type,
        verified
    ) |>
    arrange(match(tcga_project, curated$tcga_project), pdc_study_id) |>
    as.data.frame()
print(cptacProjects)

usethis::use_data(cptacProjects, overwrite = TRUE)
