## Build the `cptacProjects` crosswalk and save it as exported package data
## (data/cptacProjects.rda).
##
## Run from the package root: Rscript data-raw/cptacProjects.R
## Needs internet access to pdc.cancer.gov and api.gdc.cancer.gov.
##
## Only the choices below are curated by hand; study names, fractions, TCGA
## projects and disease labels come from the PDC API. The script stops if
## the API and the curated table disagree (e.g. a new TCGA study in the PDC,
## or a quantification matrix appearing for a study marked NA), and if the
## finished table fails the checks in `checkCrosswalk()`.
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
        lapply(progs[["allPrograms"]], `[[`, "projects"),
        recursive = FALSE
    )
    projects <- Filter(
        function(p) {
            identical(p[["project_submitter_id"]], project_submitter_id)
        },
        projects
    )
    ids <- sort(unique(unlist(
        lapply(projects, function(p) {
            vapply(p[["studies"]], `[[`, character(1L), "pdc_study_id")
        })
    )))

    data("diseaseCodes", package = "TCGAutils", envir = environment())
    bind_rows(lapply(ids, function(id) {
        bio <- .pdcBiospecimenTable(id)
        if (!any(startsWith(bio[["case_submitter_id"]], "TCGA-"))) {
            return(NULL)
        }
        s <- .pdcGQL(sprintf(
            '{ study(pdc_study_id: "%s" acceptDUA: true) {
            study_name analytical_fraction disease_type } }',
            id
        ))[["study"]][[1L]]
        diseases <- setdiff(strsplit(s[["disease_type"]], ";")[[1L]], "Other")
        codes <- diseaseCodes[["Study.Abbreviation"]][
            match(tolower(diseases), tolower(diseaseCodes[["Study.Name"]]))
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
            pdc_name = s[["study_name"]],
            analytical_fraction = s[["analytical_fraction"]],
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
    vapply(hits[["data"]][["hits"]], `[[`, character(1L), "submitter_id")
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

## ------------------------------------------------------------- checks --

## Invariants of the finished crosswalk. Returns the failed checks as a
## character vector (empty if all pass); the build stops if there are any,
## so a bad table never overwrites data/cptacProjects.rda.
checkCrosswalk <- function(xw, curated) {
    failures <- character()
    check <- function(ok, message) {
        if (!isTRUE(ok)) {
            failures <<- c(failures, message)
        }
    }
    key <- function(df) {
        paste(df[["tcga_project"]], df[["pdc_study_id"]], df[["experiment"]])
    }

    check(
        setequal(key(xw), key(curated)) && !anyDuplicated(key(xw)),
        "rows do not match `curated` one-to-one"
    )
    check(
        !anyNA(unlist(xw[c("pdc_name", "experiment", "analytical_fraction")])),
        "missing pdc_name, experiment or analytical_fraction"
    )
    check(
        all(grepl("^TCGA-[A-Z]+$", xw[["tcga_project"]])),
        "malformed tcga_project"
    )
    check(
        all(grepl("^PDC[0-9]{6}$", xw[["pdc_study_id"]])),
        "malformed pdc_study_id"
    )
    check(
        all(xw[["analytical_fraction"]] %in%
            c("Proteome", "Phosphoproteome", "Glycoproteome")),
        "unexpected analytical_fraction"
    )
    check(
        !anyDuplicated(paste(xw[["tcga_project"]], xw[["experiment"]])),
        "experiment names repeat within a TCGA project"
    )
    ## a study spanning several projects must say which disease goes where
    spanning <- xw[["pdc_study_id"]] %in%
        xw[["pdc_study_id"]][duplicated(xw[["pdc_study_id"]])]
    check(
        !anyNA(xw[["pdc_disease_type"]][spanning]) &&
            !anyDuplicated(xw[["pdc_disease_type"]][spanning]),
        paste(
            "studies spanning several TCGA projects need",
            "distinct pdc_disease_type"
        )
    )
    check(
        all(xw[["verified"]]),
        "some qualified TCGA patients are not in the GDC project"
    )
    failures
}

## Unit tests for checkCrosswalk() on small synthetic tables: a valid table
## must pass and each kind of corruption must be caught.
testCheckCrosswalk <- function() {
    ov <- tibble(
        tcga_project = c("TCGA-OV", "TCGA-OV"),
        pdc_study_id = c("PDC000113", "PDC000114"),
        pdc_name = c("study A", "study B"),
        experiment = c("CPTAC_JHU_Proteome", "CPTAC_PNNL_Proteome"),
        analytical_fraction = "Proteome",
        data_type = "log2_ratio",
        pdc_disease_type = NA_character_,
        verified = TRUE
    )
    curatedOv <- ov[
        c("tcga_project", "pdc_study_id", "experiment", "data_type")
    ]

    spanning <- tibble(
        tcga_project = c("TCGA-COAD", "TCGA-READ"),
        pdc_study_id = "PDC000111",
        pdc_name = "study C",
        experiment = "CPTAC_Proteome",
        analytical_fraction = "Proteome",
        data_type = NA_character_,
        pdc_disease_type = c("Colon Adenocarcinoma", "Rectum Adenocarcinoma"),
        verified = TRUE
    )
    curatedSpanning <- spanning[
        c("tcga_project", "pdc_study_id", "experiment", "data_type")
    ]
    spanningUntyped <- mutate(spanning, pdc_disease_type = NA_character_)

    passes <- function(xw, curated) {
        length(checkCrosswalk(xw, curated)) == 0L
    }

    stopifnot(
        passes(ov, curatedOv),
        passes(spanning, curatedSpanning),
        !passes(mutate(ov, verified = c(FALSE, TRUE)), curatedOv),
        !passes(ov[1L, ], curatedOv),
        !passes(mutate(ov, pdc_name = NA_character_), curatedOv),
        !passes(mutate(ov, pdc_study_id = c("PDC1", "PDC000114")), curatedOv),
        !passes(mutate(ov, analytical_fraction = "Genome"), curatedOv),
        !passes(
            mutate(ov, experiment = "E"),
            mutate(curatedOv, experiment = "E")
        ),
        !passes(spanningUntyped, curatedSpanning)
    )
    invisible(TRUE)
}

testCheckCrosswalk()

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
avail <- vapply(unique(cptacProjects[["pdc_study_id"]]), hasMatrix, logical(1L))
wrong <- avail[cptacProjects[["pdc_study_id"]]] !=
    !is.na(cptacProjects[["data_type"]])
if (any(wrong)) {
    stop(
        "quantDataMatrix availability changed for: ",
        paste(unique(cptacProjects[["pdc_study_id"]][wrong]), collapse = ", "),
        "; update `data_type` in `curated`."
    )
}

## verified: every qualified TCGA patient in the study is in the GDC project
gdc <- lapply(
    stats::setNames(nm = unique(cptacProjects[["tcga_project"]])),
    gdcCaseIds
)
cptacProjects[["verified"]] <- vapply(
    seq_len(nrow(cptacProjects)),
    function(i) {
        row <- cptacProjects[i, ]
        cases <- unique(.tcgaSpecimens(row)[["case_submitter_id"]])
        length(cases) > 0L &&
            all(cases %in% gdc[[row[["tcga_project"]]]])
    },
    logical(1L)
)

failures <- checkCrosswalk(cptacProjects, curated)
if (length(failures)) {
    stop(paste(c("crosswalk checks failed:", failures), collapse = "\n  "))
}

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
    arrange(match(tcga_project, curated[["tcga_project"]]), pdc_study_id) |>
    as.data.frame()
print(cptacProjects)

usethis::use_data(cptacProjects, overwrite = TRUE)
