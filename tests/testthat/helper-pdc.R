# Fixtures are real PDC responses captured by dev/make-fixtures.R; quant
# matrices are trimmed to their first 5 features.
gqlCalls <- new.env()

# Fixture file prefix for each GraphQL query name.
fixtureTypes <- c(
    biospecimenPerStudy = "bio",
    quantDataMatrix = "quant",
    paginatedCaseDemographicsPerStudy = "demo",
    paginatedCaseDiagnosesPerStudy = "diag"
)

# Stand-in for .pdcGQL(): replay the captured response for the study named in
# the query, or fail the way the PDC does when no matrix exists.
fakeGQL <- function(query, ...) {
    id <- regmatches(query, regexpr("PDC[0-9]{6}", query))
    queryName <- regmatches(query, regexpr("[A-Za-z]+(?=\\()", query,
                                           perl = TRUE))
    type <- fixtureTypes[[queryName]]
    key <- paste(type, id, sep = "_")
    gqlCalls[[key]] <- sum(gqlCalls[[key]], 1L)
    f <- test_path("fixtures", sprintf("gql_%s.rds", key))
    if (!file.exists(f))
        stop("PDC API error (HTTP 200): Matrix data not found!", call. = FALSE)
    readRDS(f)
}

# Stand-in for .pdcStudyCatalog(): one latest version per study, with
# made-up identifiers.
fakeCatalog <- function(latest_only = FALSE) {
    ids <- unique(c(
        .loadData("cptacProjects")[["pdc_study_id"]],
        .loadData("cptacCohorts")[["pdc_study_id"]]
    ))
    data.frame(
        pdc_study_id = ids,
        study_id = paste0("uuid-", ids),
        study_version = "2",
        is_latest_version = TRUE,
        stringsAsFactors = FALSE
    )
}

localPDC <- function(env = parent.frame()) {
    rm(list = ls(.pdcCache, all.names = TRUE), envir = .pdcCache)
    rm(list = ls(gqlCalls), envir = gqlCalls)
    local_mocked_bindings(
        .pdcGQL = fakeGQL,
        .pdcStudyCatalog = fakeCatalog,
        .env = env
    )
}

# Small MAE: `patients` plus one patient absent from the PDC.
mockMAE <- function(patients, disease_code) {
    patients <- c(patients, "TCGA-ZZ-ZZZZ")
    samples <- paste0(patients, "-01A")
    exps <- list(RNA = matrix(
        seq_along(samples), nrow = 1L, dimnames = list("g1", samples)
    ))
    cd <- S4Vectors::DataFrame(
        admin.disease_code = rep(disease_code, length(patients)),
        row.names = patients
    )
    MultiAssayExperiment::MultiAssayExperiment(
        experiments = exps, colData = cd,
        sampleMap = TCGAutils::generateMap(
            exps, cd, idConverter = TCGAutils::TCGAbarcode
        )
    )
}

brcaPatients <- c("TCGA-BH-A18U", "TCGA-C8-A131", "TCGA-AO-A12B",
                  "TCGA-BH-A0C0", "TCGA-A8-A08Z")
ovPatients <- c("TCGA-24-1555", "TCGA-25-1321", "TCGA-61-2612")
