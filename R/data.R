# Load a dataset from an installed package. The data are not lazy-loaded
# (LazyData is off), so each dataset is read explicitly with this helper.
.loadData <- function(name, package = "CPTACdata") {
    env <- new.env(parent = emptyenv())
    utils::data(list = name, package = package, envir = env)
    env[[name]]
}
