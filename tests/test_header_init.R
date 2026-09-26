.libPaths(c(normalizePath(".tmp/lib"), .libPaths()))
suppressPackageStartupMessages({library(jmvcore); library(prognosis); library(RProtoBuf)})
if (exists("initProtoBuf", envir = asNamespace("jmvcore"))) jmvcore:::initProtoBuf()
lung <- survival::lung
lung <- lung[!is.na(lung$ph.ecog) & lung$ph.ecog < 3, ]
lung$status_f <- factor(lung$status, 1:2, c("Alive", "Dead"))
lung$sex <- factor(lung$sex, 1:2, c("Male", "Female"))
lung$ecog <- factor(lung$ph.ecog, 0:2, c("ECOG 0", "ECOG 1", "ECOG 2"))
state <- NULL
step <- function(cls, optsCls, opts, label, tabs) {
  for (perform in c("init", "run")) {
    a <- cls$new(options = do.call(optsCls$new, opts), data = NULL, analysisId = 1, revision = 1)
    a$.setReadDatasetHeaderSource(function(vars) lung[0, unlist(vars), drop = FALSE])
    a$.setReadDatasetSource(function(vars) lung[, unlist(vars), drop = FALSE])
    a$init(noThrow = TRUE)
    if (!is.null(state)) a$results$fromProtoBuf(state$results, a$options$compProtoBuf(state$options), character())
    a$postInit(noThrow = TRUE)
    if (perform == "run") a$run(noThrow = TRUE)
    state <<- a$asProtoBuf()
  }
  cat("==", label, "\n")
  for (t in tabs) cat(sprintf("  %-12s rows=%d filled=%s\n", t, a$results[[t]]$rowCount, a$results[[t]]$isFilled()))
}
km <- list(elapsed = "time", event = "status_f", eventLevel = "Dead", group = "ecog",
           survTimes = "180, 365", tests = c("logrank", "gehan"))
kt <- c("summary", "survTable", "tests")
step(prognosis:::kmClass, prognosis:::kmOptions, km, "km first", kt)
step(prognosis:::kmClass, prognosis:::kmOptions, modifyList(km, list(ci = TRUE)), "km ci", kt)
step(prognosis:::kmClass, prognosis:::kmOptions, modifyList(km, list(group = "sex")), "km group sex", kt)
state <- NULL
cx <- list(elapsed = "time", event = "status_f", eventLevel = "Dead", factors = c("sex", "ecog"),
           covs = "age", globalTests = TRUE, ph = TRUE, interactions = list(c("sex", "ecog")))
ct <- c("modelTable", "globalTable", "coefTable", "intTable", "subTable", "phTable")
step(prognosis:::coxClass, prognosis:::coxOptions, cx, "cox first", ct)
step(prognosis:::coxClass, prognosis:::coxOptions, modifyList(cx, list(forest = TRUE)), "cox forest", ct)
