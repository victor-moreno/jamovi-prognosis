.libPaths(c(normalizePath(".tmp/lib"), .libPaths()))
suppressPackageStartupMessages({library(jmvcore); library(prognosis)})
coxClass <- prognosis:::coxClass; coxOptions <- prognosis:::coxOptions
kmClass <- prognosis:::kmClass; kmOptions <- prognosis:::kmOptions
lung <- survival::lung
lung <- lung[!is.na(lung$ph.ecog) & lung$ph.ecog < 3, ]
lung$time_m <- lung$time / 30.4375
lung$status_f <- factor(lung$status, 1:2, c("Alive", "Dead"))
lung$sex <- factor(lung$sex, 1:2, c("Male", "Female"))
lung$ecog <- factor(lung$ph.ecog, 0:2, c("ECOG 0", "ECOG 1", "ECOG 2"))
suppressPackageStartupMessages(library(RProtoBuf))
if (exists("initProtoBuf", envir = asNamespace("jmvcore"))) jmvcore:::initProtoBuf()
base <- list(elapsed = "time_m", event = "status_f", eventLevel = "Dead",
             factors = c("sex", "ecog"), covs = "age", globalTests = TRUE, ph = TRUE,
             interactions = list(c("sex", "ecog")))
restore <- function(cls, optsCls, o1, change, label, tables) {
  a1 <- cls$new(options = do.call(optsCls$new, o1), data = lung, analysisId = 1, revision = 1); a1$init(); a1$run()
  pb <- a1$asProtoBuf()
  o2 <- utils::modifyList(o1, change)
  a2 <- cls$new(options = do.call(optsCls$new, o2), data = lung, analysisId = 1, revision = 2); a2$init()
  a2$results$fromProtoBuf(pb$results, oChanges = names(change), vChanges = character())
  st <- vapply(tables, function(t) if (a2$results[[t]]$isNotFilled()) "REFILL" else "kept", "")
  a2$run()
  ok <- vapply(tables, function(t) a2$results[[t]]$isFilled(), logical(1))
  cat(sprintf("%-26s %s | after run filled: %s\n", label,
              paste(names(st), st, sep = "=", collapse = " "), all(ok)))
}
ct <- c("modelTable", "globalTable", "coefTable", "intTable", "subTable", "phTable")
restore(coxClass, coxOptions, base, list(forest = TRUE), "cox: forest on", ct)
restore(coxClass, coxOptions, base, list(colours = "set1"), "cox: colours", ct)
restore(coxClass, coxOptions, base, list(showExplanations = TRUE), "cox: explanations", ct)
restore(coxClass, coxOptions, base, list(uniMulti = TRUE), "cox: uniMulti (clears HR)", ct)
restore(coxClass, coxOptions, base, list(coefDetails = TRUE), "cox: beta/SE/z", ct)
restore(coxClass, coxOptions, base, list(covScales = list(list(var = "age", scale = "sd"))), "cox: covScales (clears HR)", ct)
restore(coxClass, coxOptions, base, list(refLevels = list(list(var = "sex", ref = "Female"))), "cox: refLevels (clears HR)", ct)
kb <- list(elapsed = "time_m", event = "status_f", eventLevel = "Dead", group = "ecog",
           survTimes = "6, 12, 60", tests = c("logrank", "gehan"))
kt <- c("summary", "survTable", "tests")
restore(kmClass, kmOptions, kb, list(ci = TRUE), "km: CI on", kt)
restore(kmClass, kmOptions, kb, list(riskTable = TRUE), "km: risk table", kt)
restore(kmClass, kmOptions, kb, list(survTimes = "12"), "km: times (clears surv)", kt)
restore(kmClass, kmOptions, kb, list(tests = c("logrank")), "km: tests (clears tests)", kt)
