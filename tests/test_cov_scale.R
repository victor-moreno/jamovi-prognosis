.libPaths(c(normalizePath(".tmp/lib"), .libPaths()))
suppressPackageStartupMessages({library(jmvcore); library(prognosis); library(survival)})
lung <- survival::lung
lung <- lung[!is.na(lung$ph.ecog) & lung$ph.ecog < 3 & !is.na(lung$ph.karno), ]
lung$status_f <- factor(lung$status, 1:2, c("Alive", "Dead"))
lung$sex <- factor(lung$sex, 1:2, c("Male", "Female"))
run <- function(...) {
  o <- list(elapsed = "time", event = "status_f", eventLevel = "Dead", factors = "sex",
            covs = c("age", "ph.karno"), interactions = list(c("sex", "age")),
            globalTests = TRUE, ph = TRUE, forest = TRUE, uniMulti = TRUE)
  a <- prognosis:::coxClass$new(options = do.call(prognosis:::coxOptions$new, replace(o, names(list(...)), list(...))),
                                data = lung, analysisId = 1, revision = 1)
  a$init(); a$run(); a
}
ref <- coxph(Surv(time, status == 2) ~ sex + age + ph.karno + sex:age, data = lung)
b <- coef(ref)
getr <- function(a, key) lapply(a$results$coefTable$getRow(rowKey = key), function(cl) cl$value)
cmp <- function(label, got, want) cat(sprintf("%-34s got %.6f want %.6f %s\n", label, got, want,
                                             if (abs(got - want) < 1e-6) "pass" else "FAIL"))
u <- run(); s <- run(covScale = "sd"); k <- run(covScale = "custom", covMult = 10)
sdAge <- sd(lung$age); sdK <- sd(lung$ph.karno)
cmp("unit: age HR", getr(u, "v2_")$hr, exp(b[["age"]]))
cmp("sd: age HR", getr(s, "v2_")$hr, exp(b[["age"]] * sdAge))
cmp("sd: karno HR", getr(s, "v3_")$hr, exp(b[["ph.karno"]] * sdK))
cmp("sd: karno upper", getr(s, "v3_")$upper, exp((b[["ph.karno"]] + qnorm(0.975) * sqrt(vcov(ref)["ph.karno", "ph.karno"])) * sdK))
cmp("x10: karno HR", getr(k, "v3_")$hr, exp(b[["ph.karno"]] * 10))
cmp("x10: karno p same as unit", getr(k, "v3_")$p, getr(u, "v3_")$p)
uni <- coxph(Surv(time, status == 2) ~ ph.karno, data = lung)
cmp("x10: karno univariable HR", getr(k, "v3_")$hr_u, exp(coef(uni) * 10))
cat("labels:", getr(u, "v3_")$level, "|", getr(s, "v3_")$level, "|", getr(k, "v3_")$level, "\n")
for (nm in c("modelTable", "globalTable", "intTable", "phTable"))
  cat(sprintf("%-12s unchanged by scaling: %s\n", nm,
      identical(capture.output(print(u$results[[nm]])), capture.output(print(k$results[[nm]])))))
cat("sub unit:\n"); print(u$results$subTable)
cat("sub x10:\n"); print(k$results$subTable)
fr <- k$results$forestPlot$state
print(fr[, c("label", "hr")])
cat("interaction label unit:", getr(u, "v1_Female:v2_")$level, "| x10:", getr(k, "v1_Female:v2_")$level, "\n")
cmp("x10: sex x age HR", getr(k, "v1_Female:v2_")$hr, exp(b[["sexFemale:age"]] * 10))
a <- run(interactions = list(c("age", "ph.karno")), covScale = "custom", covMult = 10)
a1 <- run(interactions = list(c("age", "ph.karno")))
ik <- tail(unlist(a$results$coefTable$rowKeys), 1)
cat("cont x cont key", ik, "label unit:", getr(a1, ik)$level, "| x10:", getr(a, ik)$level, "\n")
cmp("x10: age x karno HR", getr(a, ik)$hr, getr(a1, ik)$hr^100)
print(a1$results$subTable); print(a$results$subTable)
sv <- function(x) vapply(x$results$subTable$rowKeys, function(k) x$results$subTable$getCell(rowKey = k, "hr")$value, 0)
cmp("x10: subgroup HRs = unit^10", max(abs(sv(a) - sv(a1)^10)), 0)
