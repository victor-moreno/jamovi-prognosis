.libPaths(c(normalizePath(".tmp/lib"), .libPaths()))
suppressPackageStartupMessages({library(jmvcore); library(prognosis); library(survival)})
lung <- survival::lung
lung <- lung[!is.na(lung$ph.ecog) & lung$ph.ecog < 3, ]
lung$status_f <- factor(lung$status, 1:2, c("Alive", "Dead"))
lung$sex <- factor(lung$sex, 1:2, c("Male", "Female"))
lung$ecog <- factor(lung$ph.ecog, 0:2, c("ECOG 0", "ECOG 1", "ECOG 2"))
lungOrd <- lung
lungOrd$ecog <- factor(lungOrd$ecog, ordered = TRUE)

ref <- coxph(Surv(time, status == 2) ~ sex + ecog, data = lung)
r <- prognosis::cox(data = lungOrd, elapsed = "time", event = "status_f", eventLevel = "Dead",
                    factors = c("sex", "ecog"), uniMulti = TRUE, forest = TRUE)
ct <- r$coefTable$asDF
stopifnot(identical(ct$level, c("Male", "Female", "ECOG 0", "ECOG 1", "ECOG 2")))
stopifnot(isTRUE(all.equal(unname(as.numeric(ct$hr[-c(1, 3)])), unname(exp(coef(ref))))))
cat("ordinal factor: treatment contrasts, HRs match coxph\n")

r <- prognosis::cox(data = lung, elapsed = "time", event = "status_f", eventLevel = "Dead",
                    factors = c("sex", "ecog"), interactions = list(c("sex", "ecog")), forest = TRUE)
fr <- r$forestPlot$state
int <- coxph(Surv(time, status == 2) ~ sex * ecog, data = lung)
stopifnot(sum(grepl("×", fr$label) & !is.na(fr$hr)) == 2)
stopifnot(isTRUE(all.equal(fr$hr[!is.na(fr$hr) & !fr$ref], unname(exp(coef(int))))))
cat("forest plot: one row per interaction coefficient\n")
