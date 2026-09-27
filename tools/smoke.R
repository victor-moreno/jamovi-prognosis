# Smoke test, run by tools/install.sh right after installing the module into
# jamovi desktop and/or the docker container (the runner has put jamovi's base
# library and the installed module first on .libPaths() and attached it).
# Reference values: the survival package on the same data.
d <- survival::lung
d <- d[!is.na(d$ph.ecog) & d$ph.ecog < 3, ]
d$event <- factor(d$status, 1:2, c("Alive", "Dead"))
d$ecog <- factor(d$ph.ecog)
d$sex <- factor(d$sex, 1:2, c("Male", "Female"))

km <- prognosis::km(data = d, elapsed = "time", event = "event", eventLevel = "Dead", group = "ecog")
lr <- km$tests$asDF$chisq[1]
ref <- survival::survdiff(survival::Surv(time, status) ~ ecog, data = d)$chisq
stopifnot(isTRUE(all.equal(lr, ref)))
cat(sprintf("   smoke test passed: log-rank chi2 %.2f\n", lr))

# ecog as an ordinal (ordered factor, as jamovi sends it): HRs per level, not .L/.Q
d$ecog <- factor(d$ecog, ordered = TRUE)
cx <- prognosis::cox(data = d, elapsed = "time", event = "event", eventLevel = "Dead",
                     factors = c("sex", "ecog"), covs = "age")
ct <- cx$coefTable$asDF
hr <- ct$hr[!is.na(ct$lower)]    # reference-level rows have HR 1 and a blank CI
fit <- survival::coxph(survival::Surv(time, status) ~ sex + factor(ph.ecog) + age, data = d)
stopifnot(isTRUE(all.equal(unname(hr), unname(exp(coef(fit))))))
cat(sprintf("   smoke test passed: %d hazard ratios match coxph\n", length(hr)))
