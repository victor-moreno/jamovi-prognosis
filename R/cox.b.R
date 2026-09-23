
coxClass <- if (requireNamespace('jmvcore', quietly = TRUE)) R6::R6Class(
    "coxClass",
    inherit = coxBase,
    private = list(

        .run = function() {
            o <- self$options
            if (is.null(o$elapsed) || is.null(o$event) ||
                length(o$factors) + length(o$covs) == 0)
                return()

            m <- private$.prepare()
            fit <- private$.fit(m$df, m$terms, m$strata)
            if (length(fit$coefficients) == 0 || all(is.na(stats::coef(fit))))
                jmvcore::reject("The model could not be estimated")

            uni <- if (o$uniMulti) private$.uniRows(m) else NULL

            private$.fillModel(fit, m)
            private$.fillGlobal(fit)
            private$.fillCoef(fit, m, uni)
            private$.fillInteractions(fit, m)
            private$.fillPH(fit, m)

            if (o$forest) {
                r <- forestRows(fit, m$terms, m$lab, uni)
                self$results$forestPlot$setSize(700, 90 + 26 * nrow(r))
                self$results$forestPlot$setState(r)
            }
            if (o$phPlot) {
                # one panel per coefficient (the table above tests whole terms)
                zph <- survival::cox.zph(fit, terms = FALSE)
                cr <- coefRows(fit, m$lab)
                labs <- stats::setNames(ifelse(nzchar(cr$level), paste0(cr$var, ": ", cr$level), cr$var),
                                        cr$key)
                self$results$phPlot$setSize(700, 60 + 240 * ceiling(ncol(zph$y) / 2))
                self$results$phPlot$setState(list(zph = zph, labels = labs, coef = stats::coef(fit)))
            }
            private$.prepareAdjusted(fit, m)
        },

        # ---- data and model --------------------------------------------
        .prepare = function() {
            o <- self$options
            vars <- c(o$factors, o$covs)
            safe <- sprintf("v%d_", seq_along(vars))
            lab <- stats::setNames(vars, safe)
            df <- data.frame(time = timeVar(self$data[[o$elapsed]]),
                             status = eventIndicator(self$data[[o$event]], o$eventLevel))
            for (i in seq_along(vars)) {
                x <- self$data[[vars[i]]]
                df[[safe[i]]] <- if (vars[i] %in% o$factors) as.factor(x) else jmvcore::toNumeric(x)
            }
            strata <- character()
            for (i in seq_along(o$strata)) {
                strata[i] <- sprintf("s%d_", i)
                df[[strata[i]]] <- as.factor(self$data[[o$strata[i]]])
            }
            n0 <- nrow(df)
            df <- df[stats::complete.cases(df), , drop = FALSE]
            for (v in names(df)) if (is.factor(df[[v]])) df[[v]] <- droplevels(df[[v]])
            if (nrow(df) == 0) jmvcore::reject("No complete rows")
            if (sum(df$status) == 0) jmvcore::reject("There are no events")

            ints <- Filter(function(t) length(t) > 1, o$interactions)
            ints <- vapply(ints, function(t) paste(safe[match(unlist(t), vars)], collapse = ":"), "")
            list(df = df, lab = lab, main = safe, terms = c(safe, ints), ints = ints,
                 strata = strata, excluded = n0 - nrow(df))
        },

        # small formula environment: the fitted model is kept as plot state
        .formula = function(terms, strata) {
            rhs <- paste(c(terms, if (length(strata))
                                      sprintf("strata(%s)", paste(strata, collapse = ", "))),
                         collapse = " + ")
            env <- new.env(parent = baseenv())
            env$Surv <- survival::Surv
            env$strata <- survival::strata
            stats::as.formula(paste("Surv(time, status) ~", rhs), env = env)
        },

        .fit = function(df, terms, strata) {
            fit <- tryCatch(
                survival::coxph(private$.formula(terms, strata), data = df, model = TRUE),
                error = function(e) jmvcore::reject(paste("Model error:", conditionMessage(e))))
            fit
        },

        # univariable fits on the same complete cases as the multivariable model
        .uniRows = function(m) {
            do.call(rbind, lapply(m$main, function(v)
                coefRows(private$.fit(m$df, v, m$strata), m$lab)))
        },

        # ---- tables ----------------------------------------------------
        .fillModel = function(fit, m) {
            o <- self$options
            tab <- self$results$modelTable
            cc <- survival::concordance(fit)
            cse <- sqrt(cc$var)
            row <- list(n = fit$n, events = fit$nevent, cindex = cc$concordance,
                        clower = cc$concordance - 1.96 * cse,
                        cupper = cc$concordance + 1.96 * cse)
            if (o$cBoot) row$cboot <- private$.cBoot(fit, m, cc$concordance)
            tab$setRow(rowNo = 1, values = row)
            if (m$excluded > 0)
                tab$setNote("missing", sprintf("%d rows with missing values excluded", m$excluded))
            if (length(m$strata))
                tab$setNote("strata", paste("Stratified by", paste(self$options$strata, collapse = ", ")))
            if (o$showExplanations)
                tab$setNote("expl", paste(
                    "C-index: probability that, of two patients, the one who has the event first",
                    "has the higher predicted risk (0.5 = chance, 1 = perfect).",
                    if (o$cBoot) "Optimism-corrected: bootstrap estimate of the C-index expected in new patients."))
        },

        # Harrell's bootstrap optimism: C(boot model, boot data) - C(boot model, original data)
        .cBoot = function(fit, m, capp) {
            set.seed(1234)
            f <- private$.formula(m$terms, m$strata)
            df <- m$df
            tt <- stats::delete.response(stats::terms(stats::reformulate(m$terms)))
            cf <- stats::as.formula(paste("survival::Surv(time, status) ~ lp",
                                          if (length(m$strata))
                                              sprintf("+ survival::strata(%s)", paste(m$strata, collapse = ", "))
                                          else ""))
            opt <- vapply(seq_len(self$options$bootN), function(i) {
                db <- df[sample.int(nrow(df), replace = TRUE), , drop = FALSE]
                fb <- tryCatch(suppressWarnings(survival::coxph(f, data = db, model = TRUE)), error = function(e) NULL)
                if (is.null(fb)) return(NA_real_)
                # X b by hand: predict.coxph fails on new data for stratified models
                b <- stats::coef(fb); b[is.na(b)] <- 0
                lp <- tryCatch(drop(stats::model.matrix(tt, df, xlev = fb$xlevels)[, names(b), drop = FALSE] %*% b),
                               error = function(e) NULL)
                if (is.null(lp)) return(NA_real_)
                d2 <- df; d2$lp <- lp
                # warns on rebuilding the model frame of stratified fits; harmless
                suppressWarnings(survival::concordance(fb))$concordance -
                    survival::concordance(cf, data = d2, reverse = TRUE)$concordance
            }, numeric(1))
            capp - mean(opt, na.rm = TRUE)
        },

        .fillGlobal = function(fit) {
            if (!self$options$globalTests) return()
            tab <- self$results$globalTable
            s <- summary(fit)
            tests <- list(lr = list("Likelihood ratio", s$logtest),
                          wald = list("Wald", s$waldtest),
                          score = list("Score (log-rank)", s$sctest))
            for (k in names(tests))
                tab$addRow(rowKey = k, values = list(test = tests[[k]][[1]],
                    chisq = tests[[k]][[2]][["test"]], df = tests[[k]][[2]][["df"]],
                    p = tests[[k]][[2]][["pvalue"]]))
            if (self$options$showExplanations)
                tab$setNote("expl", "H0: all hazard ratios in the model are 1.")
        },

        .fillCoef = function(fit, m, uni) {
            o <- self$options
            tab <- self$results$coefTable
            cr <- coefRows(fit, m$lab)
            if (o$uniMulti) {
                for (col in c("hr", "lower", "upper", "p"))
                    tab$getColumn(col)$setSuperTitle("Multivariable")
            }
            for (i in seq_len(nrow(cr))) {
                r <- cr[i, ]
                u <- if (!is.null(uni) && r$key %in% uni$key) uni[uni$key == r$key, ] else NULL
                tab$addRow(rowKey = r$key, values = list(
                    var = r$var, level = r$level,
                    hr_u = if (is.null(u)) "" else u$hr,
                    lower_u = if (is.null(u)) "" else u$lower,
                    upper_u = if (is.null(u)) "" else u$upper,
                    p_u = if (is.null(u)) "" else u$p,
                    hr = r$hr, lower = r$lower, upper = r$upper, p = r$p,
                    beta = r$beta, se = r$se, z = r$z))
            }
            if (o$showExplanations)
                tab$setNote("expl", paste(
                    "HR > 1: higher hazard (worse prognosis) than the reference level,",
                    "or per one-unit increase of a covariate.",
                    if (length(c(o$factors, o$covs)) > 1)
                        "Multivariable HRs are adjusted for the other variables in the model.",
                    if (o$uniMulti) "Univariable: each variable alone.",
                    if (length(m$ints))
                        "With interactions, main-effect HRs apply at the reference level (or 0) of the interacting variable."))
        },

        .fillInteractions = function(fit, m) {
            if (length(m$ints) == 0) return()
            it <- self$results$intTable
            st <- self$results$subTable
            for (t in m$ints) {
                red <- private$.fit(m$df, setdiff(m$terms, t), m$strata)
                chi <- 2 * (fit$loglik[2] - red$loglik[2])
                df <- sum(!is.na(stats::coef(fit))) - sum(!is.na(stats::coef(red)))
                it$addRow(rowKey = t, values = list(term = termLabel(t, m$lab), chisq = chi,
                    df = df, p = stats::pchisq(chi, df, lower.tail = FALSE)))
                sg <- subgroupHR(fit, m$df, t, m$lab, m$terms)
                for (i in seq_len(NROW(sg)))
                    st$addRow(rowKey = paste(t, i), values = as.list(sg[i, ]))
            }
            if (self$options$showExplanations) {
                it$setNote("expl", "Likelihood-ratio test for adding the interaction; a small p suggests that the effect of one variable depends on the other.")
                st$setNote("expl", "Hazard ratio of the first variable within each level (or quartile) of the second; other variables at their reference level or median.")
            }
        },

        .fillPH = function(fit, m) {
            if (!self$options$ph) return(NULL)
            tab <- self$results$phTable
            zph <- tryCatch(survival::cox.zph(fit), error = function(e) NULL)
            if (is.null(zph)) {
                tab$setNote("err", "The proportional hazards test could not be computed")
                return(NULL)
            }
            z <- zph$table
            for (k in rownames(z))
                tab$addRow(rowKey = k, values = list(
                    term = if (k == "GLOBAL") "Global" else termLabel(k, m$lab),
                    chisq = z[k, "chisq"], df = z[k, "df"], p = z[k, "p"]))
            if (self$options$showExplanations)
                tab$setNote("expl", "Test based on Schoenfeld residuals; a small p suggests that the hazard ratio changes over time (non-proportional hazards).")
            zph
        },

        .prepareAdjusted = function(fit, m) {
            o <- self$options
            if (is.null(o$adjVar)) return()
            img <- self$results$adjPlot
            if (!o$adjVar %in% o$factors)
                return(img$setError("The variable for adjusted curves must also be one of the Factors"))
            if (length(m$strata))
                return(img$setError("Adjusted curves are not available for stratified models"))
            var <- names(m$lab)[m$lab == o$adjVar]
            img$setState(list(fit = fit, df = m$df, var = var,
                              adjusted = setdiff(unname(m$lab), o$adjVar)))
        },

        # ---- plots -----------------------------------------------------
        .forestPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(forestPlot(image$state, uniMulti = self$options$uniMulti))
            TRUE
        },

        .phPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state)) return(FALSE)
            print(schoenfeldPlot(image$state$zph, image$state$labels, image$state$coef,
                                 xlab = timeLabel(self$options$timeUnit)))
            TRUE
        },

        .adjPlot = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            print(adjustedPlot(st$fit, st$df, st$var, self$options$adjVar, st$adjusted,
                               km = self$options$adjKM, pal = self$options$palette,
                               xlab = timeLabel(self$options$timeUnit)))
            TRUE
        })
)
