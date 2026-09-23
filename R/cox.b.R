
coxClass <- if (requireNamespace('jmvcore', quietly = TRUE)) R6::R6Class(
    "coxClass",
    inherit = coxBase,
    private = list(
        .m = NULL,

        # Table rows and plot sizes are laid out here from the data; .run only
        # fills tables that were cleared (clearWith in cox.r.yaml), so options
        # that do not affect a table leave it untouched.
        .init = function() {
            o <- self$options
            r <- self$results
            if (!hideUnlessReady(r, private$.ready())) return()
            m <- tryCatch(private$.prepare(), error = function(e) NULL)
            if (is.null(m)) return()     # the error is reported by .run
            private$.m <- m

            rows <- private$.coefLayout(m)
            if (o$uniMulti)
                for (col in c("hr", "lower", "upper", "p", "ptrend"))
                    r$coefTable$getColumn(col)$setSuperTitle("Multivariable")
            for (i in seq_len(nrow(rows)))
                r$coefTable$addRow(rowKey = rows$key[i],
                                   values = list(var = rows$var[i], level = rows$level[i]))

            if (o$globalTests)
                for (k in names(globalLabels))
                    r$globalTable$addRow(rowKey = k, values = list(test = globalLabels[[k]]))
            for (t in m$ints) {
                r$intTable$addRow(rowKey = t, values = list(term = termLabel(t, m$lab)))
                plan <- subgroupPlan(m$df, t, m$lab)
                for (i in seq_len(NROW(plan)))
                    r$subTable$addRow(rowKey = plan$key[i],
                                      values = list(effect = plan$effect[i], within = plan$within[i]))
            }
            if (o$ph) {
                for (t in m$terms)
                    r$phTable$addRow(rowKey = t, values = list(term = termLabel(t, m$lab)))
                r$phTable$addRow(rowKey = "GLOBAL", values = list(term = "Global"))
            }

            nCoef <- sum(!rows$ref)
            # forest rows: table rows plus a header per factor and per interaction
            nForest <- nrow(rows) + sum(rows$ref) + length(m$ints)
            r$forestPlot$setSize(700, 90 + 26 * nForest)
            r$phPlot$setSize(700, 60 + 240 * ceiling(nCoef / 2))
        },

        .ready = function() {
            o <- self$options
            !is.null(o$elapsed) && !is.null(o$event) && !noEventLevel(o$eventLevel) &&
                length(o$factors) + length(o$covs) > 0
        },

        .run = function() {
            o <- self$options
            r <- self$results
            if (!private$.ready()) return()

            m <- if (is.null(private$.m)) private$.prepare() else private$.m
            fit <- private$.fit(m$df, m$terms, m$strata)
            if (length(fit$coefficients) == 0 || all(is.na(stats::coef(fit))))
                jmvcore::reject("The model could not be estimated")

            uni <- NULL
            getUni <- function() {
                if (is.null(uni)) uni <<- private$.uniRows(m)
                uni
            }

            if (r$modelTable$isNotFilled()) private$.fillModel(fit)
            if (o$globalTests && r$globalTable$isNotFilled()) private$.fillGlobal(fit)
            if (r$coefTable$isNotFilled())
                private$.fillCoef(fit, m, if (o$uniMulti) getUni(),
                                  if (o$trend) private$.trendP(m))
            if (length(m$ints) && (r$intTable$isNotFilled() || r$subTable$isNotFilled()))
                private$.fillInteractions(fit, m)
            if (o$ph && r$phTable$isNotFilled()) private$.fillPH(fit)
            private$.setNotes(m)

            # plot states are small data, rebuilt only when a plot was cleared
            if (o$forest && is.null(r$forestPlot$state))
                r$forestPlot$setState(forestRows(fit, m$terms, m$lab, if (o$uniMulti) getUni()))
            if (o$phPlot && is.null(r$phPlot$state))
                private$.phState(fit, m)
            if (o$adjCurves)
                private$.adjStates(fit, m)
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
            one <- vapply(df, function(x) is.factor(x) && nlevels(x) < 2, logical(1))
            if (any(one)) {
                nm <- c(lab, if (length(strata)) stats::setNames(o$strata, strata))[names(df)[one]]
                jmvcore::reject(paste("Only one level (after removing missing values):",
                                      paste(nm, collapse = ", ")))
            }

            ints <- Filter(function(t) length(t) > 1, o$interactions)
            ints <- vapply(ints, function(t) paste(safe[match(unlist(t), vars)], collapse = ":"), "")
            list(df = df, lab = lab, main = safe, terms = c(safe, ints), ints = ints,
                 strata = strata, excluded = n0 - nrow(df))
        },

        # Rows of the hazard-ratio table, from the model matrix (same names
        # and order as the Cox coefficients): reference row + one per coefficient.
        .coefLayout = function(m) {
            cols <- colnames(stats::model.matrix(stats::reformulate(m$terms), m$df))[-1]
            safe <- names(m$lab)
            parts <- lapply(cols, coefParts, safe = safe)
            term <- vapply(parts, function(p) paste(p$vars, collapse = ":"), "")
            out <- list()
            for (t in m$terms) {
                isMainFactor <- !grepl(":", t, fixed = TRUE) && is.factor(m$df[[t]])
                label <- if (grepl(":", t, fixed = TRUE)) termLabel(t, m$lab) else m$lab[[t]]
                if (isMainFactor)
                    out[[length(out) + 1]] <- data.frame(key = paste0(t, "_ref"), term = t, var = label,
                                                         level = levels(m$df[[t]])[1], ref = TRUE)
                for (i in which(term == t)) {
                    isf <- vapply(parts[[i]]$vars, function(v) is.factor(m$df[[v]]), logical(1))
                    lev <- paste(parts[[i]]$levels[isf], collapse = " × ")
                    # a non-empty level keeps jamovi's row headers aligned
                    out[[length(out) + 1]] <- data.frame(key = cols[i], term = t, var = label,
                                                         level = if (nzchar(lev)) lev else "per unit",
                                                         ref = FALSE)
                }
            }
            do.call(rbind, out)
        },

        # small formula environment (not named .formula: that is a jmvcore
        # hook used to generate the syntax)
        .coxFormula = function(terms, strata) {
            rhs <- paste(c(terms, if (length(strata))
                                      sprintf("strata(%s)", paste(strata, collapse = ", "))),
                         collapse = " + ")
            env <- new.env(parent = baseenv())
            env$Surv <- survival::Surv
            env$strata <- survival::strata
            stats::as.formula(paste("Surv(time, status) ~", rhs), env = env)
        },

        .fit = function(df, terms, strata) {
            tryCatch(
                survival::coxph(private$.coxFormula(terms, strata), data = df, model = TRUE),
                error = function(e) jmvcore::reject(paste("Model error:", conditionMessage(e))))
        },

        # univariable fits on the same complete cases as the multivariable model
        .uniRows = function(m) {
            wald <- c()
            rows <- do.call(rbind, lapply(m$main, function(v) {
                f <- private$.fit(m$df, v, m$strata)
                wald[[v]] <<- summary(f)$waldtest[["pvalue"]]
                coefRows(f, m$lab)
            }))
            attr(rows, "wald") <- wald
            rows
        },

        # Trend tests for factors: the factor enters as.numeric() (its level
        # order as scores). Univariable: alone; multivariable: replacing only
        # that factor in the model, so the adjustment matches the HRs. Skipped
        # for factors that are part of an interaction.
        .trendP = function(m) {
            fac <- m$main[vapply(m$main, function(v) is.factor(m$df[[v]]), logical(1))]
            pOf <- function(terms, v) {
                f <- private$.fit(m$df, terms, m$strata)
                summary(f)$coefficients[sprintf("as.numeric(%s)", v), "Pr(>|z|)"]
            }
            inInt <- function(v) any(vapply(strsplit(m$ints, ":", fixed = TRUE),
                                            function(t) v %in% t, logical(1)))
            list(uni = vapply(fac, function(v) pOf(sprintf("as.numeric(%s)", v), v), numeric(1)),
                 multi = vapply(fac, function(v) if (inInt(v)) NA_real_ else
                     pOf(replace(m$terms, m$terms == v, sprintf("as.numeric(%s)", v)), v), numeric(1)))
        },

        # write a row by key (created in .init), adding it if missing
        .putRow = function(tab, key, values) {
            values <- blankNA(values)
            if (key %in% tab$rowKeys) tab$setRow(rowKey = key, values = values)
            else tab$addRow(rowKey = key, values = values)
        },

        # ---- tables ----------------------------------------------------
        .fillModel = function(fit) {
            cc <- survival::concordance(fit)
            cse <- sqrt(cc$var)
            self$results$modelTable$setRow(rowNo = 1, values = list(
                n = fit$n, events = fit$nevent, cindex = cc$concordance,
                clower = cc$concordance - 1.96 * cse, cupper = cc$concordance + 1.96 * cse))
        },

        .fillGlobal = function(fit) {
            s <- summary(fit)
            tests <- list(lr = s$logtest, wald = s$waldtest, score = s$sctest)
            for (k in names(tests))
                private$.putRow(self$results$globalTable, k, list(test = globalLabels[[k]],
                    chisq = tests[[k]][["test"]], df = tests[[k]][["df"]], p = tests[[k]][["pvalue"]]))
        },

        .fillCoef = function(fit, m, uni, trend) {
            tab <- self$results$coefTable
            cr <- coefRows(fit, m$lab)
            layout <- private$.coefLayout(m)
            blank <- list(hr_u = "", lower_u = "", upper_u = "", p_u = "", ptrend_u = "",
                          hr = "", lower = "", upper = "", p = "", ptrend = "",
                          beta = "", se = "", z = "")
            for (i in seq_len(nrow(layout))) {
                l <- layout[i, ]
                vals <- list(var = l$var, level = l$level)
                if (l$ref) {
                    # reference level: HR 1, blank CI; overall/trend p of the factor
                    vals$hr <- 1
                    if (!is.null(uni)) { vals$hr_u <- 1; vals$p_u <- attr(uni, "wald")[[l$term]] }
                    if (!is.null(trend)) { vals$ptrend_u <- trend$uni[[l$term]]
                                           vals$ptrend <- trend$multi[[l$term]] }
                } else {
                    c1 <- cr[cr$key == l$key, ]
                    vals[c("hr", "lower", "upper", "p", "beta", "se", "z")] <-
                        as.list(c1[1, c("hr", "lower", "upper", "p", "beta", "se", "z")])
                    u <- if (!is.null(uni)) uni[uni$key == l$key, ] else NULL
                    if (NROW(u)) vals[c("hr_u", "lower_u", "upper_u", "p_u")] <-
                        as.list(u[1, c("hr", "lower", "upper", "p")])
                }
                private$.putRow(tab, l$key, utils::modifyList(blank, vals))
            }
        },

        .fillInteractions = function(fit, m) {
            for (t in m$ints) {
                red <- private$.fit(m$df, setdiff(m$terms, t), m$strata)
                chi <- 2 * (fit$loglik[2] - red$loglik[2])
                df <- sum(!is.na(stats::coef(fit))) - sum(!is.na(stats::coef(red)))
                private$.putRow(self$results$intTable, t, list(term = termLabel(t, m$lab),
                    chisq = chi, df = df, p = stats::pchisq(chi, df, lower.tail = FALSE)))
                plan <- subgroupPlan(m$df, t, m$lab)
                if (is.null(plan)) next
                sg <- subgroupHR(fit, m$df, plan, m$terms)
                for (i in seq_len(nrow(sg)))
                    private$.putRow(self$results$subTable, sg$key[i], list(
                        effect = plan$effect[i], within = plan$within[i],
                        hr = sg$hr[i], lower = sg$lower[i], upper = sg$upper[i], p = sg$p[i]))
            }
        },

        .fillPH = function(fit) {
            tab <- self$results$phTable
            zph <- tryCatch(survival::cox.zph(fit), error = function(e) NULL)
            if (is.null(zph)) {
                tab$setNote("err", "The proportional hazards test could not be computed")
                return()
            }
            z <- zph$table
            for (k in rownames(z))
                private$.putRow(tab, k, list(chisq = z[k, "chisq"], df = z[k, "df"], p = z[k, "p"]))
        },

        # notes follow their options on every run without refilling tables
        .setNotes = function(m) {
            o <- self$options
            r <- self$results
            on <- o$showExplanations
            r$modelTable$setNote("missing", if (m$excluded > 0)
                sprintf("%d rows with missing values excluded", m$excluded))
            r$modelTable$setNote("strata", if (length(m$strata))
                paste("Stratified by", paste(o$strata, collapse = ", ")))
            r$modelTable$setNote("expl", if (on) paste(
                "C-index: probability that, of two patients, the one who has the event first",
                "has the higher predicted risk (0.5 = chance, 1 = perfect)."))
            r$globalTable$setNote("expl", if (on) "H0: all hazard ratios in the model are 1.")
            r$coefTable$setNote("expl", if (on) paste(
                "HR > 1: higher hazard (worse prognosis) than the reference level,",
                "or per one-unit increase of a covariate.",
                if (length(m$main) > 1)
                    "Multivariable HRs are adjusted for the other variables in the model.",
                if (o$uniMulti) "Univariable: each variable alone; on the reference row, p of the Wald test for the whole variable.",
                if (o$trend) "p trend: the factor entered as a numeric score (level order); multivariable trend keeps the other variables as in the model.",
                if (length(m$ints))
                    "With interactions, main-effect HRs apply at the reference level (or 0) of the interacting variable."))
            r$intTable$setNote("expl", if (on) "Likelihood-ratio test for adding the interaction; a small p suggests that the effect of one variable depends on the other.")
            r$subTable$setNote("expl", if (on) "Hazard ratio of the first variable within each level (or quartile) of the second; other variables at their reference level or median.")
            r$phTable$setNote("expl", if (on) "Test based on Schoenfeld residuals; a small p suggests that the hazard ratio changes over time (non-proportional hazards).")
        },

        # ---- plot states -----------------------------------------------
        .phState = function(fit, m) {
            # one panel per coefficient (the table tests whole terms)
            zph <- survival::cox.zph(fit, terms = FALSE)
            cr <- coefRows(fit, m$lab)
            labs <- stats::setNames(ifelse(nzchar(cr$level), paste0(cr$var, ": ", cr$level), cr$var),
                                    cr$key)
            self$results$phPlot$setState(list(zph = zph, labels = labs, coef = stats::coef(fit)))
        },

        # one plot per factor; curves computed here so the state stays small
        .adjStates = function(fit, m) {
            arr <- self$results$adjPlots
            for (key in arr$itemKeys) {
                img <- arr$get(key = key)
                if (!is.null(img$state)) next
                if (length(m$strata)) {
                    img$setError("Adjusted curves are not available for stratified models")
                    next
                }
                var <- names(m$lab)[m$lab == key]
                ac <- adjustedCurves(fit, m$df, var)
                img$setState(list(curves = ac$curves, km = ac$km,
                                  adjusted = setdiff(unname(m$lab), key)))
            }
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
            print(adjustedPlot(st$curves, st$km, st$adjusted, showKM = self$options$adjKM,
                               pal = self$options$colours,
                               xlab = timeLabel(self$options$timeUnit)))
            TRUE
        })
)

globalLabels <- c(lr = "Likelihood ratio", wald = "Wald", score = "Score (log-rank)")
