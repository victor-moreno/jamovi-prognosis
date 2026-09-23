
kmClass <- if (requireNamespace('jmvcore', quietly = TRUE)) R6::R6Class(
    "kmClass",
    inherit = kmBase,
    private = list(

        .init = function() {
            # taller survival plot when the number-at-risk table is drawn
            if (self$options$riskTable) {
                k <- if (is.null(self$options$group)) 1
                     else max(1, nlevels(self$data[[self$options$group]]))
                self$results$kmPlot$setSize(600, 450 + 60 + 22 * k)
                self$results$cumEventsPlot$setSize(600, 450 + 60 + 22 * k)
            }
        },

        .run = function() {
            if (is.null(self$options$elapsed) || is.null(self$options$event))
                return()

            df <- private$.cleanData()
            grouped <- !is.null(self$options$group)
            form <- if (grouped) survival::Surv(time, status) ~ group
                    else survival::Surv(time, status) ~ 1
            fit <- survival::survfit(form, data = df)

            tests <- NULL
            if (grouped && nlevels(df$group) > 1)
                tests <- rankTests(df$time, df$status, df$group,
                                   union("logrank", self$options$tests))

            private$.fillSummary(fit, df, tests)
            private$.fillSurvTable(fit)
            private$.fillTests(tests)

            # plots are refitted from the data at render time
            state <- list(df = df, pval = private$.plotPval(tests))
            for (nm in c("kmPlot", "cumEventsPlot", "cumHazPlot", "logLogPlot"))
                self$results$get(nm)$setState(state)
        },

        .cleanData = function() {
            o <- self$options
            df <- data.frame(time = timeVar(self$data[[o$elapsed]]),
                             status = eventIndicator(self$data[[o$event]], o$eventLevel))
            if (!is.null(o$group))
                df$group <- droplevels(as.factor(self$data[[o$group]]))
            n0 <- nrow(df)
            df <- df[stats::complete.cases(df), , drop = FALSE]
            if (nrow(df) == 0)
                jmvcore::reject("No complete rows (time, event and group)")
            if (n0 > nrow(df))
                self$results$summary$setNote("missing",
                    sprintf("%d rows with missing values excluded", n0 - nrow(df)))
            df
        },

        .fillSummary = function(fit, df, tests) {
            tab <- self$results$summary
            st <- summary(fit)$table
            if (!is.matrix(st)) st <- t(as.matrix(st))
            keys <- if (is.null(fit$strata)) "All" else strataNames(names(fit$strata))
            expected <- if (!is.null(tests)) attr(tests, "expected") else NULL
            for (i in seq_along(keys)) {
                tab$addRow(rowKey = keys[i], values = list(
                    group = keys[i],
                    n = st[i, "n.start"],
                    events = st[i, "events"],
                    censored = st[i, "n.start"] - st[i, "events"],
                    expected = if (is.null(expected)) NaN else expected[[keys[i]]],
                    median = st[i, "median"],
                    mlower = st[i, "0.95LCL"],
                    mupper = st[i, "0.95UCL"]))
            }
            if (self$options$showExplanations)
                tab$setNote("expl", paste(
                    "Median: time at which the Kaplan-Meier survival falls to 50%;",
                    "empty if not reached.",
                    if (!is.null(expected))
                        "Expected: events expected in each group if all groups had the same survival (log-rank)."))
        },

        .fillSurvTable = function(fit) {
            times <- parseTimes(self$options$survTimes)
            if (length(times) == 0) return()
            tab <- self$results$survTable
            s <- summary(fit, times = times, extend = TRUE)
            grp <- if (is.null(s$strata)) rep("All", length(s$time))
                   else strataNames(as.character(s$strata))
            # summary() gives events per interval; accumulate within group
            cumev <- stats::ave(s$n.event, grp, FUN = cumsum)
            for (i in seq_along(s$time)) {
                tab$addRow(rowKey = i, values = list(
                    group = grp[i], time = s$time[i], nrisk = s$n.risk[i],
                    nevent = cumev[i], surv = s$surv[i],
                    lower = s$lower[i], upper = s$upper[i]))
            }
            if (any(s$n.risk == 0))
                tab$setNote("beyond", "Times with 0 at risk are beyond follow-up; the last estimate is carried forward")
            if (self$options$showExplanations)
                tab$setNote("expl", "Survival: Kaplan-Meier probability of being event-free at that time. Events: cumulative events up to that time.")
        },

        .fillTests = function(tests) {
            sel <- self$options$tests
            if (is.null(tests) || length(sel) == 0) return()
            tab <- self$results$tests
            labels <- c(logrank = "Log-rank", gehan = "Gehan-Breslow",
                        taroneware = "Tarone-Ware", petopeto = "Peto-Peto",
                        trend = "Log-rank trend")
            for (t in sel) {
                r <- tests[[t]]
                tab$addRow(rowKey = t, values = list(test = labels[[t]],
                    chisq = r[["chisq"]], df = r[["df"]], p = r[["p"]]))
            }
            if (self$options$showExplanations)
                tab$setNote("expl", paste(
                    "H0: the survival curves are identical.",
                    "Log-rank weights all event times equally; Gehan-Breslow and Tarone-Ware",
                    "give more weight to early times, Peto-Peto weights by overall survival.",
                    "The trend test (1 df) uses the group order as scores."))
        },

        .plotPval = function(tests) {
            if (!self$options$pvalPlot || is.null(tests)) return(NULL)
            sel <- self$options$tests
            t <- if (length(sel)) sel[1] else "logrank"
            labels <- c(logrank = "Log-rank", gehan = "Gehan-Breslow",
                        taroneware = "Tarone-Ware", petopeto = "Peto-Peto",
                        trend = "Log-rank trend")
            paste(labels[[t]], fmtP(tests[[t]][["p"]]))
        },

        .plot = function(image, ggtheme, theme, ...) {
            st <- image$state
            if (is.null(st)) return(FALSE)
            o <- self$options
            df <- st$df
            form <- if ("group" %in% names(df)) survival::Surv(time, status) ~ group
                    else survival::Surv(time, status) ~ 1
            fit <- survival::survfit(form, data = df)
            fun <- switch(image$name, kmPlot = "surv", cumEventsPlot = "event",
                          cumHazPlot = "cumhaz", logLogPlot = "cloglog")
            g <- kmPlot(fit, fun = fun, ci = o$ci, censor = o$censor,
                        median = o$medianLine,
                        at = if (o$markTimes) parseTimes(o$survTimes) else numeric(),
                        risk = o$riskTable && fun %in% c("surv", "event"),
                        subtitle = st$pval, pal = o$palette,
                        xlab = timeLabel(o$timeUnit), xmax = o$xmax, by = o$xby)
            drawGrob(g)
        })
)
