
testLabels <- c(logrank = "Log-rank", gehan = "Gehan-Breslow",
                taroneware = "Tarone-Ware", petopeto = "Peto-Peto",
                trend = "Log-rank trend")

kmClass <- if (requireNamespace('jmvcore', quietly = TRUE)) R6::R6Class(
    "kmClass",
    inherit = kmBase,
    private = list(

        # Rows are created here and only filled in .run when a table was
        # cleared (see clearWith), so plot options leave the tables untouched.
        .init = function() {
            o <- self$results
            opt <- self$options
            if (!hideUnlessReady(o, !is.null(opt$elapsed) && !is.null(opt$event) &&
                                    !noEventLevel(opt$eventLevel)))
                return()
            times <- parseTimes(opt$survTimes)
            o$survTable$setVisible(length(times) > 0)

            keys <- private$.groupKeys()
            for (k in keys) {
                o$summary$addRow(rowKey = k, values = list(group = k))
                for (t in times)
                    o$survTable$addRow(rowKey = paste(k, t), values = list(group = k, time = t))
            }
            if (length(keys) > 1)
                for (t in self$options$tests)
                    o$tests$addRow(rowKey = t, values = list(test = testLabels[[t]]))

            # taller survival plots when the number-at-risk table is drawn
            if (self$options$riskTable) {
                h <- 450 + 45 + 14 * max(1, length(keys))
                o$kmPlot$setSize(600, h)
                o$cumEventsPlot$setSize(600, h)
            }
        },

        # group levels present among usable rows ("All" without a group)
        .groupKeys = function() {
            o <- self$options
            if (is.null(o$group)) return("All")
            g <- self$data[[o$group]]
            ok <- !is.na(g)
            if (!is.null(o$elapsed)) ok <- ok & !is.na(self$data[[o$elapsed]])
            if (!is.null(o$event)) ok <- ok & !is.na(self$data[[o$event]])
            levels(droplevels(as.factor(g)[ok]))
        },

        .run = function() {
            o <- self$options
            if (is.null(o$elapsed) || is.null(o$event) || noEventLevel(o$eventLevel))
                return()

            df <- private$.cleanData()
            grouped <- !is.null(o$group)
            form <- if (grouped) survival::Surv(time, status) ~ group
                    else survival::Surv(time, status) ~ 1
            fit <- survival::survfit(form, data = df)

            tests <- NULL
            if (grouped && nlevels(df$group) > 1)
                tests <- rankTests(df$time, df$status, df$group,
                                   union("logrank", o$tests))

            if (self$results$summary$isNotFilled()) private$.fillSummary(fit)
            if (self$results$survTable$isNotFilled()) private$.fillSurvTable(fit)
            if (self$results$tests$isNotFilled()) private$.fillTests(tests)
            private$.setNotes()

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
            self$results$summary$setNote("missing",
                if (n0 > nrow(df)) sprintf("%d rows with missing values excluded", n0 - nrow(df)))
            df
        },

        # write a row by key, adding it if .init did not create it
        .putRow = function(tab, key, values) {
            values <- blankNA(values)
            if (key %in% tab$rowKeys) tab$setRow(rowKey = key, values = values)
            else tab$addRow(rowKey = key, values = values)
        },

        .fillSummary = function(fit) {
            tab <- self$results$summary
            st <- summary(fit)$table
            if (!is.matrix(st)) st <- t(as.matrix(st))
            keys <- if (is.null(fit$strata)) "All" else strataNames(names(fit$strata))
            for (i in seq_along(keys))
                private$.putRow(tab, keys[i], list(
                    group = keys[i],
                    n = st[i, "n.start"],
                    events = st[i, "events"],
                    censored = st[i, "n.start"] - st[i, "events"],
                    median = st[i, "median"],
                    mlower = st[i, "0.95LCL"],
                    mupper = st[i, "0.95UCL"]))
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
            for (i in seq_along(s$time))
                private$.putRow(tab, paste(grp[i], s$time[i]), list(
                    group = grp[i], time = s$time[i], nrisk = s$n.risk[i],
                    nevent = cumev[i], surv = s$surv[i],
                    lower = s$lower[i], upper = s$upper[i]))
            tab$setNote("beyond", if (any(s$n.risk == 0))
                "Times with 0 at risk are beyond follow-up; the last estimate is carried forward")
        },

        .fillTests = function(tests) {
            if (is.null(tests)) return()
            tab <- self$results$tests
            for (t in self$options$tests) {
                r <- tests[[t]]
                private$.putRow(tab, t, list(test = testLabels[[t]],
                    chisq = r[["chisq"]], df = r[["df"]], p = r[["p"]]))
            }
        },

        # notes follow showExplanations without refilling the tables
        .setNotes = function() {
            on <- self$options$showExplanations
            r <- self$results
            r$summary$setNote("expl", if (on)
                "Median: time at which the Kaplan-Meier survival falls to 50%; empty if not reached.")
            r$survTable$setNote("expl", if (on)
                "Survival: Kaplan-Meier probability of being event-free at that time. Events: cumulative events up to that time.")
            r$tests$setNote("expl", if (on) paste(
                "H0: the survival curves are identical.",
                "Log-rank weights all event times equally; Gehan-Breslow and Tarone-Ware",
                "give more weight to early times, Peto-Peto weights by overall survival.",
                "The trend test (1 df) uses the group order as scores."))
        },

        .plotPval = function(tests) {
            if (!self$options$pvalPlot || is.null(tests)) return(NULL)
            sel <- self$options$tests
            t <- if (length(sel)) sel[1] else "logrank"
            paste(testLabels[[t]], fmtP(tests[[t]][["p"]]))
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
                        pval = st$pval, pal = o$colours,
                        xlab = timeLabel(o$timeUnit), xmax = o$xmax, by = o$xby)
            drawGrob(g)
        })
)
