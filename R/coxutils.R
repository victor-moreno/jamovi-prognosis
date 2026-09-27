# Helpers for the Cox analysis. Variables are renamed internally to safe
# names ("v1_", "v2_", ...) so any jamovi variable name works in formulas;
# `lab` maps safe names back to the original names for display.

# Split a coefficient name ("v1_Female:v3_ECOG 1") into its variables and levels
coefParts <- function(nm, safe) {
    parts <- strsplit(nm, ":", fixed = TRUE)[[1]]
    v <- vapply(parts, function(p) safe[startsWith(p, safe)][1], "")
    list(vars = unname(v), levels = unname(substring(parts, nchar(v) + 1)))
}

# Level label of a coefficient: its factor levels plus the unit of any scaled
# covariate in it ("per unit" alone for an unscaled continuous term)
coefLevel <- function(levels, vars, isf, unit, . = identity) {
    u <- if (is.null(unit)) character() else setdiff(unique(unit[vars[!isf]]), "per unit")
    lev <- paste(c(levels[isf], u), collapse = " × ")
    if (nzchar(lev)) lev else .("per unit")
}

termLabel <- function(term, lab) paste(lab[strsplit(term, ":", fixed = TRUE)[[1]]], collapse = " × ")

# One row per coefficient: display labels plus estimates
coefRows <- function(fit, lab) {
    safe <- names(lab)
    b <- stats::coef(fit)
    s <- summary(fit)$coefficients
    ci <- summary(fit)$conf.int
    rows <- lapply(names(b), function(nm) {
        cp <- coefParts(nm, safe)
        isf <- cp$vars %in% names(fit$xlevels)
        lev <- paste(cp$levels[isf], collapse = " × ")
        data.frame(key = nm, term = paste(cp$vars, collapse = ":"),
                   var = paste(lab[cp$vars], collapse = " × "), level = lev,
                   hr = ci[nm, 1], lower = ci[nm, 3], upper = ci[nm, 4],
                   p = s[nm, "Pr(>|z|)"], beta = b[[nm]], se = s[nm, "se(coef)"],
                   z = s[nm, "z"], stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

# Rows for the forest plot: variable header, reference level and each level
forestRows <- function(fit, terms, lab, uni = NULL, unit = NULL, . = identity) {
    cr <- coefRows(fit, lab)
    out <- list()
    add <- function(label, r = NULL, ref = FALSE, key = NA) {
        u <- if (!is.null(uni) && !is.na(key) && key %in% uni$key) uni[uni$key == key, ] else NULL
        out[[length(out) + 1]] <<- data.frame(
            label = label, ref = ref,
            hr = if (ref) 1 else if (is.null(r)) NA else r$hr,
            lo = if (is.null(r)) NA else r$lower, hi = if (is.null(r)) NA else r$upper,
            p = if (is.null(r)) NA else r$p,
            hr_u = if (is.null(u)) NA else u$hr, lo_u = if (is.null(u)) NA else u$lower,
            hi_u = if (is.null(u)) NA else u$upper)
    }
    for (t in terms) {
        rt <- cr[cr$term == t, , drop = FALSE]
        if (!grepl(":", t, fixed = TRUE) && t %in% names(fit$xlevels)) {
            add(lab[[t]])
            add(paste0("    ", fit$xlevels[[t]][1]), ref = TRUE)
            for (i in seq_len(nrow(rt)))
                add(paste0("    ", rt$level[i]), rt[i, ], key = rt$key[i])
        } else if (!grepl(":", t, fixed = TRUE)) {
            # show the covariate unit only when HRs are scaled
            u <- if (!is.null(unit) && t %in% names(unit)) unit[[t]] else "per unit"
            add(if (u == "per unit") lab[[t]] else sprintf("%s (%s)", lab[[t]], u),
                rt[1, ], key = rt$key[1])
        } else {
            add(termLabel(t, lab))
            # braces: one forest row per interaction coefficient, not only the last
            for (i in seq_len(nrow(rt))) {
                cp <- coefParts(rt$key[i], names(lab))
                add(paste0("    ", coefLevel(cp$levels, cp$vars, cp$vars %in% names(fit$xlevels), unit, . = .)),
                    rt[i, ], key = rt$key[i])
            }
        }
    }
    do.call(rbind, out)
}

fmtHR <- function(hr, lo, hi, p, ref, . = identity) {
    ifelse(ref, .("Reference"),
    ifelse(is.na(hr), "",
           sprintf("%.2f (%.2f-%.2f)   %s", hr, lo, hi,
                   ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))))
}

forestTheme <- function() {
    theme_classic(base_size = 13) +
        theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
              axis.text.y.left = element_text(hjust = 0, colour = "black"),
              axis.text.y.right = element_text(hjust = 0, colour = "black",
                                               family = "mono", size = 11))
}

# FO4 (multivariable only) or FO5 (univariable vs multivariable)
forestPlot <- function(r, uniMulti = FALSE, col = "#1B4F8A", . = identity) {
    header <- if (uniMulti) .("HR uni | multi") else sprintf("%-19s %s", .("HR (95% CI)"), "p")
    text <- if (uniMulti)
        ifelse(r$ref, .("Reference"), ifelse(is.na(r$hr), "",
               paste0(ifelse(is.na(r$hr_u), "  -  ", sprintf("%5.2f", r$hr_u)),
                      " | ", sprintf("%.2f", r$hr))))
    else fmtHR(r$hr, r$lo, r$hi, r$p, r$ref, . = .)
    n <- nrow(r)
    y <- rev(seq_len(n))
    lab <- data.frame(y = c(n + 1, y), label = c("", r$label), text = c(header, text))

    # legend labels are translated; colours are matched by those labels
    uLab <- .("Univariable"); mLab <- .("Multivariable")
    if (uniMulti) {
        m <- data.frame(y = y - 0.15, hr = r$hr, lo = r$lo, hi = r$hi, model = mLab, ref = r$ref)
        u <- data.frame(y = y + 0.15, hr = r$hr_u, lo = r$lo_u, hi = r$hi_u, model = uLab, ref = r$ref)
        pts <- rbind(u, m)
        pts$model <- factor(pts$model, c(uLab, mLab))
        cols <- stats::setNames(c("grey55", col), c(uLab, mLab))
    } else {
        pts <- data.frame(y = y, hr = r$hr, lo = r$lo, hi = r$hi, model = mLab, ref = r$ref)
        cols <- stats::setNames(col, mLab)
    }
    refs <- data.frame(y = y[r$ref], hr = 1)
    est <- pts[!pts$ref & !is.na(pts$hr), ]

    p <- ggplot(est, aes(hr, y, colour = model)) +
        geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
        geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.8, na.rm = TRUE) +
        geom_point(shape = 15, size = 3, na.rm = TRUE) +
        geom_point(data = refs, aes(hr, y), inherit.aes = FALSE, shape = 23,
                   size = 2.5, fill = "white", colour = col) +
        scale_colour_manual(values = cols, name = NULL) +
        scale_x_log10() +
        scale_y_continuous(breaks = lab$y, labels = lab$label, limits = c(0.4, n + 1.4),
                           expand = expansion(add = 0),
                           sec.axis = dup_axis(labels = lab$text, name = NULL)) +
        labs(x = .("Hazard ratio (log scale)"), y = NULL) +
        forestTheme()
    p + if (uniMulti) theme(legend.position = "top") else theme(legend.position = "none")
}

# SC3: scaled Schoenfeld residuals per coefficient, spline smooth as cox.zph
# dashed line: the Cox estimate (constant beta under proportional hazards)
schoenfeldPlot <- function(zph, labels, coef, xlab = "Time") {
    pv <- zph$table[rownames(zph$table) != "GLOBAL", "p"]
    df <- do.call(rbind, lapply(colnames(zph$y), function(v)
        data.frame(term = sprintf("%s  (%s)", labels[[v]], fmtP(pv[[v]])),
                   time = zph$time, beta = zph$y[, v], coef = coef[[v]])))
    df$term <- factor(df$term, unique(df$term))
    ggplot(df, aes(time, beta)) +
        geom_hline(aes(yintercept = coef), linetype = 2, colour = "#C0392B") +
        geom_point(alpha = 0.35, size = 1.2) +
        geom_smooth(method = "lm", formula = y ~ splines::ns(x, 4), colour = "#1B4F8A",
                    fill = "#1B4F8A", alpha = 0.15) +
        facet_wrap(~ term, scales = "free_y", ncol = 2) +
        labs(x = xlab, y = expression(beta(t))) +
        theme_bw(base_size = 12) +
        theme(strip.background = element_blank(), strip.text = element_text(hjust = 0))
}

# AD3: direct standardisation. Curves are computed in .run (small data
# frames as plot state, not the fitted model), drawn by adjustedPlot().
adjustedCurves <- function(fit, df, var) {
    lev <- levels(df[[var]])
    curves <- do.call(rbind, lapply(lev, function(l) {
        nd <- df; nd[[var]] <- factor(l, lev)
        sf <- survival::survfit(fit, newdata = nd)
        data.frame(strata = l, time = c(0, sf$time), surv = c(1, rowMeans(as.matrix(sf$surv))))
    }))
    curves$strata <- factor(curves$strata, lev)
    km <- kmData(survival::survfit(stats::as.formula(paste("survival::Surv(time, status) ~", var)),
                                   data = df))[, c("strata", "time", "surv")]
    list(curves = curves, km = km)
}

adjustedPlot <- function(curves, km, adjusted, showKM = TRUE, pal = "jmv", xlab = "Time",
                         . = identity) {
    cols <- paletteCols(nlevels(curves$strata), pal)
    p <- ggplot(curves, aes(time, surv, colour = strata))
    if (showKM)
        p <- p + geom_step(data = km, aes(time, surv, colour = strata),
                           linetype = 2, linewidth = 0.5, alpha = 0.8)
    # parentheses are added outside .(): the extractor skips strings fully in parentheses
    p + geom_step(linewidth = 1) +
        scale_colour_manual(values = cols, name = NULL) +
        scale_y_continuous(limits = c(0, 1)) +
        labs(x = xlab, y = .("Adjusted survival probability"),
             caption = paste0(if (length(adjusted)) paste(.("Adjusted for"), paste(adjusted, collapse = ", "))
                              else .("No other variables in the model"),
                              if (showKM) paste0("   (", .("dashed: unadjusted Kaplan-Meier"), ")") else "")) +
        plotTheme() + legendInside("topright") +
        theme(plot.caption = element_text(size = 10, colour = "grey30"))
}

# Subgroup rows for a two-way interaction: the effect of the focal variable
# within each level (or quartile) of the moderator. Planned in .init (row
# labels), estimated in .run by subgroupHR().
subgroupPlan <- function(df, term, lab, unit, scale, . = identity) {
    vars <- strsplit(term, ":", fixed = TRUE)[[1]]
    if (length(vars) != 2) return(NULL)
    isf <- vapply(vars, function(v) is.factor(df[[v]]), logical(1))
    ord <- if (!isf[2] && isf[1]) rev(vars) else vars
    focal <- ord[1]; mod <- ord[2]
    modVals <- if (is.factor(df[[mod]])) levels(df[[mod]])
               # quartiles shown in original units (df holds scaled covariates);
               # unique(): tied quartiles (discrete covariates) would repeat a row
               else unique(signif(stats::quantile(df[[mod]], c(0.25, 0.5, 0.75), names = FALSE) * scale[[mod]], 3))
    focalVals <- if (is.factor(df[[focal]])) levels(df[[focal]])[-1] else NA
    g <- expand.grid(mv = seq_along(modVals), fv = seq_along(focalVals))
    data.frame(
        key = paste(term, g$fv, g$mv),
        effect = if (is.factor(df[[focal]]))
                     sprintf("%s: %s – %s", lab[[focal]], focalVals[g$fv], levels(df[[focal]])[1])
                 else sprintf("%s (%s)", lab[[focal]],
                              # "per unit" is the internal marker of an unscaled covariate
                              if (unit[[focal]] == "per unit") .("per unit") else unit[[focal]]),
        within = sprintf("%s = %s", lab[[mod]], modVals[g$mv]),
        focal = focal, mod = mod,
        fval = as.character(focalVals[g$fv]),
        mval = as.character(if (is.factor(df[[mod]])) modVals[g$mv] else modVals[g$mv] / scale[[mod]]),
        stringsAsFactors = FALSE)
}

# log HR = c'b with c the difference of model-matrix rows; other variables
# at their reference level or median
subgroupHR <- function(fit, df, plan, rhsTerms) {
    tt <- stats::delete.response(stats::terms(stats::reformulate(rhsTerms)))
    base <- df[1, , drop = FALSE]
    for (v in names(df)) {
        if (is.factor(df[[v]])) base[[v]] <- factor(levels(df[[v]])[1], levels(df[[v]]))
        else if (is.numeric(df[[v]])) base[[v]] <- stats::median(df[[v]])
    }
    b <- stats::coef(fit); b[is.na(b)] <- 0
    V <- stats::vcov(fit)
    mm <- function(nd) stats::model.matrix(tt, nd, xlev = fit$xlevels)[, names(b), drop = FALSE]
    setv <- function(nd, v, val) {
        if (is.factor(df[[v]])) nd[[v]][1] <- val else nd[[v]] <- as.numeric(val)
        nd
    }
    do.call(rbind, lapply(seq_len(nrow(plan)), function(i) {
        p <- plan[i, ]
        nd0 <- setv(base, p$mod, p$mval); nd1 <- nd0
        if (is.factor(df[[p$focal]])) nd1 <- setv(nd1, p$focal, p$fval)
        else nd1[[p$focal]] <- nd0[[p$focal]] + 1
        cv <- drop(mm(nd1) - mm(nd0))
        est <- sum(cv * b)
        se <- sqrt(drop(t(cv) %*% V %*% cv))
        data.frame(key = p$key, hr = exp(est), lower = exp(est - 1.96 * se),
                   upper = exp(est + 1.96 * se), p = 2 * stats::pnorm(-abs(est / se)))
    }))
}
