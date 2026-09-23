# Shared helpers for the prognosis analyses (km, cox).
# Plotting uses only ggplot2 + gtable (bundled with jamovi) so the module
# carries no extra dependencies.

# ---- data ---------------------------------------------------------------

# Event indicator: factor with a chosen event level, or numeric 0/1 or 1/2
# (the codings accepted by survival::Surv). No silent guessing of the level.
eventIndicator <- function(x, level) {
    if (is.factor(x) || is.character(x)) {
        if (is.null(level) || !nzchar(level))
            jmvcore::reject("Choose the level of the event variable that indicates the event")
        return(as.integer(as.character(x) == level))
    }
    x <- jmvcore::toNumeric(x)
    v <- unique(x[!is.na(x)])
    if (all(v %in% c(0, 1))) return(as.integer(x))
    if (all(v %in% c(1, 2))) return(as.integer(x == 2))
    jmvcore::reject("A numeric event variable must be coded 0/1 or 1/2; otherwise set it as nominal and choose the event level")
}

timeVar <- function(x) {
    x <- jmvcore::toNumeric(x)
    if (any(x < 0, na.rm = TRUE))
        jmvcore::reject("The time variable contains negative values")
    x
}

timeLabel <- function(unit) if (unit == "none") "Time" else sprintf("Time (%s)", unit)

parseTimes <- function(s) {
    if (is.null(s) || !nzchar(trimws(s))) return(numeric())
    v <- suppressWarnings(as.numeric(strsplit(s, "[,; ]+")[[1]]))
    sort(unique(v[!is.na(v) & v >= 0]))
}

fmtP <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)

# ---- palettes -----------------------------------------------------------

paletteCols <- function(n, pal = "jmv") {
    switch(pal,
        set1  = rep_len(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                          "#A65628", "#F781BF", "#999999"), n),
        okabe = rep_len(c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2",
                          "#D55E00", "#F0E442", "#000000"), n),
        rep_len(jmvcore::colorPalette(max(n, 1), "jmv", "color"), n))
}

# ---- weighted k-sample log-rank tests -----------------------------------
# Implemented here (instead of coin/EnvStats) to avoid dependencies.
# Weights: log-rank 1, Gehan-Breslow n(t), Tarone-Ware sqrt(n(t)),
# Peto-Peto pooled KM S(t-) (as survival::survdiff with rho = 1).
rankTests <- function(time, status, group, weights = "logrank") {
    group <- droplevels(as.factor(group))
    lev <- levels(group)
    K <- length(lev)
    tt <- sort(unique(time[status == 1]))
    nj <- sapply(lev, function(l) {
        tl <- sort(time[group == l])
        length(tl) - findInterval(tt, tl, left.open = TRUE)
    })
    dj <- sapply(lev, function(l) tabulate(match(time[group == l & status == 1], tt), length(tt)))
    nj <- matrix(nj, ncol = K); dj <- matrix(dj, ncol = K)
    n <- rowSums(nj); d <- rowSums(dj)
    expected <- colSums(nj * d / n)
    vt <- ifelse(n > 1, d * (n - d) / (n - 1), 0)
    s_left <- c(1, cumprod(1 - d / n))[seq_along(tt)]

    one <- function(w) {
        U <- colSums(w * (dj - nj * d / n))
        V <- matrix(0, K, K)
        for (i in seq_along(tt)) {
            p <- nj[i, ] / n[i]
            V <- V + w[i]^2 * vt[i] * (diag(p, K) - tcrossprod(p))
        }
        list(U = U, V = V)
    }
    out <- list()
    for (wt in weights) {
        w <- switch(wt, logrank = rep(1, length(tt)), trend = rep(1, length(tt)),
                    gehan = n, taroneware = sqrt(n), petopeto = s_left)
        uv <- one(w)
        if (wt == "trend") {
            sc <- seq_len(K)
            chi <- sum(sc * uv$U)^2 / drop(t(sc) %*% uv$V %*% sc)
            df <- 1
        } else {
            idx <- seq_len(K - 1)
            chi <- drop(t(uv$U[idx]) %*% solve(uv$V[idx, idx, drop = FALSE], uv$U[idx]))
            df <- K - 1
        }
        out[[wt]] <- c(chisq = chi, df = df,
                       p = stats::pchisq(chi, df, lower.tail = FALSE))
    }
    attr(out, "expected") <- stats::setNames(expected, lev)
    out
}

# ---- Kaplan-Meier plot data ---------------------------------------------

strataNames <- function(x) sub("^[^=]*=", "", x)

kmData <- function(fit) {
    st <- if (is.null(fit$strata)) rep("All", length(fit$time))
          else rep(strataNames(names(fit$strata)), fit$strata)
    lev <- unique(st)
    df <- data.frame(strata = st, time = fit$time, surv = fit$surv,
                     lower = fit$lower, upper = fit$upper,
                     n.censor = fit$n.censor, cumhaz = fit$cumhaz)
    df0 <- data.frame(strata = lev, time = 0, surv = 1, lower = 1, upper = 1,
                      n.censor = 0, cumhaz = 0)
    df <- rbind(df0, df)
    df$strata <- factor(df$strata, lev)
    df[order(df$strata, df$time), ]
}

# geom_ribbon interpolates linearly; expand to a step shape for KM CIs
stepRibbon <- function(df) {
    do.call(rbind, lapply(split(df, df$strata), function(g) {
        n <- nrow(g)
        i <- rep(seq_len(n), each = 2)[-1]
        j <- rep(seq_len(n), each = 2)[-2 * n]
        data.frame(strata = g$strata[1], time = g$time[i],
                   ylo = g$ylo[j], yhi = g$yhi[j])
    }))
}

xBreaks <- function(xmax, by) {
    if (by > 0) return(seq(0, xmax, by = by))
    b <- pretty(c(0, xmax))
    b[b <= xmax]
}

# Legend inside the panel without a box (the KM5 look); corner depends on
# whether the curves go down (survival) or up (incidence, hazard).
legendInside <- function(corner = c("topright", "topleft")) {
    corner <- match.arg(corner)
    pos <- if (corner == "topright") c(0.98, 0.98) else c(0.02, 0.98)
    theme(legend.position = "inside", legend.position.inside = pos,
          legend.justification = c(if (corner == "topright") 1 else 0, 1),
          legend.background = element_blank(), legend.key = element_blank(),
          legend.title = element_blank(), legend.key.width = grid::unit(1.6, "lines"))
}

plotTheme <- function() {
    theme_classic(base_size = 13) +
        theme(plot.subtitle = element_text(size = 11, colour = "grey25"))
}

# Survival-type plot (KM6/KM9 style). fun: surv, event, cumhaz, cloglog.
kmPlot <- function(fit, fun = "surv", ci = FALSE, censor = TRUE, median = FALSE,
                   at = numeric(), risk = FALSE, subtitle = NULL, pal = "jmv",
                   xlab = "Time", xmax = 0, by = 0) {
    df  <- kmData(fit)
    lev <- levels(df$strata)
    cols <- paletteCols(length(lev), pal)
    if (xmax <= 0) xmax <- max(df$time)
    brks <- xBreaks(xmax, by)

    if (fun == "surv") {
        df$y <- df$surv; df$ylo <- df$lower; df$yhi <- df$upper
        ylab <- "Survival probability"
    } else if (fun == "event") {
        df$y <- 1 - df$surv; df$ylo <- 1 - df$upper; df$yhi <- 1 - df$lower
        ylab <- "Cumulative incidence"
    } else if (fun == "cumhaz") {
        df$y <- df$cumhaz; df$ylo <- -log(df$upper); df$yhi <- -log(df$lower)
        ylab <- "Cumulative hazard"
    } else {
        df <- df[df$time > 0 & df$surv > 0 & df$surv < 1, ]
        df$y <- log(-log(df$surv)); df$ylo <- NA; df$yhi <- NA
        ci <- FALSE; median <- FALSE; at <- numeric(); risk <- FALSE
        ylab <- "log(-log survival)"
        xlab <- paste(xlab, "- log scale")
    }

    p <- ggplot(df, aes(time, y, colour = strata))
    if (ci)
        p <- p + geom_ribbon(data = stepRibbon(df),
                             aes(time, ymin = ylo, ymax = yhi, fill = strata),
                             inherit.aes = FALSE, alpha = 0.15, na.rm = TRUE,
                             show.legend = FALSE)
    p <- p + geom_step(linewidth = 0.9)
    if (censor)
        p <- p + geom_point(data = df[df$n.censor > 0, ], shape = 3, size = 1.8,
                            show.legend = FALSE)

    if (median && fun %in% c("surv", "event")) {
        med <- summary(fit)$table
        med <- if (is.matrix(med)) med[, "median"] else med["median"]
        med <- med[!is.na(med)]
        if (length(med))
            p <- p +
                annotate("segment", x = 0, xend = max(med), y = 0.5, yend = 0.5,
                         linetype = 2, colour = "grey45") +
                annotate("segment", x = med, xend = med, y = 0.5, yend = 0,
                         linetype = 2, colour = "grey45")
    }

    if (length(at) && fun %in% c("surv", "event")) {
        s <- summary(fit, times = at, extend = TRUE)
        sd <- data.frame(strata = factor(if (is.null(s$strata)) "All"
                                         else strataNames(as.character(s$strata)), lev),
                         time = s$time,
                         y = if (fun == "surv") s$surv else 1 - s$surv)
        p <- p +
            geom_vline(xintercept = at, linetype = 3, colour = "grey40") +
            geom_point(data = sd, size = 2.4, show.legend = FALSE) +
            geom_label(data = sd, aes(label = scales::percent(y, 1)),
                       hjust = -0.15, size = 3.2, label.size = 0,
                       fill = scales::alpha("white", 0.8), show.legend = FALSE)
    }

    p <- p +
        scale_x_continuous(breaks = if (fun == "cloglog") waiver() else brks,
                           trans = if (fun == "cloglog") "log10" else "identity",
                           expand = expansion(mult = c(0.04, 0.02))) +
        scale_colour_manual(values = cols) +
        scale_fill_manual(values = cols) +
        labs(x = xlab, y = ylab, subtitle = subtitle) +
        plotTheme()
    p <- p + if (fun %in% c("surv", "event"))
                 scale_y_continuous(limits = c(0, 1))
    if (fun != "cloglog")
        p <- p + coord_cartesian(xlim = c(0, xmax))
    p <- p + if (length(lev) == 1) theme(legend.position = "none")
             else legendInside(if (fun == "surv") "topright" else "topleft")

    if (!risk) return(ggplotGrob(p))

    # risk table: second ggplot on the same x scale, stacked with gtable
    s <- summary(fit, times = brks, extend = TRUE)
    rt <- data.frame(strata = factor(if (is.null(s$strata)) "All"
                                     else strataNames(as.character(s$strata)), lev),
                     time = s$time, n = s$n.risk)
    tp <- ggplot(rt, aes(time, strata, label = n, colour = strata)) +
        geom_text(size = 3.5, show.legend = FALSE) +
        scale_y_discrete(limits = rev(lev)) +
        scale_x_continuous(breaks = brks, expand = expansion(mult = c(0.04, 0.02))) +
        coord_cartesian(xlim = c(0, xmax), clip = "off") +
        scale_colour_manual(values = cols) +
        labs(x = NULL, y = NULL, title = "Number at risk") +
        theme_minimal(base_size = 12) +
        theme(panel.grid = element_blank(), axis.text.x = element_blank(),
              plot.title = element_text(size = 11, face = "bold"),
              plot.title.position = "plot")
    g <- rbind(ggplotGrob(p), ggplotGrob(tp), size = "max")
    panels <- g$layout$t[grepl("^panel", g$layout$name)]
    g$heights[panels] <- grid::unit(c(4, 0.28 * length(lev) + 0.3), "null")
    g
}

drawGrob <- function(g) {
    grid::grid.newpage()
    grid::grid.draw(g)
    TRUE
}
