# prognosis: survival analysis for jamovi

A jamovi module with compact output for Kaplan-Meier and Cox regression. It is aimed at teaching (medical students) and reporting. Explanations are optional and hidden by default.

It needs no packages beyond those bundled with jamovi: `survival`, `ggplot2`, `gtable`, `scales`, `jmvcore` and `R6`.

## Analyses (menu: Prognosis)

**Kaplan-Meier**
- Survival summary: N, events, censored, expected events (log-rank), and median with 95% CI.
- Survival at given times (e.g. `12, 36, 60`): number at risk, cumulative events, and S(t) with 95% CI.
- Tests: log-rank, Gehan-Breslow, Tarone-Ware, Peto-Peto, and log-rank trend for ordered groups.
- Plots: survival curve, plus optional cumulative incidence (1 − S), cumulative hazard and log-log plots.
- Plot options: confidence bands, censor marks, number-at-risk table, median lines, survival marks at the given times, test p-value, palette (jamovi / Set1 / Okabe-Ito), and x-axis end and step.

**Cox regression**
- Factors, covariates, interaction terms, and strata.
- Model summary: N, events, and C-index with 95% CI, with an optional bootstrap optimism-corrected C-index.
- Global tests: likelihood ratio, Wald and score.
- Hazard ratios with 95% CI and p. Options: univariable and multivariable side by side, and β, SE and z columns.
- Interactions: likelihood-ratio test per term, and hazard ratios within subgroups.
- Proportional hazards test (Schoenfeld) with an optional residual plot.
- Forest plot (multivariable, or univariable vs multivariable).
- Adjusted survival curves: direct standardisation, with an optional unadjusted KM overlay.

## Coding of the event
- A **nominal** event variable: choose the event level.
- A **continuous** event variable must be coded 0/1 or 1/2 (1 = censored, 2 = event, as in `survival::Surv`).

## Install
In jamovi, go to Modules (⊕) → Sideload and choose `prognosis_0.1.0.jmo`.

## Build (developers)
```sh
# regenerate R/*.h.R after editing jamovi/*.yaml
Rscript -e 'jmvtools::prepare(".")'
# build the .jmo (node and jmvtools required)
JMC=$(Rscript -e 'cat(jmvtools:::jmcPath())' | tr -d '"')
node "$JMC" --build . --home /Applications/jamovi.app --rpath "$(Rscript -e 'cat(R.home("bin"))')"
```
`dev/gallery.R` renders the plot-style gallery that was used to choose the plot designs. It needs the comparison packages installed in `.R/`.

## Notes on methods
- The weighted log-rank tests are implemented in `R/utils.R` using hypergeometric variances. Log-rank and Peto-Peto match `survival::survdiff` (rho = 0, 1). Gehan-Breslow and Tarone-Ware match `survMisc::comp`.
- Subgroup hazard ratios come from contrasts of the interaction model. Other variables are held at their reference level or median.
- The optimism-corrected C-index uses Harrell's bootstrap with a fixed seed (1234).

## Licence
GPL (>= 3). Ideas come from the jamovi modules deathwatch (AGPL-3; no code copied), jYS and jsurvival (GPL ≥ 2). The interaction-builder JavaScript is adapted from jsurvival.
