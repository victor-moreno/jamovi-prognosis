# prognosis: survival analysis for jamovi

A jamovi module with compact output for Kaplan-Meier and Cox regression. It is aimed at teaching (medical students) and reporting. Explanations are optional and hidden by default.

It needs no packages beyond those bundled with jamovi: `survival`, `ggplot2`, `gtable`, `scales`, `jmvcore` and `R6`.

## Analyses (menu: Prognosis)

**Kaplan-Meier**
- Survival summary: N, events, censored, expected events (log-rank), and median with 95% CI.
- Survival at given times (e.g. `12, 36, 60`): number at risk, cumulative events, and S(t) with 95% CI.
- Tests: log-rank, Gehan-Breslow, Tarone-Ware, Peto-Peto, and log-rank trend for ordered groups.
- Plots: survival curve, plus optional cumulative incidence (1 − S), cumulative hazard and log-log plots.
- Plot options: confidence bands, censor marks, number-at-risk table, median lines, survival marks at the given times, test p-value (under the legend), median values, palette (jamovi / Set1 / Okabe-Ito), and x-axis end and step.

**Cox regression**
- Factors, covariates, interaction terms, and strata.
- Model summary: N, events, and C-index with 95% CI.
- Global tests: likelihood ratio, Wald and score.
- Hazard ratios with 95% CI and p, with reference-level rows. Options: univariable and multivariable side by side (with the univariable Wald p per variable), p for trend, and β, SE and z columns.
- Covariate HRs per unit, per 1 SD, or per a given value (e.g. per 10 years). The covariates are divided by the scale before fitting, so the HR table, univariable fits, subgroup HRs and forest plot use the same unit; model tests and p-values do not change.
- Interactions: likelihood-ratio test per term, and hazard ratios within subgroups.
- Proportional hazards test (Schoenfeld) with an optional residual plot.
- Forest plot (multivariable, or univariable vs multivariable).
- Adjusted survival curves: one plot per factor, by direct standardisation, with an optional unadjusted KM overlay.

## Coding of the event
The event variable must be nominal or ordinal; choose the level that marks the event. A variable coded 0/1 can be used once its measure type is set to nominal. The analysis stays empty until an event level is chosen.

## Install
In jamovi, go to Modules (⊕) → Sideload and choose `prognosis_0.1.0.jmo`.

## Build and install (developers)
```sh
bash tools/install.sh                 # desktop and Docker (container "jamovi")
bash tools/install.sh desktop         # jamovi.app, via jmvtools::install()
bash tools/install.sh docker [name]   # jmc --install in a running container, then a smoke test
bash tools/build-jmo.sh current       # release: dist/*_current_R4.6.0_{x64,arm64}.jmo
bash tools/build-jmo.sh solid         # release: dist/*_solid_R4.5.0_{x64,arm64}.jmo
```
Manual steps:
```sh
# regenerate R/*.h.R after editing jamovi/*.yaml
Rscript -e 'jmvtools::prepare(".")'
# build the .jmo (node and jmvtools required)
JMC=$(Rscript -e 'cat(jmvtools:::jmcPath())' | tr -d '"')
node "$JMC" --build . --home /Applications/jamovi.app --rpath "$(Rscript -e 'cat(R.home("bin"))')"
```
`dev/gallery.R` renders the plot-style gallery that was used to choose the plot designs. It needs the comparison packages installed in `.R/`.

## Results that do not refresh (developer notes)
jamovi re-runs an analysis whenever any option changes. It restores each result from the previous run unless an option listed in that item's `clearWith` (in `jamovi/*.r.yaml`) has changed. Both analyses follow this pattern:

1. **`clearWith`** lists only the options that change a table's content. For example, the Cox HR table clears with the model variables, `uniMulti` and `trend`, but not with plot options or `coefDetails`. The β/SE/z columns are always filled and only shown or hidden.
2. **`.init`** lays out every table: it computes rows from the options and data (factor levels, model-matrix columns, interaction subgroups) and adds them with their labels. It also sets image sizes and hides all results until the required variables are selected.
3. **`.run`** fills a table only if `table$isNotFilled()`, i.e. it was cleared. Plot states are rebuilt only when `image$state` is `NULL`.
4. **No `NA` in cells.** jmvcore treats an `NA` cell as not filled, so the table would be recomputed on every run. Missing values are written as `""` (`blankNA()`), which counts as filled and renders empty.
5. **Notes** are set on every run, so explanations toggle without refilling the tables.
6. **Keep plot states small.** Store the data the plot needs (for example the adjusted curves), not fitted models.
7. **Reserved names:** never name an option `palette` or a private method `.formula`, because jmvcore uses both.

`tests/test_freeze.R` (run after `R CMD INSTALL -l .tmp/lib .`) checks this by serialising the results, changing one option, restoring them the way the jamovi engine does, and reporting which tables need refilling.

## Notes on methods
- The weighted log-rank tests are implemented in `R/utils.R` using hypergeometric variances. Log-rank and Peto-Peto match `survival::survdiff` (rho = 0, 1). Gehan-Breslow and Tarone-Ware match `survMisc::comp`.
- Subgroup hazard ratios come from contrasts of the interaction model. Other variables are held at their reference level or median.

## Licence
GPL (>= 3). Ideas come from the jamovi modules deathwatch (AGPL-3; no code copied), jYS and jsurvival (GPL ≥ 2). The interaction-builder JavaScript is adapted from jsurvival.
