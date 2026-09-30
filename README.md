# prognosis: survival analysis for jamovi

A jamovi module with compact output for Kaplan-Meier and Cox regression. It is aimed at teaching (medical students) and reporting. Explanations are optional and hidden by default.

It needs no packages beyond those bundled with jamovi: `survival`, `ggplot2`, `gtable`, `scales`, `jmvcore` and `R6`.

## Analyses (menu: Prognosis)

**Kaplan-Meier**
- Survival summary: N, events, censored, and median survival with 95% CI.
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

## Tutorial
A [website for medical students](https://victor-moreno.github.io/jamovi-prognosis/) (in `docs/`): time-to-event data and censoring, the Kaplan-Meier estimator and the log-rank test (each worked by hand on ten patients), Cox regression, proportional hazards, interactions, adjusted curves and reporting. The worked examples run this module on a teaching dataset (`docs/data/lung.omv`, from `survival::lung`) that students open in jamovi, so every table and plot is what jamovi shows. The site is rendered with Quarto from a source kept outside this repository.

## Languages
The options, tables, notes and plots are available in English, Spanish and Catalan. jamovi uses the language chosen in its settings. The translations are in `jamovi/i18n/*.po`.

## Coding of the event
The event variable must be nominal or ordinal; choose the level that marks the event. A variable coded 0/1 can be used once its measure type is set to nominal. The analysis stays empty until an event level is chosen.

## Install
Prebuilt `.jmo` files are attached to the [Releases](../../releases) page — four per version, one for each combination of jamovi series and CPU:

| your jamovi | Apple silicon | Intel / AMD |
| --- | --- | --- |
| **current** (bundles R 4.6.0) | `prognosis_<version>_current_R4.6.0_arm64.jmo` | `prognosis_<version>_current_R4.6.0_x64.jmo` |
| **solid** (bundles R 4.5.0) | `prognosis_<version>_solid_R4.5.0_arm64.jmo` | `prognosis_<version>_solid_R4.5.0_x64.jmo` |

The same file works on macOS, Windows and Linux: jamovi checks the R version and the CPU, not the operating system (see Help → About). In jamovi, go to Modules (⊕) → Sideload and choose the `.jmo`.

## Build and install (developers)
```sh
bash tools/install.sh                 # desktop and Docker (container "jamovi")
bash tools/install.sh desktop         # jamovi.app, via jmvtools::install()
bash tools/install.sh docker [name]   # jamovi >= 28.4 container, via ../jamovi-src's compiler
bash tools/build-jmo.sh current       # release: dist/*_current_R4.6.0_{x64,arm64}.jmo
bash tools/build-jmo.sh solid         # release: dist/*_solid_R4.5.0_{x64,arm64}.jmo
bash tools/release.sh                 # build all four and publish one GitHub release
```
After each install, `tools/smoke.R` checks the installed module against `survival` (log-rank χ², Cox HRs with an ordinal factor). The other checks in `tests/` run against a scratch install: `R CMD INSTALL -l .tmp/lib .`, then `Rscript tests/<file>.R`.
Manual steps:
```sh
# regenerate R/*.h.R after editing jamovi/*.yaml
Rscript -e 'jmvtools::prepare(".")'
# build the .jmo (node and jmvtools required)
JMC=$(Rscript -e 'cat(jmvtools:::jmcPath())' | tr -d '"')
node "$JMC" --build . --home /Applications/jamovi.app --rpath "$(Rscript -e 'cat(R.home("bin"))')"
```

## Notes on methods
- The weighted log-rank tests are implemented in `R/utils.R` using hypergeometric variances. Log-rank and Peto-Peto match `survival::survdiff` (rho = 0, 1). Gehan-Breslow and Tarone-Ware match `survMisc::comp`.
- Subgroup hazard ratios come from contrasts of the interaction model. Other variables are held at their reference level or median.

## Acknowledgements
This module was designed and built together with Claude Code, over several rounds of design, implementation and testing in jamovi.

Ideas come from the jamovi modules deathwatch (AGPL-3; no code copied), jYS and jsurvival (GPL ≥ 2). The interaction-builder JavaScript is adapted from jsurvival. Kaplan-Meier estimates and Cox models use the R package `survival` by Terry Therneau.

This module has been developed with support of the Instituto de Salud Carlos III (ISCIII), “Programa FORTALECE del Ministerio de Ciencia e Innovación”, through the project number FORT23/00032 and the Consortium for Biomedical Research in Epidemiology and Public Health (CIBERESP), action Genrisk.

## Licence
GPL (>= 3).
