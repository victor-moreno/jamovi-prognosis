#!/usr/bin/env bash
# Build and install prognosis into jamovi desktop and/or a running jamovi Docker
# container.
#
#   bash tools/install.sh              both targets, whichever are available
#   bash tools/install.sh desktop
#   bash tools/install.sh docker [container]     (default container: jamovi)
#
# The desktop target uses whichever R `Rscript` resolves to (respecting
# ~/.Rprofile, which appends jamovi.app's bundled module library). R here is
# managed by rig: pick the active R version with `rig default <version>` before
# running this if needed; it must match the R version jamovi.app itself
# bundles, or jmvcore segfaults on load.
set -euo pipefail

TARGET="${1:-both}"
CONTAINER="${2:-jamovi}"

# prognosis is the module itself (no subfolder as in jmvplus)
HERE="$(cd "$(dirname "$0")/.." && pwd)"
MODULE="$(awk -F': *' '$1 == "Package" { print $2; exit }' "$HERE/DESCRIPTION")"
VERSION="$(awk -F': *' '$1 == "Version" { print $2; exit }' "$HERE/DESCRIPTION")"
ARTIFACT="$HERE/${MODULE}_${VERSION}.jmo"

# ── desktop ──────────────────────────────────────────────────────────────────
install_desktop() {
  local APP APP_R LOG
  APP=/Applications/jamovi.app
  APP_R="$APP/Contents/Frameworks/R.framework/Versions/Current/Resources/bin/R"
  [ -x "$APP_R" ] || { echo "error: no R inside $APP" >&2; return 1; }

  echo ">> desktop: building $MODULE with $(Rscript -e 'cat(R.version.string)')"
  cd "$HERE"

  # log kept inside the project, not the system /tmp
  mkdir -p "$HERE/.tmp"
  LOG="$(mktemp "$HERE/.tmp/install.XXXXXX")"
  Rscript -e 'jmvtools::install()' 2>&1 | tee "$LOG" | grep -vE '^\s*$' || true

  # jmvtools::install() can report errors on stdout while exiting successfully.
  # It can also claim installation succeeded after a SingletonLock failure.
  [ -f "$ARTIFACT" ] || {
    echo "error: jmvtools did not produce $ARTIFACT" >&2
    rm -f "$LOG"; return 1
  }
  if grep -q 'SingletonLock' "$LOG"; then
    echo
    echo "!! jamovi.app could not be driven (SingletonLock denied)."
    echo "!! The .jmo was still built. Install it by hand:"
    echo "!!   jamovi -> Modules -> Install from file -> $ARTIFACT"
    rm -f "$LOG"
    return 0
  fi
  if ! grep -q 'Module installed successfully' "$LOG"; then
    echo "error: jmvtools::install() did not install the module (see above)" >&2
    rm -f "$LOG"; return 1
  fi
  rm -f "$LOG"

  local MODDIR="$HOME/Library/Application Support/jamovi/modules/$MODULE"
  # jamovi unpacks the module shortly after jmvtools reports success; wait for it
  local i
  for i in $(seq 1 20); do
    [ -d "$MODDIR" ] && break
    sleep 1
  done
  if [ -d "$MODDIR" ]; then
    echo ">> desktop: installed at $MODDIR"
  else
    echo "!! desktop: install reported success but $MODDIR does not exist."
    echo "!! Install $ARTIFACT by hand (Modules -> Install from file)."
    return 1
  fi
}

# ── docker ───────────────────────────────────────────────────────────────────
install_docker() {
  if ! docker ps --format '{{.Names}}' 2>/dev/null | grep -qx "$CONTAINER"; then
    echo "!! docker: container '$CONTAINER' is not running — skipping"
    return 0
  fi

  if ! docker exec "$CONTAINER" sh -c 'command -v jmc >/dev/null 2>&1'; then
    echo "!! docker: jmc is not in the container." >&2
    echo "!! Install the jamovi compiler in the image before using this target." >&2
    return 1
  fi

  echo ">> docker: copying source into $CONTAINER"
  # --no-mac-metadata/--no-xattrs: AppleDouble ._ files otherwise land in the
  # container and jmc tries to compile them.
  tar --no-mac-metadata --no-xattrs -C "$HERE" -cf - DESCRIPTION NAMESPACE R jamovi \
    | docker exec -i "$CONTAINER" sh -c \
        "rm -rf /tmp/$MODULE-src && mkdir -p /tmp/$MODULE-src && tar -C /tmp/$MODULE-src -xf -"

  echo ">> docker: jmc --install"
  docker exec -i -e MODULE="$MODULE" "$CONTAINER" bash -s <<'INCONTAINER'
set -euo pipefail
source /usr/lib/jamovi/bin/env.conf 2>/dev/null || true
RHOME="${R_HOME:-$(R RHOME 2>/dev/null || true)}"
[ -n "$RHOME" ] || { echo "   error: no R in the container" >&2; exit 1; }
RLIBS=/usr/lib/jamovi/modules/base/R

jmc --install "/tmp/$MODULE-src" \
    --to /usr/lib/jamovi/modules \
    --rhome "$RHOME" \
    --rlibs "$RLIBS" \
    --patch-version --skip-deps

[ -f "/usr/lib/jamovi/modules/$MODULE/jamovi.yaml" ] || {
  echo "   error: jmc did not install $MODULE" >&2; exit 1; }
INCONTAINER

  echo ">> docker: restarting $CONTAINER to load the module"
  docker restart "$CONTAINER" >/dev/null
  echo ">> docker: smoke test (km and cox against survival)"
  docker exec -i "$CONTAINER" bash -s <<'INCONTAINER'
set -euo pipefail
Rscript --vanilla -e '
    .libPaths(c(
        "/usr/lib/jamovi/modules/base/R",
        "/usr/lib/jamovi/modules/prognosis/R",
        .libPaths()
    ))
    library(prognosis)

    d <- survival::lung
    d <- d[!is.na(d$ph.ecog) & d$ph.ecog < 3, ]
    d$event <- factor(d$status, 1:2, c("Alive", "Dead"))
    d$ecog <- factor(d$ph.ecog)
    d$sex <- factor(d$sex, 1:2, c("Male", "Female"))

    km <- prognosis::km(data = d, elapsed = "time", event = "event",
                        eventLevel = "Dead", group = "ecog")
    lr <- km$tests$asDF$chisq[1]
    ref <- survival::survdiff(survival::Surv(time, status) ~ ecog, data = d)$chisq
    stopifnot(isTRUE(all.equal(lr, ref)))

    cx <- prognosis::cox(data = d, elapsed = "time", event = "event",
                         eventLevel = "Dead", factors = c("sex", "ecog"), covs = "age")
    hr <- cx$coefTable$asDF$hr
    fit <- survival::coxph(survival::Surv(time, status) ~ sex + ecog + age, data = d)
    stopifnot(isTRUE(all.equal(unname(hr), unname(exp(coef(fit))))))
    cat(sprintf("   smoke test passed: log-rank chi2 %.2f, %d hazard ratios\n", lr, length(hr)))
'
INCONTAINER
  echo ">> docker: installed $MODULE; open Prognosis -> Kaplan-Meier to verify"
}

case "$TARGET" in
  desktop) install_desktop ;;
  docker)  install_docker ;;
  both)    install_desktop || true; echo; install_docker || true ;;
  *)       echo "usage: install.sh [desktop|docker|both] [container]" >&2; exit 1 ;;
esac
