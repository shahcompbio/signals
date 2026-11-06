#!/usr/bin/env bash
# Run the signals test suite against the current working tree.
#
# The package needs ~20 dependencies and compiles src/, which is why the suite
# cannot simply be run with Rscript on a dev machine. The signals Docker image
# already carries every dependency plus testthat, so this installs the working
# tree into a throwaway container and runs the real suite there.
#
#   dev/run-tests.sh                # whole suite
#   dev/run-tests.sh phasing        # only test files matching "phasing"
#   SIGNALS_TEST_REPORTER=summary dev/run-tests.sh    # also list warnings
#
# The repo is mounted read-only and the package is built from a copy, so
# compilation artefacts never land in the working tree. test_data/ is excluded
# (it is several GB and not part of the package).
set -euo pipefail
# Resolve the repo from git when possible, so the script still works when it is
# run from outside a checkout that carries dev/ (e.g. kept to one side while
# working on another branch). SIGNALS_REPO overrides both.
REPO="${SIGNALS_REPO:-$(git rev-parse --show-toplevel 2>/dev/null)}"
REPO="${REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
IMAGE="${SIGNALS_TEST_IMAGE:-signals:latest}"
FILTER="${1:-}"

docker run --rm --platform linux/x86_64 \
  -v "$REPO":/src:ro -e FILTER="$FILTER" \
  -e REPORTER="${SIGNALS_TEST_REPORTER:-}" "$IMAGE" bash -c '
set -e
mkdir -p /build
cp -a /src/DESCRIPTION /src/NAMESPACE /build/
for d in R src man tests data; do [ -e "/src/$d" ] && cp -a "/src/$d" /build/; done
# a stale shared object from the host would shadow the sources under test
rm -f /build/src/*.o /build/src/*.so
cd /build
echo "installing working tree..."
if ! R CMD INSTALL --no-docs --no-multiarch . > /tmp/install.log 2>&1; then
  echo "--- install failed ---"; tail -40 /tmp/install.log; exit 1
fi
cd /build/tests
Rscript -e '"'"'
  library(testthat); library(signals)
  f <- Sys.getenv("FILTER"); if (!nzchar(f)) f <- NULL
  r <- Sys.getenv("REPORTER"); if (!nzchar(r)) r <- NULL
  cat("signals", as.character(packageVersion("signals")),
      if (is.null(f)) "- full suite" else paste("- filter:", f), "\n\n")
  if (is.null(r)) test_check("signals", filter = f) else
    test_check("signals", filter = f, reporter = r)
'"'"'
'
