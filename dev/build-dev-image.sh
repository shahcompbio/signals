#!/usr/bin/env bash
# Build a dev image with the current working tree's signals installed, on top of
# the dependency layers of signals:latest. No GitHub token needed.
set -euo pipefail
REPO="${SIGNALS_REPO:-$(git rev-parse --show-toplevel 2>/dev/null)}"
REPO="${REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
TAG="${1:-signals:dev}"
cd "$REPO"
# .dockerignore keeps the multi-GB test datasets out of the build context
docker build --platform linux/x86_64 -f dev/Dockerfile.dev -t "$TAG" .
echo "built $TAG"
