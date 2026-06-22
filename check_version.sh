#!/usr/bin/env bash
# Fails if README.md's stated version drifts from APP_VERSION in
# chimera_functions.R (the single source of truth, also read by app.R
# and chimera_cli.R). Run directly, or via the .git/hooks/pre-commit hook.

set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

code_version=$(grep -m1 -oE 'APP_VERSION <- "[0-9]+\.[0-9]+\.[0-9]+"' "$DIR/chimera_functions.R" | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')
readme_version=$(grep -m1 -oE '_version [0-9]+\.[0-9]+\.[0-9]+_' "$DIR/README.md" | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')

if [[ -z "$code_version" || -z "$readme_version" ]]; then
  echo "check_version.sh: could not find a version string in chimera_functions.R or README.md" >&2
  exit 1
fi

if [[ "$code_version" != "$readme_version" ]]; then
  echo "Version mismatch: APP_VERSION is $code_version but README.md says $readme_version" >&2
  exit 1
fi

echo "Version OK: $code_version"
