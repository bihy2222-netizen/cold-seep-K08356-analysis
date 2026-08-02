#!/usr/bin/env bash
set -euo pipefail

# Compatibility entry point. Separate habitat-specific MAG references are not
# comparable across habitats, so all new runs use the common MAG48 workflow.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

printf '%s\n' \
  'NOTICE: habitat-specific derep references are deprecated for comparison.' \
  'Forwarding to the unified 48-MAG reference workflow.' >&2

exec "${SCRIPT_DIR}/02_run_unified_MAG48_coverm_tpm.sh" "$@"
