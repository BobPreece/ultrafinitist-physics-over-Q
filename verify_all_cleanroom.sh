#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS=(
  "Run_1_standalone_ultrafinitist_qft.py"
  "Run_2_standalone_rational_mechanics.py"
  "Run_3_standalone_ultrafinitist_turbulence.py"
  "Run_4_standalone_ultrafinitist_gravity.py"
  "Run_5_standalone_ultrafinitist_statmech.py"
)

echo "=== Zero-Trust /tmp Clean-Room Verification (5 Pure-Q Exhibits) ==="
for script in "${SCRIPTS[@]}"; do
  SANDBOX="$(mktemp -d /tmp/cleanroom_ultrafinitist_XXXXXX)"
  cp "$ROOT_DIR/$script" "$SANDBOX/$script"
  echo "-> Running /usr/bin/python3 -I $script in $SANDBOX ..."
  (cd "$SANDBOX" && /usr/bin/python3 -I "$script")
  rm -rf "$SANDBOX"
  echo "   [PASS] $script"
done
echo "=== ALL 5 EXHIBITS PASSED CLEAN-ROOM ISOLATION ==="
