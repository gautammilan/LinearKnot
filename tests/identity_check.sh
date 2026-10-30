#!/bin/bash
# V and G outside identity (--out-id) at the paper preset with b_V and b_G uncapped:
# sum_{p',q'} C_pre(p'-1) G(p',q') Sbar_pre(q'+1) = Qb_b(x,y) for every kept pair, same for V.
# Fails if either max relative error exceeds TOL or a line is missing.
#   tests/identity_check.sh BIN PARAM_FILE [TOL]
set -euo pipefail
BIN=$1
PF=$2
TOL=${3:-1e-9}
FA=$(dirname "$0")/data/RNaseP_D.radiodurans_1-150.fa
[ -f "$PF" ] || { echo "identity check: parameter file $PF not found (set PARAM_FILE=...)"; exit 2; }
ERR=$(mktemp)
trap 'rm -f "$ERR"' EXIT
"$BIN" --preset accurate --material rna1999 --dangles some --param-file "$PF" \
  --bV 10000000 --bG 10000000 --out-id --no-bpp < "$FA" > /dev/null 2> "$ERR"
status=0
for t in V G-states; do
  line=$(grep -E "^\[out-id\] $t pairs" "$ERR" | head -1 || true)
  if [ -z "$line" ]; then
    echo "identity check: no [out-id] $t line"
    status=1
    continue
  fi
  rel=$(echo "$line" | sed -E 's/.*max_rel ([^ ]+).*/\1/')
  if awk -v r="$rel" -v t="$TOL" 'BEGIN { exit !(r + 0 <= t + 0) }'; then
    echo "identity check: $line  ok (<= $TOL)"
  else
    echo "identity check: $line  FAIL (> $TOL)"
    status=1
  fi
done
exit $status
