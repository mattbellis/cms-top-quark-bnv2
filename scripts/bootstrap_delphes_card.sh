#!/usr/bin/env bash
#
# Vendor the CMS Delphes card into the repo so it is version controlled
# alongside everything else.
#
# We copy from your Delphes install rather than shipping a card, because the
# tcl syntax and the available modules track the Delphes version. Copy it
# once, commit it, and from then on the repo -- not your install -- is the
# source of truth for the detector configuration.

set -euo pipefail

: "${BNV2_REPO:?source env.sh first}"
: "${DELPHES_PATH:?source env.sh first}"

SRC="$DELPHES_PATH/cards/delphes_card_CMS.tcl"
DEST="$BNV2_REPO/cards/delphes/delphes_card_CMS.tcl"

[[ -f "$SRC" ]] || { echo "not found: $SRC" >&2; ls "$DELPHES_PATH/cards" >&2; exit 1; }

mkdir -p "$(dirname "$DEST")"
if [[ -f "$DEST" ]]; then
    echo "already vendored: $DEST"
    diff -q "$SRC" "$DEST" || echo "  (differs from $SRC -- keeping the repo copy)"
else
    cp -v "$SRC" "$DEST"
    {
        echo "# Vendored from: $SRC"
        echo "# Delphes version: $("$DELPHES_PATH/DelphesHepMC3" --version 2>&1 | head -1 || echo unknown)"
        echo "# Date: $(date -Is)"
    } > "$(dirname "$DEST")/PROVENANCE.txt"
fi

echo
echo "Things to check in the card before you trust the plots:"
echo "  - JetFinder ParameterR   (0.4 for CMS anti-kT; must match drjj in the run card)"
echo "  - BTagging efficiency formula and the working point it represents"
echo "  - Jet pT threshold in the JetFinder and in TreeWriter"
