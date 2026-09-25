#!/usr/bin/env bash
#
# Run the MadGraph5 -> Pythia8 -> Delphes chain for one sample.
#
#   ./scripts/run_mpd.sh <sample> [nevents] [run_tag]
#
# e.g.
#   ./scripts/run_mpd.sh ttbar_semilep 100
#   ./scripts/run_mpd.sh ttbar_semilep_2j 100 test_v2
#
# The process directory is created once in $BNV2_WORK and reused. Delete it
# by hand if you change the proc card (the script cannot tell that the
# generated Fortran is stale).

set -euo pipefail

SAMPLE="${1:?usage: run_mpd.sh <sample> [nevents] [run_tag]}"
NEVENTS="${2:-100}"
RUN_TAG="${3:-run_$(date +%Y%m%d_%H%M%S)}"

: "${BNV2_REPO:?source env.sh first}"
: "${BNV2_WORK:?source env.sh first}"
: "${MG5_ROOT:?source env.sh first}"
: "${DELPHES_PATH:?source env.sh first}"

CARDDIR="$BNV2_REPO/cards/$SAMPLE"
PROCDIR="$BNV2_WORK/$SAMPLE"
MG5="$MG5_ROOT/bin/mg5_aMC"

[[ -d "$CARDDIR" ]] || { echo "no such sample: $CARDDIR" >&2; exit 1; }
[[ -x "$MG5" ]]     || { echo "mg5_aMC not executable at $MG5" >&2; exit 1; }

DELPHES_CARD="$BNV2_REPO/cards/delphes/delphes_card_CMS.tcl"
[[ -f "$DELPHES_CARD" ]] || {
    echo "Delphes card missing. Run ./scripts/bootstrap_delphes_card.sh first." >&2
    exit 1
}

mkdir -p "$BNV2_WORK"

# ---------------------------------------------------------------------------
# 1. Generate the process directory if it does not exist.
#    `output <name>` in the proc card is relative to the CWD, hence the cd.
# ---------------------------------------------------------------------------
if [[ ! -d "$PROCDIR" ]]; then
    echo "==> generating process directory $PROCDIR"
    ( cd "$BNV2_WORK" && "$MG5" "$CARDDIR/proc_card.dat" )
else
    echo "==> reusing existing process directory $PROCDIR"
fi

# ---------------------------------------------------------------------------
# 2. Copy (never symlink) the versioned cards into the process directory.
#    MG5 rewrites these in place during the run.
# ---------------------------------------------------------------------------
echo "==> installing cards"
for card in run_card.dat pythia8_card.dat param_card.dat madspin_card.dat; do
    if [[ -f "$CARDDIR/$card" ]]; then
        cp -v "$CARDDIR/$card" "$PROCDIR/Cards/$card"
    fi
done
cp -v "$DELPHES_CARD" "$PROCDIR/Cards/delphes_card.dat"

# nevents is the one thing that changes run to run, so patch it here rather
# than editing the committed card.
sed -i -E "s|^[[:space:]]*[0-9]+[[:space:]]*=[[:space:]]*nevents|  ${NEVENTS} = nevents|" \
    "$PROCDIR/Cards/run_card.dat"
grep -E "=[[:space:]]*nevents" "$PROCDIR/Cards/run_card.dat"

# ---------------------------------------------------------------------------
# 3. Launch with the shower and detector switches on.
#
#    The switch lines below are what turn Pythia8 and Delphes on
#    non-interactively. If your MG5 build rejects them, run
#        $PROCDIR/bin/generate_events
#    once interactively, toggle "shower" and "detector" in the menu, and
#    confirm afterwards that Events/<tag>/tag_1_delphes_events.root exists.
#    Silently skipping Delphes is the classic failure mode here.
# ---------------------------------------------------------------------------
LAUNCH_SCRIPT="$(mktemp)"
cat > "$LAUNCH_SCRIPT" <<EOF
launch $PROCDIR -n $RUN_TAG
shower=Pythia8
detector=Delphes
analysis=OFF
madspin=OFF
done
EOF

echo "==> launching"
cat "$LAUNCH_SCRIPT"
"$MG5" "$LAUNCH_SCRIPT"
rm -f "$LAUNCH_SCRIPT"

# ---------------------------------------------------------------------------
# 4. Provenance. The MG5 banner records the cards; this records the repo.
# ---------------------------------------------------------------------------
OUTDIR="$PROCDIR/Events/$RUN_TAG"
if [[ -d "$OUTDIR" ]]; then
    {
        echo "sample    : $SAMPLE"
        echo "run_tag   : $RUN_TAG"
        echo "nevents   : $NEVENTS"
        echo "date      : $(date -Is)"
        echo "host      : $(hostname)"
        echo "repo_sha  : $(git -C "$BNV2_REPO" rev-parse HEAD 2>/dev/null || echo unknown)"
        echo "repo_dirty: $(git -C "$BNV2_REPO" status --porcelain 2>/dev/null | wc -l) modified files"
        echo "mg5       : $MG5_ROOT"
        echo "delphes   : $DELPHES_PATH"
        echo "--- card checksums ---"
        md5sum "$CARDDIR"/*.dat "$DELPHES_CARD" 2>/dev/null || true
    } > "$OUTDIR/PROVENANCE.txt"
    echo "==> wrote $OUTDIR/PROVENANCE.txt"

    if compgen -G "$OUTDIR/*delphes_events.root" > /dev/null; then
        echo "==> Delphes output:"
        ls -lh "$OUTDIR"/*delphes_events.root
    else
        echo "!!! No Delphes ROOT file in $OUTDIR -- the detector step did not run." >&2
        ls -l "$OUTDIR" >&2
        exit 2
    fi
else
    echo "!!! Expected output directory $OUTDIR does not exist." >&2
    exit 2
fi
