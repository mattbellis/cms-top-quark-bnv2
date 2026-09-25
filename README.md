# cms-top-quark-bnv2

MC studies for baryon-number-violating (|ΔB| = 2) top quark decays, stage 1:
MadGraph5_aMC@NLO → Pythia8 → Delphes (MPD).

## What this repo is and is not

**In the repo (version controlled):**

- `cards/` — every input card that defines a sample: proc cards, run cards,
  Pythia8 cards, param cards, and the Delphes card.
- `scripts/` — thin drivers that glue the cards to an MG5 installation.
- `bnv2/` — the analysis package (uproot + awkward + matplotlib).
- `tests/` — pytest unit tests on mock data.
- `notebooks/` — interactive analysis.

**Not in the repo (build artifacts, gitignored):**

- MG5 process directories (`ttbar_semilep/`, `SubProcesses/`, `Events/`, …)
- LHE / HepMC / Delphes ROOT files
- The MG5 and Delphes installations themselves

An MG5 process directory is generated output — several hundred MB of
auto-written Fortran plus the events. It is reproducible from `cards/` plus a
git SHA, so committing it buys nothing and costs a lot.

## The layout on disk

```
$BNV2_REPO/            # this repo, wherever you cloned it
$BNV2_WORK/            # scratch; MG5 process dirs land here
    ttbar_semilep/
        Cards/         # copies of $BNV2_REPO/cards/ttbar_semilep/*
        Events/run_01/ # unweighted_events.lhe.gz, tag_1_delphes_events.root
    ttbar_semilep_2j/
$MG5_ROOT/             # MadGraph5_aMC@NLO 3.6.2 install
$DELPHES_PATH/         # Delphes install
```

Copy `env.sh.example` to `env.sh` (gitignored), edit the four paths, and
`source env.sh`.

## Answering the three questions directly

**Do we store cards in the repo and link them into the MG5 directory?**

Store them in the repo, yes. **Copy** them into `PROCDIR/Cards/`, do not
symlink. MG5 rewrites cards in place: `generate_events` normalizes the run
card, writes back defaults for anything you omitted, stamps the seed, and
regenerates `param_card.dat` from the model. If `Cards/run_card.dat` is a
symlink into your working tree, every run silently mutates your canonical
card and your git status goes dirty in a way that is very easy to commit by
accident. Copy in, and let the run banner in `Events/run_01/` be the
authoritative record of what was actually used.

**Do we call `bin/mg5_aMC` from the repo?**

You call `$MG5_ROOT/bin/mg5_aMC` and hand it a proc card that lives in the
repo. The one thing to know: `output <name>` in a proc card is resolved
relative to the *current working directory* of `mg5_aMC`, not to the proc
card's location. So `scripts/run_mpd.sh` cds into `$BNV2_WORK` first, which
is why the proc cards can just say `output ttbar_semilep` with no absolute
paths and stay portable across machines.

**So the full loop is:**

```bash
source env.sh
./scripts/bootstrap_delphes_card.sh          # once, to vendor the CMS Delphes card
./scripts/run_mpd.sh ttbar_semilep    100    # generate + shower + detector
./scripts/run_mpd.sh ttbar_semilep_2j 100
```

`run_mpd.sh` does: generate the process dir (only if missing) → copy cards in
→ override `nevents` → launch with Pythia8 and Delphes on → write a
`PROVENANCE.txt` recording the repo SHA and the card checksums next to the
events.

## Samples

| Sample | Process | Final state | Purpose |
|---|---|---|---|
| `signal_bnv` | pp → tt̄, t → ūd̄d̄d̄d̄, t̄ → b̄ℓ⁻ν̄ | 5j + ℓ + ν, **1 b** | signal (needs the BNV UFO; see `cards/signal_bnv/README.md`) |
| `ttbar_semilep` | pp → tt̄ semileptonic | 4j + ℓ + ν, 2 b | irreducible-ish, below signal in jet multiplicity |
| `ttbar_semilep_2j` | pp → tt̄ + 2j semileptonic | 6j + ℓ + ν, 2 b | brackets the signal from above |

The two backgrounds bracket the signal's jet multiplicity, which is the
point. The handle that actually separates them is **b-tag multiplicity**: the
BNV vertex t → ūd̄d̄d̄d̄ has no b, so the signal has exactly one true b-jet
against two for both backgrounds. MET is *not* a discriminant here — all
three have one real neutrino — so do not expect the sanity plots to show
separation there, and be suspicious if they do.

## Caveats you should decide on before scaling up

- **`ttbar_semilep_2j` is exclusive, not merged.** It is a fixed-order 2-jet
  matrix element with a `ptj` cut, no MLM matching. That is the right thing
  for a shape comparison at fixed jet multiplicity and the wrong thing for
  anything resembling a rate. If you later want an inclusive tt̄+jets sample,
  you need 0/1/2-jet merged with `ickkw=1`, `xqcut`, and matching settings in
  the Pythia8 card.
- **Taus are excluded** from `l+`/`l-` (MG5's default multiparticle). Real
  τ → ℓνν feed-down into the signal region is not modeled.
- **`cut_decays = False`** in both run cards, so the `ptj`/`drjj` cuts apply
  only to the matrix-element jets and not to W decay products. Without this
  the 2-jet sample would be cutting into the hadronic W and distorting exactly
  the distributions you want to compare.
- **Delphes b-tagging is a parameterization**, not a tagger. Fine for stage 1;
  do not tune a selection to the third decimal place on it.

## Analysis

```bash
pip install -e .
pytest
jupyter lab notebooks/01_delphes_sanity.ipynb
```

## Developing before you have events

You do not need to wait for MG5. `scripts/make_mock_delphes.py` writes fake
Delphes-format ROOT files into `$BNV2_WORK` with the right branch names and
jagged structure, so the notebook and the analysis code can be developed and
debugged immediately. Delete `work/` before you run the real pipeline.
