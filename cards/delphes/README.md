# Delphes card

This directory holds the vendored detector configuration. It starts empty,
which is why it may not have survived the unzip — archives and git both drop
directories with no files in them.

Populate it with:

```bash
source env.sh
./scripts/bootstrap_delphes_card.sh
```

That copies `$DELPHES_PATH/cards/delphes_card_CMS.tcl` to
`delphes_card_CMS.tcl` here and writes a `PROVENANCE.txt` recording the
source path and Delphes version. Commit both. From that point the repo — not
your Delphes install — is the source of truth for the detector
configuration, and `run_mpd.sh` copies this file into
`PROCDIR/Cards/delphes_card.dat` on every run.

The card is copied from your install rather than shipped with the repo
because the tcl syntax and the available modules track the Delphes version.

Before trusting any plots, check three things in it:

- `JetFinder` `ParameterR` — 0.4 for CMS anti-kT, and it must match `drjj` in
  `cards/ttbar_semilep_2j/run_card.dat`
- the `BTagging` efficiency formula, and which working point it represents
- the jet pT threshold in both `JetFinder` and `TreeWriter`
