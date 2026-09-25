"""Write fake Delphes-format ROOT files into $BNV2_WORK.

Purely for developing the analysis code before (or instead of) waiting on a
real MG5 run. The branch names and the jagged structure match real Delphes
output; the kinematics are made up and mean nothing.

    python scripts/make_mock_delphes.py

The multiplicities are rigged to the expected pattern -- signal 5 jets / 1 b,
ttbar_semilep 4 jets / 2 b, ttbar_semilep_2j 6 jets / 2 b -- so the sanity
notebook produces the shapes you are hoping to see. That is exactly why you
must not confuse these files with real output: delete work/ before running
the pipeline for real.
"""
import numpy as np, awkward as ak, uproot
from pathlib import Path

rng = np.random.default_rng(7)
import os
WORK = Path(os.environ.get("BNV2_WORK", Path(__file__).resolve().parent.parent / "work"))

def sample(n, njet_mean, nb_true):
    njet = rng.poisson(njet_mean, n).clip(1, 12)
    pt   = ak.unflatten(rng.exponential(60, njet.sum()) + 25, njet)
    eta  = ak.unflatten(rng.normal(0, 1.2, njet.sum()), njet)
    phi  = ak.unflatten(rng.uniform(-np.pi, np.pi, njet.sum()), njet)
    mass = ak.unflatten(rng.exponential(8, njet.sum()), njet)
    flat_b = np.zeros(njet.sum(), dtype=np.int32)
    off = np.concatenate([[0], np.cumsum(njet)])
    for i in range(n):
        k = min(nb_true, njet[i])
        idx = rng.choice(njet[i], size=k, replace=False)
        for j in idx:
            if rng.random() < 0.7:
                flat_b[off[i] + j] = 1
    btag = ak.unflatten(flat_b, njet)
    nmu = rng.integers(0, 2, n)
    nel = np.where(nmu == 0, 1, 0)
    def lep(nl):
        return (ak.unflatten(rng.exponential(40, nl.sum()) + 20, nl),
                ak.unflatten(rng.normal(0, 1.0, nl.sum()), nl),
                ak.unflatten(rng.uniform(-np.pi, np.pi, nl.sum()), nl),
                ak.unflatten(rng.choice([-1, 1], nl.sum()).astype(np.int32), nl))
    mpt, meta, mphi, mq = lep(nmu)
    ept, eeta, ephi, eq = lep(nel)
    met = rng.exponential(45, n) + 5
    return {
        "Jet.PT": pt, "Jet.Eta": eta, "Jet.Phi": phi, "Jet.Mass": mass,
        "Jet.BTag": btag, "Jet.TauTag": ak.zeros_like(btag), "Jet.Flavor": ak.zeros_like(btag),
        "Muon.PT": mpt, "Muon.Eta": meta, "Muon.Phi": mphi, "Muon.Charge": mq,
        "Electron.PT": ept, "Electron.Eta": eeta, "Electron.Phi": ephi, "Electron.Charge": eq,
        "MissingET.MET": ak.unflatten(met, np.ones(n, dtype=int)),
        "MissingET.Eta": ak.unflatten(np.zeros(n), np.ones(n, dtype=int)),
        "MissingET.Phi": ak.unflatten(rng.uniform(-np.pi, np.pi, n), np.ones(n, dtype=int)),
        "ScalarHT.HT": ak.unflatten(ak.to_numpy(ak.sum(pt, axis=1)), np.ones(n, dtype=int)),
    }

for name, njm, nb in [("signal_bnv", 5, 1), ("ttbar_semilep", 4, 2), ("ttbar_semilep_2j", 6, 2)]:
    d = WORK / name / "Events" / "run_test"
    d.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(d / "tag_1_delphes_events.root") as f:
        f["Delphes"] = sample(500, njm, nb)
    print("wrote", d / "tag_1_delphes_events.root")
