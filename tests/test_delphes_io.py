"""Unit tests for bnv2.delphes_io.

Two levels of mock:

1. A plain dict keyed by ``Prefix.Member``, which is what ``zip_collections``
   actually consumes. Cheap, and covers the naming/flattening logic.
2. A real (tiny) ROOT file written with uproot, with Delphes-style dotted
   branch names. Covers the ``load_delphes`` path end to end without needing
   a 100 MB Delphes output to hand.
"""

import awkward as ak
import numpy as np
import pytest
import uproot

from bnv2 import delphes_io as dio


# --------------------------------------------------------------------------
# mock data
# --------------------------------------------------------------------------

N_EVENTS = 4

MOCK_BRANCHES = {
    # 4 events with 3, 0, 5, 4 jets
    "Jet.PT": ak.Array([[90.0, 50.0, 25.0], [], [120.0, 80.0, 60.0, 40.0, 22.0],
                        [75.0, 55.0, 35.0, 21.0]]),
    "Jet.Eta": ak.Array([[0.1, -1.2, 2.0], [], [0.0, 0.5, -0.5, 1.5, -2.1],
                         [0.3, -0.3, 1.1, -1.1]]),
    "Jet.Phi": ak.Array([[0.0, 1.5, -2.0], [], [0.2, 1.0, -1.0, 2.5, -2.5],
                         [0.0, 3.0, -0.7, 1.9]]),
    "Jet.Mass": ak.Array([[10.0, 7.0, 4.0], [], [15.0, 9.0, 8.0, 5.0, 3.0],
                          [11.0, 8.0, 6.0, 4.0]]),
    "Jet.BTag": ak.Array([[1, 0, 0], [], [1, 1, 0, 0, 0], [0, 1, 0, 0]]),
    "Jet.TauTag": ak.Array([[0, 0, 0], [], [0, 0, 0, 0, 0], [0, 0, 0, 0]]),
    "Jet.Flavor": ak.Array([[5, 21, 1], [], [5, 5, 21, 2, 3], [21, 5, 1, 2]]),
    # one muon in events 0 and 2, none elsewhere
    "Muon.PT": ak.Array([[45.0], [], [33.0], []]),
    "Muon.Eta": ak.Array([[0.4], [], [-1.1], []]),
    "Muon.Phi": ak.Array([[2.2], [], [0.6], []]),
    "Muon.Charge": ak.Array([[1], [], [-1], []]),
    # one electron in event 1
    "Electron.PT": ak.Array([[], [61.0], [], []]),
    "Electron.Eta": ak.Array([[], [0.9], [], []]),
    "Electron.Phi": ak.Array([[], [-1.4], [], []]),
    "Electron.Charge": ak.Array([[], [-1], [], []]),
    # Delphes stores MissingET and ScalarHT as length-1 arrays per event
    "MissingET.MET": ak.Array([[30.0], [12.0], [95.0], [44.0]]),
    "MissingET.Eta": ak.Array([[0.0], [0.0], [0.0], [0.0]]),
    "MissingET.Phi": ak.Array([[-1.0], [2.7], [0.3], [1.6]]),
    "ScalarHT.HT": ak.Array([[165.0], [61.0], [322.0], [186.0]]),
}


@pytest.fixture
def mock_file(tmp_path):
    path = tmp_path / "mock_delphes_events.root"
    with uproot.recreate(path) as f:
        f["Delphes"] = {k: v for k, v in MOCK_BRANCHES.items()}
    return str(path)


# --------------------------------------------------------------------------
# branch_names
# --------------------------------------------------------------------------

def test_branch_names_are_prefix_dot_member():
    names = dio.branch_names()
    assert "Jet.PT" in names
    assert "MissingET.MET" in names
    assert all("." in n for n in names)


def test_branch_names_respects_a_custom_spec():
    spec = {"jet": ("Jet", ("PT",))}
    assert dio.branch_names(spec) == ["Jet.PT"]


# --------------------------------------------------------------------------
# zip_collections
# --------------------------------------------------------------------------

def test_zip_collections_field_names_are_lowercased():
    objs = dio.zip_collections(MOCK_BRANCHES)
    assert set(ak.fields(objs["jet"])) == {
        "pt", "eta", "phi", "mass", "btag", "tautag", "flavor"
    }


def test_zip_collections_preserves_jaggedness_and_values():
    objs = dio.zip_collections(MOCK_BRANCHES)
    assert ak.to_list(ak.num(objs["jet"], axis=1)) == [3, 0, 5, 4]
    assert objs["jet"].pt[2][0] == pytest.approx(120.0)
    assert objs["muon"].charge[0][0] == 1


def test_zip_collections_flattens_singletons():
    objs = dio.zip_collections(MOCK_BRANCHES)
    # met is one value per event, not a length-1 list per event
    assert objs["met"].met.ndim == 1
    assert ak.to_list(objs["met"].met) == pytest.approx([30.0, 12.0, 95.0, 44.0])
    assert ak.to_list(objs["scalar_ht"].ht) == pytest.approx([165.0, 61.0, 322.0, 186.0])


def test_zip_collections_skips_absent_collections():
    partial = {k: v for k, v in MOCK_BRANCHES.items() if not k.startswith("Muon.")}
    objs = dio.zip_collections(partial)
    assert "muon" not in objs
    assert "jet" in objs


def test_zip_collections_raises_on_partially_present_collection():
    broken = {k: v for k, v in MOCK_BRANCHES.items() if k != "Jet.Eta"}
    with pytest.raises(KeyError, match="Jet.Eta"):
        dio.zip_collections(broken)


def test_zip_collections_accepts_an_awkward_array():
    arr = ak.Array({k: v for k, v in MOCK_BRANCHES.items()})
    objs = dio.zip_collections(arr)
    assert ak.to_list(ak.num(objs["jet"], axis=1)) == [3, 0, 5, 4]


# --------------------------------------------------------------------------
# load_delphes (round trip through a real ROOT file)
# --------------------------------------------------------------------------

def test_load_delphes_round_trip(mock_file):
    objs = dio.load_delphes(mock_file)
    assert set(objs) == {"jet", "electron", "muon", "met", "scalar_ht"}
    assert ak.to_list(ak.num(objs["jet"], axis=1)) == [3, 0, 5, 4]
    assert ak.to_list(objs["met"].met) == pytest.approx([30.0, 12.0, 95.0, 44.0])
    assert ak.to_list(ak.sum(objs["jet"].btag, axis=1)) == [1, 0, 2, 1]


def test_load_delphes_entry_stop(mock_file):
    objs = dio.load_delphes(mock_file, entry_stop=2)
    assert len(objs["jet"]) == 2


def test_load_delphes_custom_collections(mock_file):
    spec = {"jet": ("Jet", ("PT", "BTag"))}
    objs = dio.load_delphes(mock_file, collections=spec)
    assert set(objs) == {"jet"}
    assert set(ak.fields(objs["jet"])) == {"pt", "btag"}


def test_load_delphes_raises_when_nothing_matches(tmp_path):
    path = tmp_path / "wrong.root"
    with uproot.recreate(path) as f:
        f["Delphes"] = {"Something.Else": np.arange(3, dtype=float)}
    with pytest.raises(KeyError, match="none of the requested branches"):
        dio.load_delphes(str(path))


def test_open_tree_accepts_colon_syntax(mock_file):
    tree = dio.open_tree(mock_file + ":Delphes")
    assert "Jet.PT" in tree.keys()
