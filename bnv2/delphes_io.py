"""Read Delphes ROOT output into awkward arrays.

Delphes writes each object collection as a split TClonesArray, so the
branches are flat and named ``Jet.PT``, ``Muon.Charge``, ``MissingET.MET``
and so on. uproot reads these directly with no ROOT dictionary, which is the
whole reason to use it here.

The job of this module is small: pull the branches we care about and zip them
into per-event records with lowercase field names, so downstream code says
``jets.pt`` instead of ``arrays["Jet.PT"]``.

    >>> objs = load_delphes("tag_1_delphes_events.root")
    >>> objs["jet"].pt
    >>> objs["met"].met          # one entry per event, already flattened
"""

from __future__ import annotations

from typing import Any, Iterable, Mapping, Sequence

import awkward as ak
import uproot

__all__ = [
    "DEFAULT_COLLECTIONS",
    "SINGLETON_COLLECTIONS",
    "branch_names",
    "zip_collections",
    "load_delphes",
    "open_tree",
]


#: Collection name -> (Delphes branch prefix, members to read).
#: Trim this per-analysis; reading every member of every collection is the
#: usual reason a "quick look" notebook takes a minute to start.
DEFAULT_COLLECTIONS: dict[str, tuple[str, tuple[str, ...]]] = {
    "jet": ("Jet", ("PT", "Eta", "Phi", "Mass", "BTag", "TauTag", "Flavor")),
    "electron": ("Electron", ("PT", "Eta", "Phi", "Charge")),
    "muon": ("Muon", ("PT", "Eta", "Phi", "Charge")),
    "met": ("MissingET", ("MET", "Eta", "Phi")),
    "scalar_ht": ("ScalarHT", ("HT",)),
}

#: Collections Delphes stores as a one-element array per event. We flatten
#: these so ``met.met`` is a flat array rather than a length-1 jagged one --
#: forgetting this is a reliable source of confusing broadcast errors.
SINGLETON_COLLECTIONS: frozenset[str] = frozenset({"met", "scalar_ht"})


def branch_names(
    collections: Mapping[str, tuple[str, Sequence[str]]] = DEFAULT_COLLECTIONS,
) -> list[str]:
    """Flat list of Delphes branch names implied by a collection spec.

    >>> "Jet.PT" in branch_names()
    True
    """
    names: list[str] = []
    for prefix, members in collections.values():
        names.extend(f"{prefix}.{member}" for member in members)
    return names


def zip_collections(
    arrays: Mapping[str, Any] | ak.Array,
    collections: Mapping[str, tuple[str, Sequence[str]]] = DEFAULT_COLLECTIONS,
    singletons: Iterable[str] = SINGLETON_COLLECTIONS,
) -> dict[str, ak.Array]:
    """Zip flat ``Prefix.Member`` arrays into per-collection records.

    Takes anything indexable by branch name -- an ``ak.Array`` from
    ``TTree.arrays`` or a plain dict -- which keeps this unit testable
    without a ROOT file.

    Field names are lowercased: ``Jet.PT`` -> ``jet.pt``, ``Jet.BTag`` ->
    ``jet.btag``.

    Missing collections are skipped rather than raising, because Delphes
    cards routinely omit collections (no ``Photon`` block, say) and a
    diagnostic notebook should not die over it. A collection that is present
    but missing *some* of its members does raise -- that is a spec bug.
    """
    singletons = set(singletons)
    out: dict[str, ak.Array] = {}

    for name, (prefix, members) in collections.items():
        full = [f"{prefix}.{m}" for m in members]
        present = [f for f in full if _has(arrays, f)]
        if not present:
            continue
        if len(present) != len(full):
            missing = sorted(set(full) - set(present))
            raise KeyError(f"collection {name!r} is missing branches: {missing}")

        record = ak.zip(
            {member.lower(): arrays[f"{prefix}.{member}"] for member in members},
            depth_limit=2,
        )
        out[name] = ak.firsts(record) if name in singletons else record

    return out


def open_tree(path: str, treename: str = "Delphes"):
    """Open a Delphes file and return the TTree.

    ``path`` may include uproot's ``file.root:Delphes`` colon syntax, in
    which case ``treename`` is ignored.
    """
    if ":" in path.rsplit("/", 1)[-1]:
        return uproot.open(path)
    return uproot.open(path)[treename]


def load_delphes(
    path: str,
    treename: str = "Delphes",
    collections: Mapping[str, tuple[str, Sequence[str]]] = DEFAULT_COLLECTIONS,
    entry_stop: int | None = None,
) -> dict[str, ak.Array]:
    """Load one Delphes ROOT file into a dict of awkward record arrays.

    Only the branches named in ``collections`` are read. ``entry_stop``
    limits the number of events, which is what you want the first time you
    point this at an unfamiliar file.
    """
    tree = open_tree(path, treename)
    available = set(tree.keys())
    wanted = [b for b in branch_names(collections) if b in available]
    if not wanted:
        raise KeyError(
            f"none of the requested branches are in {path!r}; "
            f"tree has e.g. {sorted(available)[:10]}"
        )
    arrays = tree.arrays(wanted, entry_stop=entry_stop, library="ak")
    return zip_collections(arrays, collections)


def _has(arrays: Mapping[str, Any] | ak.Array, key: str) -> bool:
    if isinstance(arrays, ak.Array):
        return key in ak.fields(arrays)
    return key in arrays
