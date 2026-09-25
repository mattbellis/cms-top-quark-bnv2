"""Kinematic helpers on (pt, eta, phi, mass) objects.

Everything here is written against numpy ufuncs, so it works unchanged on
plain numpy arrays and on jagged awkward arrays.

A note on tooling: ``vector`` (scikit-hep) does all of this properly, with
Lorentz boosts, coordinate conversions and awkward behaviours registered on
the records. If this module starts growing boosts or rest-frame variables,
switch to ``vector`` rather than extending it -- the point of keeping it
hand-rolled is that stage-1 sanity plots then have no dependency beyond
uproot/awkward/numpy and the tests stay trivially reproducible.
"""

from __future__ import annotations

import awkward as ak
import numpy as np

__all__ = [
    "to_cartesian",
    "invariant_mass",
    "collection_mass",
    "delta_phi",
    "delta_r",
    "transverse_mass",
    "scalar_ht",
]


def to_cartesian(pt, eta, phi, mass=0.0):
    """(pt, eta, phi, mass) -> dict of (px, py, pz, e).

    Works elementwise, so ``pt`` etc. may be scalars, numpy arrays, or
    jagged awkward arrays of any depth.
    """
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    e = np.sqrt(px**2 + py**2 + pz**2 + mass**2)
    return {"px": px, "py": py, "pz": pz, "e": e}


def _mass_from_cartesian(px, py, pz, e):
    m2 = e**2 - (px**2 + py**2 + pz**2)
    # Rounding can push a massless object's m^2 a few ulp negative.
    m2 = np.where(m2 > 0.0, m2, 0.0)
    return np.sqrt(m2)


def invariant_mass(*objects):
    """Invariant mass of the vector sum of any number of objects.

    Each object is a mapping (or record array) with ``pt``, ``eta``, ``phi``
    and optionally ``mass``. Objects are summed elementwise, so this is for
    combining *different* objects (lepton + jet + jet), not for summing a
    jagged collection -- use :func:`collection_mass` for that.
    """
    if not objects:
        raise ValueError("invariant_mass needs at least one object")

    px = py = pz = e = None
    for obj in objects:
        c = to_cartesian(obj["pt"], obj["eta"], obj["phi"], _mass_of(obj))
        if px is None:
            px, py, pz, e = c["px"], c["py"], c["pz"], c["e"]
        else:
            px = px + c["px"]
            py = py + c["py"]
            pz = pz + c["pz"]
            e = e + c["e"]
    return _mass_from_cartesian(px, py, pz, e)


def collection_mass(collection, axis: int = 1):
    """Invariant mass of all objects in a jagged collection, per event.

    >>> # mass of all jets in each event
    >>> collection_mass(objs["jet"])
    """
    c = to_cartesian(
        collection["pt"], collection["eta"], collection["phi"], _mass_of(collection)
    )
    return _mass_from_cartesian(
        ak.sum(c["px"], axis=axis),
        ak.sum(c["py"], axis=axis),
        ak.sum(c["pz"], axis=axis),
        ak.sum(c["e"], axis=axis),
    )


def delta_phi(phi1, phi2):
    """Signed phi difference wrapped into (-pi, pi]."""
    d = phi1 - phi2
    return (d + np.pi) % (2 * np.pi) - np.pi


def delta_r(eta1, phi1, eta2, phi2):
    """Angular separation sqrt(deta^2 + dphi^2)."""
    return np.sqrt((eta1 - eta2) ** 2 + delta_phi(phi1, phi2) ** 2)


def transverse_mass(pt1, phi1, pt2, phi2):
    """Transverse mass of two massless objects, e.g. lepton and MET.

    mT = sqrt(2 pt1 pt2 (1 - cos dphi))
    """
    arg = 2.0 * pt1 * pt2 * (1.0 - np.cos(delta_phi(phi1, phi2)))
    arg = np.where(arg > 0.0, arg, 0.0)
    return np.sqrt(arg)


def scalar_ht(jet_pt, axis: int = 1):
    """Scalar pT sum per event.

    Note this is *your* HT over whatever jets you pass in, which is not the
    same thing as the Delphes ``ScalarHT.HT`` branch (that one includes
    leptons and photons). Compare the two as a sanity check, but do not
    expect them to agree.
    """
    return ak.sum(jet_pt, axis=axis)


def _mass_of(obj):
    """Mass field if present, else 0. Delphes leptons have no mass branch."""
    try:
        fields = ak.fields(obj)
    except (TypeError, AttributeError):
        fields = list(obj.keys()) if hasattr(obj, "keys") else []
    return obj["mass"] if "mass" in fields else 0.0
