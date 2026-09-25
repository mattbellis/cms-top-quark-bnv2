"""Unit tests for bnv2.kinematics, on hand-computable mock data.

Every expected value here is something you can check on paper, which is the
point -- these catch sign errors and eta/theta confusion, not physics.
"""

import awkward as ak
import numpy as np
import pytest

from bnv2 import kinematics as kin


# --------------------------------------------------------------------------
# mock data
# --------------------------------------------------------------------------

def make_obj(pt, eta, phi, mass=None):
    d = {"pt": np.asarray(pt, dtype=float),
         "eta": np.asarray(eta, dtype=float),
         "phi": np.asarray(phi, dtype=float)}
    if mass is not None:
        d["mass"] = np.asarray(mass, dtype=float)
    return d


@pytest.fixture
def two_jets():
    """Three events of two jets each, jagged, with known per-event sums."""
    return ak.zip(
        {
            "pt": ak.Array([[50.0, 50.0], [100.0, 20.0], [30.0, 30.0]]),
            "eta": ak.Array([[0.0, 0.0], [1.0, -1.0], [0.0, 0.0]]),
            "phi": ak.Array([[0.0, np.pi], [0.0, np.pi], [0.0, 0.0]]),
            "mass": ak.Array([[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]),
        }
    )


# --------------------------------------------------------------------------
# to_cartesian
# --------------------------------------------------------------------------

def test_to_cartesian_along_x():
    c = kin.to_cartesian(pt=10.0, eta=0.0, phi=0.0, mass=0.0)
    assert c["px"] == pytest.approx(10.0)
    assert c["py"] == pytest.approx(0.0)
    assert c["pz"] == pytest.approx(0.0)
    assert c["e"] == pytest.approx(10.0)


def test_to_cartesian_massive_object_energy():
    # pt=0 is not expressible in (pt,eta,phi); use eta=0 so p = pt.
    c = kin.to_cartesian(pt=3.0, eta=0.0, phi=np.pi / 2, mass=4.0)
    assert c["px"] == pytest.approx(0.0, abs=1e-12)
    assert c["py"] == pytest.approx(3.0)
    assert c["e"] == pytest.approx(5.0)  # 3-4-5


def test_to_cartesian_pz_uses_sinh_eta():
    c = kin.to_cartesian(pt=10.0, eta=1.0, phi=0.0)
    assert c["pz"] == pytest.approx(10.0 * np.sinh(1.0))


# --------------------------------------------------------------------------
# invariant_mass
# --------------------------------------------------------------------------

def test_invariant_mass_single_object_returns_its_own_mass():
    obj = make_obj([20.0], [0.7], [1.1], [5.0])
    assert kin.invariant_mass(obj)[0] == pytest.approx(5.0)


def test_invariant_mass_back_to_back_massless_pair():
    # two 50 GeV massless objects, eta=0, back to back -> m = 100
    a = make_obj([50.0], [0.0], [0.0])
    b = make_obj([50.0], [0.0], [np.pi])
    assert kin.invariant_mass(a, b)[0] == pytest.approx(100.0)


def test_invariant_mass_collinear_massless_pair_is_zero():
    a = make_obj([50.0], [0.3], [0.9])
    b = make_obj([20.0], [0.3], [0.9])
    assert kin.invariant_mass(a, b)[0] == pytest.approx(0.0, abs=1e-6)


def test_invariant_mass_three_objects():
    # three massless objects at 120 degrees in phi, eta=0, pt=10 each.
    # sum of 3-momenta vanishes, so m = sum(E) = 30.
    phis = [0.0, 2 * np.pi / 3, 4 * np.pi / 3]
    objs = [make_obj([10.0], [0.0], [p]) for p in phis]
    assert kin.invariant_mass(*objs)[0] == pytest.approx(30.0)


def test_invariant_mass_is_boost_invariant_along_z():
    """Same pair, both objects shifted by the same eta -> same mass."""
    a = make_obj([50.0], [0.0], [0.0])
    b = make_obj([70.0], [0.0], [1.2])
    m0 = kin.invariant_mass(a, b)[0]

    deta = 0.8
    a2 = make_obj([50.0], [deta], [0.0])
    b2 = make_obj([70.0], [deta], [1.2])
    assert kin.invariant_mass(a2, b2)[0] == pytest.approx(m0)


def test_invariant_mass_requires_arguments():
    with pytest.raises(ValueError):
        kin.invariant_mass()


def test_invariant_mass_missing_mass_field_treated_as_massless():
    a = make_obj([50.0], [0.0], [0.0])            # no "mass" key
    b = make_obj([50.0], [0.0], [np.pi], [0.0])   # explicit zero mass
    assert kin.invariant_mass(a, a)[0] == pytest.approx(0.0, abs=1e-6)
    assert kin.invariant_mass(a, b)[0] == pytest.approx(100.0)


# --------------------------------------------------------------------------
# collection_mass
# --------------------------------------------------------------------------

def test_collection_mass_per_event(two_jets):
    m = kin.collection_mass(two_jets)
    assert len(m) == 3
    assert m[0] == pytest.approx(100.0)   # back to back, 50+50
    assert m[2] == pytest.approx(0.0, abs=1e-6)  # collinear


def test_collection_mass_handles_empty_and_single_jet_events():
    jets = ak.zip(
        {
            "pt": ak.Array([[], [40.0]]),
            "eta": ak.Array([[], [0.0]]),
            "phi": ak.Array([[], [0.0]]),
            "mass": ak.Array([[], [0.0]]),
        }
    )
    m = kin.collection_mass(jets)
    assert m[0] == pytest.approx(0.0)
    assert m[1] == pytest.approx(0.0, abs=1e-6)


# --------------------------------------------------------------------------
# delta_phi / delta_r
# --------------------------------------------------------------------------

@pytest.mark.parametrize(
    "phi1, phi2, expected",
    [
        (0.1, -0.1, 0.2),
        (0.0, 0.0, 0.0),
        (3.0, -3.0, -(2 * np.pi - 6.0)),   # wraps the short way, negative
        (-3.0, 3.0, 2 * np.pi - 6.0),
        (np.pi, -np.pi, 0.0),
    ],
)
def test_delta_phi_wrapping(phi1, phi2, expected):
    assert kin.delta_phi(phi1, phi2) == pytest.approx(expected, abs=1e-9)


def test_delta_phi_always_in_range():
    rng = np.random.default_rng(42)
    phi1 = rng.uniform(-20, 20, 1000)
    phi2 = rng.uniform(-20, 20, 1000)
    d = kin.delta_phi(phi1, phi2)
    assert np.all(d > -np.pi - 1e-12)
    assert np.all(d <= np.pi + 1e-12)


def test_delta_r_pythagorean():
    # deta = 4, dphi = 3 (no wrap, |3| < pi) -> 3-4-5
    assert kin.delta_r(0.0, 0.0, 4.0, 3.0) == pytest.approx(5.0)


def test_delta_r_uses_wrapped_phi():
    # phi separation of 0.2 the short way, not 2pi - 0.2
    assert kin.delta_r(0.0, 3.1, 0.0, -3.1) == pytest.approx(2 * np.pi - 6.2)


def test_delta_r_on_jagged_arrays():
    eta = ak.Array([[0.0, 1.0], [2.0]])
    phi = ak.Array([[0.0, 0.0], [0.0]])
    dr = kin.delta_r(eta, phi, 0.0, 0.0)
    assert ak.to_list(ak.flatten(dr)) == pytest.approx([0.0, 1.0, 2.0])


# --------------------------------------------------------------------------
# transverse_mass
# --------------------------------------------------------------------------

def test_transverse_mass_back_to_back():
    # mT = sqrt(2 * 40 * 40 * (1 - cos(pi))) = 80
    assert kin.transverse_mass(40.0, 0.0, 40.0, np.pi) == pytest.approx(80.0)


def test_transverse_mass_collinear_is_zero():
    assert kin.transverse_mass(40.0, 0.7, 25.0, 0.7) == pytest.approx(0.0, abs=1e-9)


def test_transverse_mass_symmetric_in_its_arguments():
    a = kin.transverse_mass(31.0, 0.3, 47.0, 2.2)
    b = kin.transverse_mass(47.0, 2.2, 31.0, 0.3)
    assert a == pytest.approx(b)


# --------------------------------------------------------------------------
# scalar_ht
# --------------------------------------------------------------------------

def test_scalar_ht_per_event(two_jets):
    ht = kin.scalar_ht(two_jets.pt)
    assert ak.to_list(ht) == pytest.approx([100.0, 120.0, 60.0])


def test_scalar_ht_empty_event_is_zero():
    ht = kin.scalar_ht(ak.Array([[], [10.0, 20.0]]))
    assert ak.to_list(ht) == pytest.approx([0.0, 30.0])
