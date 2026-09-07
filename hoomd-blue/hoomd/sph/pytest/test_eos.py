"""Unit tests for the SPH equations of state."""

import pytest

import hoomd.sph as sph


@pytest.fixture(params=["Tait", "Linear"])
def eos(request):
    e = getattr(sph.eos, request.param)()
    e.set_params(1000.0, 0.01)   # rho0, background-pressure factor
    e.set_speedofsound(15.0)
    return e


@pytest.mark.parametrize("rho", [990.0, 1000.0, 1005.0, 1020.0])
def test_pressure_density_roundtrip(eos, rho):
    p = eos.cpp_stateequation.Pressure(rho)
    rho_back = eos.cpp_stateequation.Density(p)
    assert rho_back == pytest.approx(rho, rel=1e-10)


def test_pressure_at_rest_density_is_background(eos):
    # p(rho0) = bp = bpfactor * rho0 * c^2
    p0 = eos.cpp_stateequation.Pressure(1000.0)
    assert p0 == pytest.approx(0.01 * 1000.0 * 15.0**2, rel=1e-12)


@pytest.mark.parametrize("rho", [995.0, 1000.0, 1010.0])
def test_dpressure_ddensity_matches_finite_difference(eos, rho):
    eps = 1e-3
    fd = (eos.cpp_stateequation.Pressure(rho + eps)
          - eos.cpp_stateequation.Pressure(rho - eps)) / (2 * eps)
    an = eos.cpp_stateequation.dPressuredDensity(rho)
    assert an == pytest.approx(fd, rel=1e-6)


def test_pressure_monotonically_increasing(eos):
    prev = eos.cpp_stateequation.Pressure(950.0)
    for rho in [980.0, 1000.0, 1020.0, 1050.0]:
        p = eos.cpp_stateequation.Pressure(rho)
        assert p > prev
        prev = p
