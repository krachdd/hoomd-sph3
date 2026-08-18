"""Unit tests for the SPH smoothing kernels.

Checks, for every registered kernel:
  * 3D normalization: 4 pi int_0^{kappa h} W(r) r^2 dr == 1
  * consistency of the self-density value W(0) with w0(h)
  * compact support: W == 0 and dW/dr == 0 outside kappa h
  * the analytic derivative matches a central finite difference
  * dW/dr <= 0 everywhere (monotonically decreasing kernels)
"""

import numpy as np
import pytest

from .conftest import KERNELS

import hoomd.sph.kernel as skernel


@pytest.fixture(params=KERNELS)
def kernel(request):
    return skernel.Kernels[request.param]()


@pytest.mark.parametrize("h", [0.5, 1.0, 0.0085])
def test_normalization(kernel, h):
    """4 pi int W r^2 dr over the support must equal 1."""
    kappa = kernel.Kappa()
    r = np.linspace(0.0, kappa * h, 20001)
    w = np.array([kernel.EvalKernel(h, ri) for ri in r])
    integral = 4.0 * np.pi * np.trapezoid(w * r**2, r)
    assert integral == pytest.approx(1.0, rel=1e-4)


def test_self_density_consistency(kernel):
    h = 0.7
    assert kernel.cpp_smoothingkernel.w0(h) == pytest.approx(
        kernel.EvalKernel(h, 0.0), rel=1e-12)


def test_compact_support(kernel):
    h = 1.0
    kappa = kernel.Kappa()
    for r in [kappa * h * 1.0001, kappa * h * 2.0, 10.0 * h]:
        assert kernel.EvalKernel(h, r) == 0.0
        assert kernel.EvalKernelDerivative(h, r) == 0.0


def test_derivative_matches_finite_difference(kernel):
    h = 1.0
    kappa = kernel.Kappa()
    eps = 1e-6
    # stay away from the piecewise-spline knots of the quintic kernel
    for q in [0.11, 0.42, 0.77, 1.23, 1.71]:
        r = q * h * kappa / 2.0  # sample within the support for kappa=2 and 3
        if r + eps >= kappa * h:
            continue
        fd = (kernel.EvalKernel(h, r + eps)
              - kernel.EvalKernel(h, r - eps)) / (2 * eps)
        an = kernel.EvalKernelDerivative(h, r)
        assert an == pytest.approx(fd, rel=1e-4, abs=1e-8)


def test_derivative_nonpositive(kernel):
    h = 1.0
    kappa = kernel.Kappa()
    r = np.linspace(0.0, kappa * h, 500)
    dw = np.array([kernel.EvalKernelDerivative(h, ri) for ri in r])
    assert np.all(dw <= 1e-14)
