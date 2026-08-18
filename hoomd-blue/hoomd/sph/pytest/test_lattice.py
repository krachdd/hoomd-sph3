"""Kernel-consistency tests on a periodic cubic lattice.

On a uniform lattice with spacing dx and h = OptimalH * dx:
  * partition of unity: sum_j V_j W_ij ~= 1 (kernel completeness)
  * gradient consistency: sum_j V_j grad W_ij == 0 (by symmetry)
These are the discrete conditions that make the density summation and the
pressure-gradient operators zeroth/first-order consistent in the bulk.
"""

import numpy as np
import pytest

from .conftest import KERNELS

import hoomd.sph.kernel as skernel


@pytest.mark.parametrize("name", KERNELS)
def test_partition_of_unity_and_gradient(name):
    kern = skernel.Kernels[name]()
    dx = 0.1
    h = skernel.OptimalH[name] * dx
    kappa = kern.Kappa()
    rcut = kappa * h
    V = dx**3

    # neighbor offsets within the support (excluding self)
    nmax = int(np.ceil(rcut / dx))
    rng = np.arange(-nmax, nmax + 1)
    ii, jj, kk = np.meshgrid(rng, rng, rng, indexing="ij")
    offsets = np.column_stack([ii.ravel(), jj.ravel(), kk.ravel()]) * dx
    r = np.linalg.norm(offsets, axis=1)
    inside = (r > 1e-12) & (r < rcut)
    offsets = offsets[inside]
    r = r[inside]

    # partition of unity, including the self contribution W(0)
    wsum = kern.EvalKernel(h, 0.0) * V
    grad = np.zeros(3)
    for d, ri in zip(offsets, r):
        w = kern.EvalKernel(h, float(ri))
        dw = kern.EvalKernelDerivative(h, float(ri))
        wsum += V * w
        grad += V * dw / ri * d

    assert wsum == pytest.approx(1.0, rel=2e-2)
    assert np.allclose(grad, 0.0, atol=1e-10)
