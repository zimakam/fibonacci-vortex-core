"""
Tests for fibonacci_vortex_core.

Run:
    pytest tests/ -v
"""
import math
import pytest
import numpy as np

from fibonacci_vortex_core import (
    FibonacciVortex, VortexCore, PHI, INV_PHI, FIB,
    lamb_oseen_gamma, rankine_gamma,
)


class TestFibonacciVortexBasics:
    """Basic properties of the φ-modulated vortex."""

    def test_gamma_at_zero(self):
        """Γ(0) = Γ₀ for any K."""
        for K in (1, 3, 5, 10):
            fv = FibonacciVortex([0, 0, 0], circulation=13.0,
                                  core_radius=0.5, k_layers=K)
            assert abs(fv.gamma_of_r(0.0) - 13.0) < 1e-9

    def test_gamma_decays(self):
        """Γ(r) is strictly decreasing for r > 0."""
        fv = FibonacciVortex([0, 0, 0], 13.0, 0.5, k_layers=5)
        rs = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0]
        gammas = [fv.gamma_of_r(r) for r in rs]
        for i in range(len(gammas) - 1):
            assert gammas[i] > gammas[i + 1]

    def test_gamma_at_infinity(self):
        """Γ(∞) = 0 (exponentially decays)."""
        fv = FibonacciVortex([0, 0, 0], 13.0, 0.5, k_layers=5)
        assert abs(fv.gamma_of_r(1000.0)) < 1e-9

    def test_k_layers_effect(self):
        """Different K → different values in between."""
        fv1 = FibonacciVortex([0, 0, 0], 13.0, 0.5, k_layers=1)
        fv5 = FibonacciVortex([0, 0, 0], 13.0, 0.5, k_layers=5)
        assert abs(fv1.gamma_of_r(0.5) - fv5.gamma_of_r(0.5)) > 1e-6

    def test_invalid_params(self):
        with pytest.raises(ValueError):
            FibonacciVortex([0, 0, 0], 13.0, 0.5, k_layers=0)
        with pytest.raises(ValueError):
            FibonacciVortex([0, 0, 0], 13.0, -1.0)
        with pytest.raises(ValueError):
            FibonacciVortex([0, 0, 0], 13.0, 0.5, k_layers=100)


class TestSelfSimilarity:
    """φ-scaling properties."""

    def test_self_similarity_phi(self):
        """Γ(φ·r)/Γ(r) ≈ 1/φ for K ≥ 5."""
        fv = FibonacciVortex([0, 0, 0], 13.0, 0.5, k_layers=5)
        s = fv.self_similarity(0.5)
        assert abs(s - INV_PHI) < 0.15

    def test_self_similarity_range(self):
        """Ratio is always positive and less than 1."""
        fv = FibonacciVortex([0, 0, 0], 13.0, 0.5, k_layers=5)
        for r in (0.1, 0.5, 1.0, 2.0):
            s = fv.self_similarity(r)
            assert 0.0 < s < 1.0


class TestVelocity:
    """Velocity field properties."""

    def test_velocity_tangent(self):
        """Velocity is perpendicular to the axis."""
        fv = FibonacciVortex([0, 0, 0], 13.0, 0.5)
        u = fv.velocity_at([1.0, 0.0, 0.0])
        # axis = z, so u_z should be ~0
        assert abs(u[2]) < 1e-12
        assert math.sqrt(u[0]**2 + u[1]**2) > 0

    def test_velocity_at_center(self):
        """Velocity at the axis is zero."""
        fv = FibonacciVortex([0, 0, 0], 13.0, 0.5)
        u = fv.velocity_at([0.0, 0.0, 0.0])
        assert abs(u[0]) < 1e-12
        assert abs(u[1]) < 1e-12

    def test_vorticity_positive(self):
        """Vorticity is non-negative."""
        fv = FibonacciVortex([0, 0, 0], 13.0, 0.5)
        for r in (0.1, 0.5, 1.0, 2.0):
            assert fv.vorticity_magnitude([r, 0, 0]) >= 0


class TestComparison:
    """Comparison with classical models."""

    def test_lamb_oseen_different(self):
        """Fibonacci and Lamb-Oseen have different forms."""
        fv = FibonacciVortex([0, 0, 0], 13.0, 0.5, k_layers=5)
        lo = lamb_oseen_gamma(0.5, 13.0, 0.5)
        fib = fv.gamma_of_r(0.5)
        assert abs(lo - fib) > 0.1

    def test_rankine_flat_outside(self):
        """Rankine is constant Γ₀ outside the core."""
        g1 = rankine_gamma(0.5, 13.0, 0.5)
        g2 = rankine_gamma(2.0, 13.0, 0.5)
        assert abs(g1 - g2) < 1e-12

    def test_rankine_grows_inside(self):
        """Rankine Γ ∝ r² inside the core."""
        g1 = rankine_gamma(0.25, 13.0, 0.5)
        g2 = rankine_gamma(0.5, 13.0, 0.5)
        assert g2 > g1


class TestZeckendorf:
    """Zeckendorf address of Γ₀."""

    def test_zeckendorf_13(self):
        fv = FibonacciVortex([0, 0, 0], 13.0, 0.5)
        # 13 × 1000 = 13000 = sum of Zeckendorf decomposition
        addr = fv.zeckendorf_address()
        assert sum(addr) == 13000

    def test_zeckendorf_zero(self):
        fv = FibonacciVortex([0, 0, 0], 0.0, 0.5)
        assert fv.zeckendorf_address() == []


class TestAxisOrientation:
    """Vortex axis invariance."""

    def test_axis_z(self):
        fv = FibonacciVortex([0, 0, 0], 13.0, 0.5,
                              axis=(0, 0, 1))
        assert abs(np.linalg.norm(fv.axis) - 1.0) < 1e-12

    def test_axis_x(self):
        fv = FibonacciVortex([0, 0, 0], 13.0, 0.5,
                              axis=(1, 0, 0))
        u = fv.velocity_at([0.0, 1.0, 0.0])
        # axis = x, so u_x should be ~0
        assert abs(u[0]) < 1e-12


if __name__ == "__main__":
    pytest.main([__file__, "-v"])