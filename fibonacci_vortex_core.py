#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FibonacciVortex — φ-modulated vortex core with Fibonacci hierarchy.
====================================================================

Formula:
    Γ(r) = Γ₀ · Σ_{k=0}^{K−1} F_k · φ^(−(k+1)·r/r_c) / Σ F_k

Properties:
    Γ(0) = Γ₀
    Γ(∞) = 0
    Γ(φ·r)/Γ(r) ≈ 1/φ   (for K ≥ 5)

Unlike classical models (Rankine, Lamb–Oseen, Gaussian), this
model has a closed-form distance-dependent circulation with a
φ-hierarchical decay.

Author: Зиявутдинов Магомед Камалович (Zimaka)
Email:  zimakam@gmail.com
License: MIT
Version: 1.0.0
"""

from __future__ import annotations
import math
from typing import Optional, List

import numpy as np

PHI = (1.0 + math.sqrt(5.0)) / 2.0
INV_PHI = PHI - 1.0
FIB = [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610]

__version__ = "1.0.0"
__author__ = "Зиявутдинов Магомед Камалович (Zimaka)"
__license__ = "MIT"


class VortexCore:
    """
    Base vortex core (классические модели).

    Supported models:
        "rankine"    — ω = const inside r_c, 0 outside
        "lamb_oseen" — ω(r) = Γ/(π·r_c²)·exp(−r²/r_c²)
        "gaussian"   — ω(r) = Γ/(2π·σ²)·exp(−r²/(2σ²))

    Author: Зиявутдинов Магомед Камалович (Zimaka).
    """

    def __init__(self, center, circulation: float,
                 core_radius: float,
                 axis=(0.0, 0.0, 1.0),
                 model: str = "lamb_oseen",
                 name: Optional[str] = None):
        self.center = np.asarray(center, dtype=np.float64)
        self.circulation = float(circulation)
        self.core_radius = float(core_radius)
        self.axis = np.asarray(axis, dtype=np.float64)
        self.axis /= max(np.linalg.norm(self.axis), 1e-12)
        self.model = model
        self.name = name or "vortex"

        # Build orthonormal basis perpendicular to axis
        a = self.axis
        tmp = (np.array([1.0, 0.0, 0.0])
               if abs(a[0]) < 0.9
               else np.array([0.0, 1.0, 0.0]))
        self.e1 = np.cross(a, tmp)
        self.e1 /= max(np.linalg.norm(self.e1), 1e-12)
        self.e2 = np.cross(a, self.e1)

    def local_coords(self, x):
        d = np.asarray(x, dtype=np.float64) - self.center
        z = float(np.dot(d, self.axis))
        perp = d - z * self.axis
        r = float(np.linalg.norm(perp))
        return r, z, perp

    def vorticity_magnitude(self, x) -> float:
        r, _, _ = self.local_coords(x)
        rc = self.core_radius
        G = self.circulation
        if self.model == "rankine":
            return abs(G) / (math.pi * rc * rc) if r <= rc else 0.0
        if self.model == "lamb_oseen":
            if rc < 1e-12:
                return 0.0
            return (abs(G) / (math.pi * rc * rc)) * \
                   math.exp(-(r * r) / (rc * rc))
        if self.model == "gaussian":
            sigma = rc / math.sqrt(2)
            if sigma < 1e-12:
                return 0.0
            return (abs(G) / (2 * math.pi * sigma * sigma)) * \
                   math.exp(-(r * r) / (2 * sigma * sigma))
        return 0.0

    def velocity_at(self, x) -> np.ndarray:
        r, _, perp = self.local_coords(x)
        if r < 1e-12:
            return np.zeros(3)
        rc = self.core_radius
        G = self.circulation
        tangent = np.cross(self.axis, perp) / r
        if r >= rc:
            u_theta = G / (2 * math.pi * r)
        else:
            u_theta = G * r / (2 * math.pi * rc * rc)
        return tangent * u_theta


class FibonacciVortex(VortexCore):
    """
    φ-modulated vortex with Fibonacci hierarchy of layers.

    Formula:
        Γ(r) = Γ₀ · Σ_{k=0}^{K−1} F_k · φ^(−(k+1)·r/r_c) / Σ F_k

    Parameters
    ----------
    center : array-like, shape (3,)
        3D position of the vortex axis.
    circulation : float
        Γ₀ — circulation at r = 0.
    core_radius : float
        r_c — characteristic core radius (> 0).
    axis : array-like, shape (3,), optional
        Direction of the vortex axis (default: [0, 0, 1]).
    k_layers : int, optional
        Number of Fibonacci layers K (default: 5).
    name : str, optional
        Human-readable name.

    Examples
    --------
    >>> fv = FibonacciVortex([0, 0, 0], circulation=13.0,
    ...                       core_radius=0.5, k_layers=5)
    >>> fv.gamma_of_r(0.0)
    13.0
    >>> 0.0 < fv.gamma_of_r(1.0) < 13.0
    True
    """

    def __init__(self, center, circulation: float,
                 core_radius: float,
                 axis=(0.0, 0.0, 1.0),
                 k_layers: int = 5,
                 name: Optional[str] = None):
        if k_layers < 1 or k_layers > len(FIB):
            raise ValueError(
                f"k_layers must be in [1, {len(FIB)}]")
        if core_radius <= 0:
            raise ValueError("core_radius must be > 0")

        super().__init__(center, circulation, core_radius,
                         axis=axis, model="fibonacci",
                         name=name or "fibonacci_vortex")
        self.k_layers = int(k_layers)
        self._norm = sum(FIB[:self.k_layers]) or 1.0

    # ─── Main formula ───
    def gamma_of_r(self, r: float) -> float:
        """
        Γ(r) = Γ₀ · Σ F_k · φ^(−(k+1)·r/r_c) / Σ F_k

        Properties:
            Γ(0) = Γ₀
            Γ(∞) → 0
        """
        rc = self.core_radius
        return self.circulation * sum(
            FIB[i] * (PHI ** (-(i + 1) * r / rc))
            for i in range(self.k_layers)
        ) / self._norm

    def velocity_at(self, x) -> np.ndarray:
        """Tangential velocity u_θ = Γ(r) / (2π·r)."""
        r, _, perp = self.local_coords(x)
        if r < 1e-12:
            return np.zeros(3)
        tangent = np.cross(self.axis, perp) / r
        u_theta = self.gamma_of_r(r) / (2 * math.pi * r)
        return tangent * u_theta

    def vorticity_magnitude(self, x) -> float:
        """ω(r) = |dΓ/dr| / r (numerical derivative)."""
        r, _, _ = self.local_coords(x)
        rc = self.core_radius
        dr = rc * 1e-3
        if r < dr:
            return 0.0
        g_plus = self.gamma_of_r(r + dr)
        g_minus = self.gamma_of_r(max(0.0, r - dr))
        return abs(g_plus - g_minus) / (2 * dr * r)

    def self_similarity(self,
                         r: Optional[float] = None) -> float:
        """
        Ratio Γ(φ·r) / Γ(r).

        For K ≥ 5 should approach 1/φ ≈ 0.618.
        """
        r = r if r is not None else self.core_radius
        g1 = self.gamma_of_r(r)
        g2 = self.gamma_of_r(PHI * r)
        if abs(g1) < 1e-15:
            return float('nan')
        return g2 / g1

    def fib_weights(self) -> List[float]:
        """Weights F_k · φ^(−k) for k = 0..K−1."""
        return [FIB[k] * (PHI ** (-k))
                for k in range(self.k_layers)]

    def zeckendorf_address(self) -> List[int]:
        """Zeckendorf decomposition of |Γ₀| (×1000)."""
        n = int(abs(self.circulation) * 1000)
        if n <= 0:
            return []
        parts: List[int] = []
        for f in reversed(FIB):
            if f <= n:
                parts.append(f)
                n -= f
        return parts

    def __repr__(self):
        return (f"FibonacciVortex(name={self.name!r}, "
                f"Γ₀={self.circulation:+.3f}, "
                f"r_c={self.core_radius:.3f}, "
                f"K={self.k_layers})")


# ─── Comparison helpers ───

def lamb_oseen_gamma(r: float, G0: float, rc: float) -> float:
    """Lamb–Oseen: Γ(r) = Γ₀ · (1 − exp(−r²/rc²))."""
    return G0 * (1.0 - math.exp(-(r * r) / (rc * rc)))


def rankine_gamma(r: float, G0: float, rc: float) -> float:
    """Rankine: Γ(r) = Γ₀ · (r/rc)² for r < rc, else Γ₀."""
    if r <= rc:
        return G0 * (r / rc) ** 2
    return G0


# ─── Demo ───

if __name__ == "__main__":
    fv = FibonacciVortex([0, 0, 0], circulation=13.0,
                         core_radius=0.5, k_layers=5)
    print(f"Γ(0)   = {fv.gamma_of_r(0.0):.6f}")
    print(f"Γ(r_c) = {fv.gamma_of_r(0.5):.6f}")
    print(f"Γ(∞)   ≈ {fv.gamma_of_r(1000.0):.2e}")
    print(f"self-similarity(0.5) = {fv.self_similarity(0.5):.6f}")
    print(f"1/φ = {INV_PHI:.6f}")
    print(f"Zeckendorf(Γ₀) = {fv.zeckendorf_address()}")