# Theoretical notes on FibonacciVortex

**Author:** Зиявутдинов Магомед Камалович (Zimaka)

This document collects the mathematical formulation, provable
properties, and honest limitations of the FibonacciVortex model.

---

## 1. Formula

```
Γ(r) = Γ₀ · Σ_{k=0}^{K−1} F_k · φ^(−(k+1)·r/r_c) / Σ F_k
```

where:
- `Γ₀` — circulation at `r = 0`
- `r_c > 0` — characteristic core radius
- `F_k` — Fibonacci numbers: `F_0 = 1, F_1 = 1, F_2 = 2, F_3 = 3, ...`
- `φ = (1 + √5)/2 ≈ 1.618...` — golden ratio
- `K ≥ 1` — number of layers

## 2. Provable properties

### Theorem 1. `Γ(0) = Γ₀`

*Proof.* At `r = 0`, `φ^(−(k+1)·0) = 1` for every `k`. Hence

```
Γ(0) = Γ₀ · Σ F_k / Σ F_k = Γ₀.  ∎
```

### Theorem 2. `Γ(r) → 0` as `r → ∞`

*Proof.* Every term `F_k · φ^(−(k+1)·r/r_c)` decays exponentially
as `r → ∞` since `φ > 1` and `(k+1) ≥ 1`. Finite sum of decaying
terms → 0.  ∎

### Theorem 3. Monotonicity

`Γ(r)` is strictly decreasing on `[0, ∞)` for `Γ₀ > 0`.

*Proof.* Each term is strictly decreasing, all terms have the same
sign (positive), and weights `F_k > 0`. Sum of strictly decreasing
nonnegative terms is strictly decreasing.  ∎

### Observation (φ-scaling)

For `K ≥ 5`, in a wide range of `r`:

```
Γ(φ·r) / Γ(r) ≈ 1/φ
```

**This is not a theorem** — it is a numerical observation confirmed
in `tests/test_fibonacci_vortex.py::TestSelfSimilarity`. The
approximation quality depends on `K`, `r_c`, and `r/r_c`.

## 3. Velocity and vorticity

### Tangential velocity

```
u_θ(r) = Γ(r) / (2π·r)
```

### Vorticity

For a 2D axisymmetric field: `ω = (1/r)·d(r·u_θ)/dr`. We use the
numerical approximation:

```
ω(r) ≈ |dΓ/dr| / r
```

computed by finite difference with `dr = 10⁻³·r_c`.

## 4. Comparison with classical models

### Rankine (1865)

```
ω(r) = Γ/(π·r_c²)   for r < r_c
ω(r) = 0            for r > r_c
```

Sharp cutoff. Discontinuous vorticity.

### Lamb–Oseen (1911)

```
ω(r) = Γ/(π·r_c²) · exp(−r²/r_c²)
```

Smooth Gaussian-like core. Derived from viscous diffusion of a
point vortex.

### Gaussian

```
ω(r) = Γ/(2π·σ²) · exp(−r²/2σ²)
```

Standard Gaussian profile with variance `σ²`.

### FibonacciVortex

Distance-dependent circulation through a φ-weighted Fibonacci
sum. No differential equation derivation — the model is
phenomenological.

## 5. Why Fibonacci + φ?

The motivation is a hierarchical decay of circulation at
multiple scales. Each layer `k` decays with rate `φ^(k+1)` and
weight `F_k`. This produces:

- a **self-similar** structure (`r → φ·r` ≈ `1/φ` decay)
- a **closed-form** circulation for any `r`
- controllable complexity via `K`

## 6. Physical interpretation (cautious)

FibonacciVortex can be viewed as a **model** for vortices with
multi-scale structure (e.g., turbulent cores, plasma vortices).
The φ-modulation is a mathematical construct, not a physical law.

## 7. What is NOT claimed

- ✗ That FibonacciVortex matches real measured vortices
- ✗ That `φ`-scaling is a fundamental property of nature
- ✗ That this model replaces Rankine / Lamb–Oseen / Gaussian
- ✗ That the formula follows from any physical equations

## 8. What is claimed

- ✓ `Γ(0) = Γ₀` (exact)
- ✓ `Γ(∞) = 0` (exact)
- ✓ Monotonicity (exact)
- ✓ Numerical φ-scaling for `K ≥ 5` (approximate, verified)
- ✓ Model is a valid closed-form construction

## 9. Numerical verification

| Property | Test |
|---|---|
| `Γ(0) = Γ₀` | `test_gamma_at_zero` |
| Monotonicity | `test_gamma_decays` |
| `Γ(∞) = 0` | `test_gamma_at_infinity` |
| φ-scaling | `test_self_similarity_phi` |
| K-dependence | `test_k_layers_effect` |
| Velocity ⟂ axis | `test_velocity_tangent` |
| Vorticity ≥ 0 | `test_vorticity_positive` |
| Comparison with Lamb–Oseen | `test_lamb_oseen_different` |
| Zeckendorf address | `test_zeckendorf_13` |
| Axis orientation | `test_axis_x`, `test_axis_z` |

## 10. Reproducibility

Environment: Python 3.9+, NumPy ≥ 1.21.
Random seeds not used in model itself — properties are
deterministic and reproducible exactly.

## 11. References

- Rankine, W. J. M. (1865). *On the thermodynamic theory of waves
  of finite longitudinal disturbance.*
- Lamb, H. (1932). *Hydrodynamics* (6th ed.). Cambridge University
  Press.
- Oseen, C. W. (1911). *Über die Wirbelbewegung in einer reibenden
  Flüssigkeit.* Ark. Mat. Astro. Fys. 7.
- Saffman, P. G. (1992). *Vortex Dynamics.* Cambridge University
  Press.
- Kauffman, L. H. (2001). *Knots and Physics* (3rd ed.).
  World Scientific.

## 12. Contact

**Зиявутдинов Магомед Камалович (Zimaka)**
zimakam@gmail.com