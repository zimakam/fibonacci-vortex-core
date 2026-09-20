# fibonacci-vortex-core
φ-modulated vortex core with Fibonacci hierarchy of layers. Γ(r) = Γ₀ · Σ F_k · φ^(−(k+1)·r/r_c) / Σ F_k. Properties: Γ(0)=Γ₀, Γ(∞)=0, Γ(φr)/Γ(r)≈1/φ. Research prototype (TRL 3). MIT.
# fibonacci-vortex-core

**φ-modulated vortex core with Fibonacci hierarchy of layers.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/downloads/)

Author: **Зиявутдинов Магомед Камалович (Zimaka)** · zimakam@gmail.com

---

## Formula

```
Γ(r) = Γ₀ · Σ_{k=0}^{K−1} F_k · φ^(−(k+1)·r/r_c) / Σ F_k
```

where:
- `Γ₀` — circulation at `r = 0`
- `r_c` — characteristic core radius
- `F_k` — Fibonacci numbers
- `φ = (1 + √5)/2` — golden ratio
- `K` — number of layers (default: 5)

## Properties

| Property | Value |
|---|---|
| `Γ(0)` | `Γ₀` (exactly, for any `K ≥ 1`) |
| `Γ(∞)` | `0` (exponentially) |
| `Γ(φ·r)/Γ(r)` | `≈ 1/φ` (for `K ≥ 5`) |

## Compared to classical models

| Model | Formula | Character |
|---|---|---|
| Rankine | `ω = const` inside `r_c`, else 0 | hard core |
| Lamb–Oseen | `ω(r) = (Γ/πr_c²)·exp(−r²/r_c²)` | diffusion |
| Gaussian | `ω(r) = (Γ/2πσ²)·exp(−r²/2σ²)` | smooth core |
| **FibonacciVortex** | `Γ(r)` — φ-hierarchical sum | self-similar decay |

Unlike the classical models, here circulation is a
distance-dependent Fibonacci-weighted φ-exponential sum.

## Quick start

```bash
pip install numpy
```

```python
from fibonacci_vortex_core import FibonacciVortex

fv = FibonacciVortex(
    center=[0, 0, 0],
    circulation=13.0,   # Γ₀
    core_radius=0.5,    # r_c
    k_layers=5,         # K
)

print(fv.gamma_of_r(0.0))         # 13.0
print(fv.gamma_of_r(0.5))         # intermediate
print(fv.self_similarity(0.5))    # ≈ 0.618 (1/φ)
```

## Tests

```bash
pip install pytest
pytest tests/ -v
```

Covers:
- `Γ(0) = Γ₀` for any `K`
- monotonic decay
- `Γ(∞) = 0`
- φ-scaling `Γ(φr)/Γ(r) ≈ 1/φ`
- velocity / vorticity
- comparison with Rankine and Lamb–Oseen
- axis orientation
- Zeckendorf address

## Visualisation

```bash
pip install matplotlib
python examples/self_similarity.py
```

Creates `fibonacci_vortex.png` with:
- `Γ(r)` curve
- self-similarity ratio `Γ(φr)/Γ(r)` vs `r`

## Repository structure

```
fibonacci-vortex-core/
├── fibonacci_vortex_core.py        # module
├── tests/
│   └── test_fibonacci_vortex.py    # pytest
├── examples/
│   └── self_similarity.py          # matplotlib demo
├── README.md                       # this file
├── THEORY.md                       # formulas and proofs
├── LICENSE                         # MIT
├── CITATION.cff                    # machine-readable citation
├── .zenodo.json                    # Zenodo metadata
├── pyproject.toml                  # package config
└── .gitignore
```

## Citation

```bibtex
@software{ziyavutdinov_fibonacci_vortex_core_2026,
  author    = {Зиявутдинов, Магомед Камалович},
  title     = {fibonacci-vortex-core: φ-modulated vortex core
               with Fibonacci hierarchy},
  year      = {2026},
  version   = {1.0.0},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.XXXXXXX},
  license   = {MIT}
}
```

## Status

Research prototype (**TRL 3**). Not a validated physical theory.
Numerical properties (`Γ(0)=Γ₀`, φ-scaling) are verified in tests.

## Related projects

- **Quantum Master V2** — multi-layer system with φ-Chain and HyperNav
- **ΔM-Antigravity** — ΔM antigravity model

## License

MIT — see [LICENSE](LICENSE).