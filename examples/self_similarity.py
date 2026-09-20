"""
Демонстрация φ-скейлинга FibonacciVortex.

Запуск:
    pip install matplotlib numpy
    python examples/self_similarity.py

Создаёт fibonacci_vortex.png
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fibonacci_vortex_core import (
    FibonacciVortex, PHI, INV_PHI,
    lamb_oseen_gamma, rankine_gamma,
)


def main():
    fv = FibonacciVortex([0, 0, 0], 13.0, 0.5, k_layers=5)

    rs = np.linspace(0.01, 5.0, 500)
    gammas = np.array([fv.gamma_of_r(r) for r in rs])
    sims = np.array([fv.self_similarity(r) for r in rs])
    lo = np.array([lamb_oseen_gamma(r, 13.0, 0.5) for r in rs])

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Γ(r) — сравнение
    axes[0].plot(rs, gammas, "b-", lw=2,
                  label="FibonacciVortex")
    axes[0].plot(rs, lo, "r--", lw=2,
                  label="Lamb–Oseen")
    axes[0].axhline(0, color="k", lw=0.5)
    axes[0].set_xlabel("r / r_c")
    axes[0].set_ylabel("Γ(r)")
    axes[0].set_title(
        r"FibonacciVortex vs Lamb–Oseen")
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    # Self-similarity
    axes[1].plot(rs, sims, "b-", lw=2,
                  label=r"$\Gamma(\varphi r)/\Gamma(r)$")
    axes[1].axhline(INV_PHI, color="k", ls="--",
                     label=f"1/φ = {INV_PHI:.4f}")
    axes[1].set_xlabel("r / r_c")
    axes[1].set_ylabel(r"$\Gamma(\varphi r)/\Gamma(r)$")
    axes[1].set_title(
        r"Self-similarity approaches $1/\varphi$")
    axes[1].legend()
    axes[1].grid(alpha=0.3)

    plt.tight_layout()
    out = "fibonacci_vortex.png"
    plt.savefig(out, dpi=120)
    plt.close()
    print(f"Сохранено: {out}")


if __name__ == "__main__":
    main()