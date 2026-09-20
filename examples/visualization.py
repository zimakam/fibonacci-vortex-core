"""
Visualization of FibonacciVortex.

Создаёт fibonacci_vortex.png с 4 графиками:
    1. Γ(r) — сравнение с Rankine и Lamb–Oseen
    2. Self-similarity Γ(φr)/Γ(r) vs 1/φ
    3. Tangential velocity u_θ(r)
    4. Vorticity ω(r)

Run:
    pip install matplotlib numpy
    python examples/visualization.py

Output:
    fibonacci_vortex.png
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
    # ─── Параметры ───
    G0 = 13.0
    rc = 0.5
    K = 5
    fv = FibonacciVortex([0, 0, 0], G0, rc, k_layers=K)

    # ─── Сетка r ───
    rs = np.linspace(0.01, 5.0 * rc, 500)

    # ─── Вычисления ───
    gamma_fib = np.array([fv.gamma_of_r(r) for r in rs])
    gamma_lo = np.array([lamb_oseen_gamma(r, G0, rc) for r in rs])
    gamma_rank = np.array([rankine_gamma(r, G0, rc) for r in rs])
    sims = np.array([fv.self_similarity(r) for r in rs])

    # Velocity u_θ = Γ/(2πr)
    u_theta = gamma_fib / (2 * np.pi * rs)

    # Vorticity (numerical)
    omega = np.array([
        fv.vorticity_magnitude([r, 0, 0]) for r in rs])

    # ─── Фигура: 2×2 ───
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        f"FibonacciVortex — φ-modulated vortex core "
        f"(Γ₀={G0}, r_c={rc}, K={K})",
        fontsize=14, fontweight="bold")

    # ─── 1. Γ(r) — сравнение ───
    ax = axes[0, 0]
    ax.plot(rs / rc, gamma_fib, "b-", lw=2.5,
             label="FibonacciVortex")
    ax.plot(rs / rc, gamma_lo, "r--", lw=2,
             label="Lamb–Oseen")
    ax.plot(rs / rc, gamma_rank, "g:", lw=2,
             label="Rankine")
    ax.axhline(G0, color="k", lw=0.5, alpha=0.3)
    ax.axhline(0, color="k", lw=0.5)
    ax.axvline(1.0, color="gray", lw=0.5, ls="--",
                alpha=0.5, label="r = r_c")
    ax.set_xlabel("r / r_c")
    ax.set_ylabel("Γ(r)")
    ax.set_title("Circulation Γ(r)")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.3)

    # ─── 2. Self-similarity ───
    ax = axes[0, 1]
    ax.plot(rs / rc, sims, "b-", lw=2.5,
             label=r"$\Gamma(\varphi r)/\Gamma(r)$")
    ax.axhline(INV_PHI, color="r", ls="--", lw=2,
                label=f"1/φ = {INV_PHI:.6f}")
    ax.axvline(1.0, color="gray", lw=0.5, ls="--",
                alpha=0.5)
    ax.set_xlabel("r / r_c")
    ax.set_ylabel(r"$\Gamma(\varphi r) / \Gamma(r)$")
    ax.set_title(r"Self-similarity approaches $1/\varphi$")
    ax.legend(loc="lower right")
    ax.grid(alpha=0.3)
    ax.set_ylim(0.4, 0.8)

    # ─── 3. Tangential velocity ───
    ax = axes[1, 0]
    ax.plot(rs / rc, u_theta, "b-", lw=2.5,
             label=r"$u_\theta(r) = \Gamma(r)/(2\pi r)$")
    # Пик скорости
    idx_max = int(np.argmax(u_theta))
    ax.axvline(rs[idx_max] / rc, color="r", ls=":",
                lw=1.5,
                label=f"max at r/r_c = {rs[idx_max]/rc:.2f}")
    ax.set_xlabel("r / r_c")
    ax.set_ylabel(r"$u_\theta(r)$")
    ax.set_title("Tangential velocity")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.3)

    # ─── 4. Vorticity ───
    ax = axes[1, 1]
    ax.plot(rs / rc, omega, "b-", lw=2.5,
             label=r"$\omega(r) \approx |d\Gamma/dr| / r$")
    ax.axvline(1.0, color="gray", lw=0.5, ls="--",
                alpha=0.5, label="r = r_c")
    ax.set_xlabel("r / r_c")
    ax.set_ylabel(r"$\omega(r)$")
    ax.set_title("Vorticity magnitude")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # ─── Сохранение ───
    out = "fibonacci_vortex.png"
    plt.savefig(out, dpi=130)
    plt.close()
    print(f"Сохранено: {out}")
    print(f"  Γ(0)   = {fv.gamma_of_r(0.0):.6f}")
    print(f"  Γ(r_c) = {fv.gamma_of_r(rc):.6f}")
    print(f"  self-similarity(r_c) = "
          f"{fv.self_similarity(rc):.6f}")
    print(f"  1/φ    = {INV_PHI:.6f}")


if __name__ == "__main__":
    main()