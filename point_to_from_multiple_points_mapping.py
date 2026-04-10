"""
Shamir's Secret Sharing — Demo
================================
Demo 1: Generate a random seed, SHA-256 hash it, split with Shamir SSS,
        plot the (x, y) share points and verify recovery.

Demo 2: Given N random (x, y) points, use Lagrange interpolation (same math
        as Shamir reconstruction) to combine them into a single hash at x=0.

All arithmetic is performed in GF(p) where p = 2^127 - 1 (Mersenne prime),
so reconstruction is exact — no floating point.
"""

import hashlib
import secrets
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ---------------------------------------------------------------------------
# Field arithmetic over GF(PRIME)
# ---------------------------------------------------------------------------

PRIME = 2**127 - 1  # 12th Mersenne prime


def mod_inverse(k: int, p: int = PRIME) -> int:
    """Modular inverse via Fermat's little theorem."""
    return pow(k, p - 2, p)


def horner(x: int, coeffs: list[int], p: int = PRIME) -> int:
    """Evaluate polynomial at x using Horner's method."""
    result = 0
    for c in reversed(coeffs):
        result = (result * x + c) % p
    return result


# ---------------------------------------------------------------------------
# Shamir Secret Sharing
# ---------------------------------------------------------------------------

def make_shares(secret_int: int, k: int, n: int, p: int = PRIME) -> list[tuple[int, int]]:
    """
    Split secret_int into n shares with threshold k.

    Constructs a random degree-(k-1) polynomial f over GF(p) with
    f(0) = secret_int, then returns [(1, f(1)), ..., (n, f(n))].
    """
    if not (0 <= secret_int < p):
        raise ValueError("secret_int must be in [0, p)")
    if k < 2:
        raise ValueError("Threshold k must be >= 2")
    if n < k:
        raise ValueError("Total shares n must be >= threshold k")

    coeffs = [secret_int] + [secrets.randbelow(p) for _ in range(k - 1)]
    return [(x, horner(x, coeffs, p)) for x in range(1, n + 1)]


def recover_secret(shares: list[tuple[int, int]], p: int = PRIME) -> int:
    """
    Recover the secret from k shares using Lagrange interpolation at x=0.

    Requires exactly k (threshold) shares — more is fine, fewer will give
    a wrong answer silently (by design of the scheme).
    """
    x_s = [s[0] for s in shares]
    y_s = [s[1] for s in shares]
    total = 0
    for i in range(len(x_s)):
        num = den = 1
        for j in range(len(x_s)):
            if i != j:
                num = (num * (0 - x_s[j])) % p
                den = (den * (x_s[i] - x_s[j])) % p
        total = (total + y_s[i] * num * mod_inverse(den, p)) % p
    return total


# ---------------------------------------------------------------------------
# Demo 1: Seed → SHA-256 → Shamir shares → plot → verify recovery
# ---------------------------------------------------------------------------

def demo1_shamir_plot(k: int = 2, n: int = 8) -> dict:
    """
    Generate a random seed, hash it, split into Shamir shares, plot, verify.

    Returns a dict with all intermediate values.
    """
    # 1. Random seed & hash
    seed = secrets.token_hex(32)
    hashed = hashlib.sha256(seed.encode()).hexdigest()
    secret_int = int(hashed, 16) % PRIME

    # 2. Split
    shares = make_shares(secret_int, k, n)

    # 3. Recover with exactly k shares
    recovered = recover_secret(shares[:k])
    assert recovered == secret_int, "Recovery failed!"

    return {
        "seed": seed,
        "sha256": hashed,
        "secret_int": secret_int,
        "k": k,
        "n": n,
        "shares": shares,
        "recovered": recovered,
        "verified": recovered == secret_int,
    }


# ---------------------------------------------------------------------------
# Demo 2: Random (x, y) points → Lagrange → single combined value at x=0
# ---------------------------------------------------------------------------

def demo2_combine_points(n: int = 6) -> dict:
    """
    Generate n random (x, y) points and combine them via Lagrange interpolation.

    The n points define a unique degree-(n-1) polynomial over GF(p).
    Evaluating at x=0 produces a deterministic "combined hash".
    """
    xs = list(range(1, n + 1))
    ys = [secrets.randbelow(PRIME) for _ in range(n)]
    shares = list(zip(xs, ys))
    combined = recover_secret(shares)
    return {
        "n": n,
        "points": shares,
        "combined_int": combined,
        "combined_hex": hex(combined)[2:].zfill(32),
    }


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

VSCALE = 10**14  # scale for visual display only


def build_plot(d1: dict, d2: dict, output_path: str = "shamir_demo.png") -> None:
    """Render both demos into a single figure and save to output_path."""
    fig = plt.figure(figsize=(18, 10), facecolor='#08080f')
    gs = GridSpec(2, 3, figure=fig, hspace=0.48, wspace=0.32,
                  left=0.06, right=0.97, top=0.92, bottom=0.07)

    ax1 = fig.add_subplot(gs[0, :2])
    ax2 = fig.add_subplot(gs[0, 2])
    ax3 = fig.add_subplot(gs[1, :])

    _plot_demo1_shares(ax1, d1)
    _plot_demo1_info(ax2, d1)
    _plot_demo2_combine(ax3, d2)

    fig.suptitle("Shamir's Secret Sharing  —  Polynomial Interpolation over GF(p)",
                 color='white', fontsize=15, fontweight='bold')
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='#08080f')
    print(f"Plot saved → {output_path}")


def _plot_demo1_shares(ax, d1: dict) -> None:
    shares = d1["shares"]
    k, n = d1["k"], d1["n"]

    xs_plot = [s[0] for s in shares]
    ys_plot = [s[1] % VSCALE for s in shares]

    # Visual line: use first 2 shares for slope (k=2 is linear)
    x0, y0 = shares[0][0], shares[0][1] % VSCALE
    x1, y1 = shares[1][0], shares[1][1] % VSCALE
    slope = (y1 - y0) / (x1 - x0)
    secret_vis = y0 - slope * x0

    x_line = np.linspace(-0.3, n + 0.5, 300)
    y_line = secret_vis + slope * x_line
    ax.plot(x_line, y_line, '--', color='#7f5af0', alpha=0.55, lw=1.8,
            label=f'Degree-{k-1} polynomial (k={k})')

    colors = plt.cm.plasma(np.linspace(0.15, 0.95, n))
    for i in range(n):
        ax.scatter(xs_plot[i], ys_plot[i], s=200, color=colors[i], zorder=5,
                   edgecolors='white', linewidths=1.1)
        ax.annotate(f' S{xs_plot[i]}', (xs_plot[i], ys_plot[i]),
                    fontsize=8.5, color=colors[i], fontfamily='monospace')

    ax.scatter(0, secret_vis, s=320, marker='*', color='#2cb67d', zorder=6,
               edgecolors='white', linewidths=1.5, label='Secret @ x=0 (recovered)')

    ax.set_facecolor('#0d0d20')
    ax.set_title(f'Demo 1 — {n} Shamir Shares (k={k} threshold)  |  '
                 f'Each share = point on polynomial', color='white', fontsize=11, pad=10)
    ax.set_xlabel('x  (share index)', color='#999')
    ax.set_ylabel('y  (value mod 10¹⁴, visual scale)', color='#999')
    ax.tick_params(colors='#777')
    ax.spines[:].set_color('#2a2a4a')
    ax.legend(facecolor='#15153a', labelcolor='white', fontsize=8.5, framealpha=0.8)
    ax.grid(alpha=0.08, color='white')


def _plot_demo1_info(ax, d1: dict) -> None:
    ax.set_facecolor('#0d0d20')
    ax.axis('off')
    lines = [
        ("SEED",       d1["seed"][:20] + "..."),
        ("SHA-256",    d1["sha256"][:20] + "..."),
        ("SECRET INT", str(d1["secret_int"])[:18] + "..."),
        ("SCHEME",     f"{d1['k']}-of-{d1['n']}  (k={d1['k']}, n={d1['n']})"),
        ("ANY k SHARES", "→ recover secret"),
        ("RECOVERY",   "✓  VERIFIED" if d1["verified"] else "✗  FAILED"),
    ]
    y_pos = 0.93
    for label, val in lines:
        ax.text(0.05, y_pos, label, transform=ax.transAxes, color='#7f5af0',
                fontsize=7.5, fontfamily='monospace', fontweight='bold')
        ax.text(0.05, y_pos - 0.065, val, transform=ax.transAxes, color='#c8ffa0',
                fontsize=7.5, fontfamily='monospace')
        y_pos -= 0.155
    ax.set_title('Parameters', color='white', fontsize=10, pad=8)
    ax.spines[:].set_color('#2a2a4a')


def _plot_demo2_combine(ax, d2: dict) -> None:
    n = d2["n"]
    xs_f = np.array([p[0] for p in d2["points"]], dtype=float)
    ys_f = np.array([p[1] % VSCALE for p in d2["points"]], dtype=float)

    coeffs_np = np.polyfit(xs_f, ys_f, n - 1)
    x_fine = np.linspace(0, n + 0.3, 500)
    y_fit = np.polyval(coeffs_np, x_fine)
    combined_vis = float(np.polyval(coeffs_np, 0))

    colors2 = plt.cm.cool(np.linspace(0.05, 0.95, n))
    ax.set_facecolor('#0d0d20')
    ax.plot(x_fine, y_fit, '-', color='#ff6b6b', alpha=0.5, lw=1.8,
            label=f'Degree-{n-1} polynomial through all {n} points')

    for i in range(n):
        ax.scatter(xs_f[i], ys_f[i], s=210, color=colors2[i], zorder=5,
                   edgecolors='white', linewidths=1.1)
        ax.annotate(f' P{int(xs_f[i])}', (xs_f[i], ys_f[i]),
                    fontsize=8.5, color=colors2[i], fontfamily='monospace')

    ax.scatter(0, combined_vis, s=370, marker='D', color='#ffd700', zorder=6,
               edgecolors='white', linewidths=1.5,
               label=f"Combined @ x=0  →  0x{d2['combined_hex'][:16]}...")

    ax.set_title(f'Demo 2 — {n} Random (x, y) Points → Lagrange Interpolation → Single Value at x=0\n'
                 f'Any unique set of {n} points defines a unique degree-{n-1} polynomial',
                 color='white', fontsize=11, pad=10)
    ax.set_xlabel('x', color='#999')
    ax.set_ylabel('y  (scaled mod 10¹⁴, visual)', color='#999')
    ax.tick_params(colors='#777')
    ax.spines[:].set_color('#2a2a4a')
    ax.legend(facecolor='#15153a', labelcolor='white', fontsize=9, framealpha=0.8)
    ax.grid(alpha=0.08, color='white')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("=" * 60)
    print("Demo 1: Seed → SHA-256 → Shamir Shares → Recover")
    print("=" * 60)
    d1 = demo1_shamir_plot(k=2, n=8)
    print(f"  Seed:        {d1['seed'][:32]}...")
    print(f"  SHA-256:     {d1['sha256'][:32]}...")
    print(f"  Secret int:  {d1['secret_int']}")
    print(f"  Shares (k={d1['k']}, n={d1['n']}):")
    for x, y in d1["shares"]:
        print(f"    ({x}, {y})")
    print(f"  Recovered:   {d1['recovered']}")
    print(f"  Verified:    {'✓ PASS' if d1['verified'] else '✗ FAIL'}")

    print()
    print("=" * 60)
    print("Demo 2: Random (x,y) Points → Lagrange → Combined Hash")
    print("=" * 60)
    d2 = demo2_combine_points(n=6)
    print(f"  {d2['n']} random points:")
    for x, y in d2["points"]:
        print(f"    ({x}, {y})")
    print(f"  Combined int: {d2['combined_int']}")
    print(f"  Combined hex: 0x{d2['combined_hex']}")

    print()
    print("Building plot...")
    build_plot(d1, d2, output_path="shamir_demo.png")

"""
============================================================
Demo 1: Seed → SHA-256 → Shamir Shares → Recover
============================================================
  Seed:        9c21b40d2f7858f6a0cc0d65bdd0f4af...
  SHA-256:     ec61365b89c4db3065c407a2a5e5988f...
  Secret int:  72826343057239839408698787436990594475
  Shares (k=2, n=8):
    (1, 131694633145287217187313320160138456005)
    (2, 20421739772865363234240549167402211808)
    (3, 79290029860912741012855081890550073338)
    (4, 138158319948960118791469614613697934868)
    (5, 26885426576538264838396843620961690671)
    (6, 85753716664585642617011376344109552201)
    (7, 144622006752633020395625909067257413731)
    (8, 33349113380211166442553138074521169534)
  Recovered:   72826343057239839408698787436990594475
  Verified:    ✓ PASS

============================================================
Demo 2: Random (x,y) Points → Lagrange → Combined Hash
============================================================
  6 random points:
    (1, 112054735000064433340246958303494305328)
    (2, 57063433438975627894794084655782831575)
    (3, 61567976799816452412058293302800451615)
    (4, 13541281826927087349587510032966600161)
    (5, 10947961206126852512525672061659602648)
    (6, 156444934026489010880362842802506391629)
  Combined int: 73295316376570096884965670249648189579
  Combined hex: 0x37242a1a012736ae4965625583c1708b

Building plot...
Plot saved → shamir_demo.png

"""