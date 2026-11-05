# /// script
# dependencies = ["matplotlib"]
# ///
"""Plot the reference FER curves in sim/logs against the curves of Shen et al., JSAC 2025.

Run from anywhere: uv run scripts/plot_reproduction.py (uv installs matplotlib on the fly),
or python scripts/plot_reproduction.py with matplotlib installed.
Writes media/reproduction_results.png and prints, per simulated point, FER(sim)/FER(paper)
with the number of frame errors in parentheses.

The paper curves were extracted once from the vector drawings of Figs. 4, 6, 8 and 9 (the solid
red PCM-II line, or the PCM-I markers of Fig. 4) and are stored below. Polar curves are against
CRC-inclusive Eb/N0, as in the paper.
"""
import csv, math
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parent.parent
C_PAPER, C_SIM, C_ALT = "#c8442a", "#1f5fa8", "#2a9d5c"

PAPER = {  # (Eb/N0 dB, FER)
    "RM(32,16) PCM-I": [(2, 0.3788), (2.5, 0.312), (3, 0.2008), (3.5, 0.1521), (4, 0.08358), (4.5, 0.04561), (5, 0.02192), (5.5, 0.01008), (6, 0.004065), (6.5, 0.00112), (7, 0.0003345), (7.5, 7.199e-05), (8, 1.7e-05)],
    "RM(32,16) PCM-II": [(2, 0.1209), (2.5, 0.07236), (3, 0.03355), (3.5, 0.01384), (4, 0.005878), (4.5, 0.001841), (5, 0.000445), (5.5, 0.000108), (6, 2.02e-05), (6.5, 2.619e-06), (7, 2.709e-07)],
    "BCH(63,45) PCM-II": [(2, 0.2272), (2.5, 0.1325), (3, 0.05222), (3.5, 0.01656), (4, 0.005114), (4.5, 0.001171), (5, 0.000202), (5.5, 2.499e-05), (6, 3.098e-06)],
    "eBCH(64,24) PCM-II": [(2, 0.05584), (2.5, 0.02338), (3, 0.008445), (3.5, 0.00245), (4, 0.000524), (4.5, 8.599e-05), (5, 1.14e-05), (5.5, 1.02e-06)],
    "RM(128,64)": [(2, 0.4568), (2.5, 0.2704), (3, 0.1224), (3.5, 0.03886), (4, 0.008695), (4.5, 0.001304), (5, 0.000127), (5.5, 6.327e-06), (6, 2.389e-07)],
    "Polar(128,20) + CRC24": [(0, 0.03351), (0.5, 0.01084), (1, 0.002982), (1.5, 0.000571), (2, 7.898e-05), (2.5, 6.209e-06), (3, 4.799e-07)],
    "Polar(128,40) + CRC24": [(0, 0.4039), (0.5, 0.2271), (1, 0.1075), (1.5, 0.04444), (2, 0.01195), (2.5, 0.002761), (3, 0.000511), (3.5, 5.299e-05), (4, 5.087e-06)],
    "Polar(128,72) + CRC24": [(1, 0.5695), (1.5, 0.3463), (2, 0.1641), (2.5, 0.05468), (3, 0.01379), (3.5, 0.002088), (4, 0.0002177)],
}

def sim(code, shift=0.0):
    with open(ROOT / "sim" / "logs" / f"{code}.reference.csv") as f:
        rows = [r for r in csv.DictReader(f) if int(r["Frame Errors"]) > 0]
    return [(float(r["Eb/N0"]) - shift, float(r["FER"]), int(r["Frame Errors"])) for r in rows]

sh = lambda k: 10*math.log10((k+24)/k)   # payload-normalized → CRC-inclusive Eb/N0

# title, iterations, paper note, [(code, note, Eb/N0 shift)], x range
PANELS = [
    ("RM(32,16) PCM-I", 6, "Fig. 4, 32×48", [("rm.32.16.pcm1", "28×44", 0)], (2, 7)),
    ("RM(32,16) PCM-II", 6, "Fig. 4, 48×38", [("rm.32.16.pcm2", "48×38", 0)], (2, 7)),
    ("BCH(63,45) PCM-II", 20, "Fig. 9, 63×72", [("bch.63.45.pcm2", "63×72", 0)], (2, 6.5)),
    ("eBCH(64,24) PCM-II", 20, "Fig. 6, 272×81",
     [("ebch.64.24.pcm2", "272×81", 0), ("ebch.64.24.pcm2.eta52", "336×85, η̄=52", 0)], (2, 6)),
    ("RM(128,64)", 20, "Fig. 9, PCM-II 141×141",
     [("rm.128.64.pcm2", "PCM-II 147×161", 0), ("rm.128.64.geom.pcm1", "geom. PCM-I 144×137", 0)], (2, 6)),
] + [(f"Polar(128,{k}) + CRC24", 200, "Fig. 8, " + s, [(f"polar.128.{k}.crc24c.pcm2", t, sh(k))],
      (0 if k < 72 else 1, 3 if k == 20 else 4))
     for k, s, t in [(20, "323×215", "323×215"), (40, "329×283", "329×283"), (72, "409×420", "472×458")]]

def draw(ax, title, iters, pnote, sims, xr):
    paper = PAPER[title]
    ys = [y for x, y in paper if xr[0] <= x <= xr[1]]
    ax.semilogy(*zip(*paper), "-s", color=C_PAPER, mfc="white", mew=1.4, lw=1.8, ms=5.5, zorder=3)
    for (code, note, s), col, mk in zip(sims, [C_SIM, C_ALT], ["o", "^"]):
        x, y, e = zip(*sim(code, s))
        lo = [f*max(1-1.96/math.sqrt(n), 0.05) for f, n in zip(y, e)]   # ponytail: normal approx, fine at >=50 errors
        hi = [f*(1+1.96/math.sqrt(n)) for f, n in zip(y, e)]
        ax.fill_between(x, lo, hi, color=col, alpha=0.18, lw=0, zorder=1)
        ax.semilogy(x, y, "--" + mk, color=col, lw=1.5, ms=4.5, zorder=4)
        ys += list(y)
    ax.set_xlim(*xr)
    ax.set_ylim(10**math.floor(math.log10(min(ys)) - 0.15), 1)
    ax.set_xticks(range(math.ceil(xr[0]), math.floor(xr[1]) + 1))
    ax.grid(True, which="major", color="#d0d0d0", lw=0.8)
    ax.grid(True, which="minor", axis="y", color="#ececec", lw=0.5)
    ax.set_title(title, fontsize=11.5, fontweight="bold", loc="left")
    ax.set_title(f"NMS, {iters} it.", fontsize=9, color="#666", loc="right")
    notes = [(f"paper: {pnote}", C_PAPER)] + [(f"sim: {n}", c) for (_, n, _), c in zip(sims, [C_SIM, C_ALT])]
    for i, (t, c) in enumerate(reversed(notes)):
        ax.text(0.03, 0.04 + 0.075*i, t, transform=ax.transAxes, fontsize=8.5, color=c)
    ax.set_xlabel("$E_b/N_0$ [dB]" + ("  (CRC-inclusive)" if "Polar" in title else ""), fontsize=9.5)
    ax.tick_params(labelsize=9)

if __name__ == "__main__":
    plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#555"})
    fig, axs = plt.subplots(2, 4, figsize=(19, 9.6))
    for ax, p in zip(axs.flat, PANELS): draw(ax, *p)
    for ax in axs[:, 0]: ax.set_ylabel("FER", fontsize=10)
    handles = [Line2D([], [], color=C_PAPER, marker="s", mfc="white", lw=1.8, label="Paper (extracted from vector figure)"),
               Line2D([], [], color=C_SIM, marker="o", ls="--", label="Simulation"),
               Line2D([], [], color=C_ALT, marker="^", ls="--", label="Simulation, alternative PCM"),
               Patch(color="#888", alpha=0.3, label="95% confidence band")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=10.5, bbox_to_anchor=(0.5, 0))
    fig.suptitle("Frame error rate: optimized PCMs vs. Shen et al., JSAC 2025", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(ROOT / "media" / "reproduction_results.png", dpi=150)
    for title, _, _, sims, _ in PANELS:
        P = dict(PAPER[title])
        for code, _, s in sims:
            r = [(x, y/P[round(x*2)/2], e) for x, y, e in sim(code, s) if abs(x*2-round(x*2)) < 1e-6 and round(x*2)/2 in P]
            print(f"{code:26s}", " ".join(f"{x:g}:{q:.2f}({e})" for x, q, e in r))
