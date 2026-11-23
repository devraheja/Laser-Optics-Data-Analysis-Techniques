import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

DATA = ROOT / "data"
FIG = ROOT / "figures"
RESULTS = ROOT / "results"
FIG.mkdir(exist_ok=True)
RESULTS.mkdir(exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 130,
    "font.size": 10,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "legend.fontsize": 9,
})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIG / name)
    plt.close(fig)


def write_results(name, lines):
    (RESULTS / name).write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
