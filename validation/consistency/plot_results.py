"""Plot saved audit outputs. Requires Python, numpy and matplotlib; no refitting.

    python3 validation/consistency/plot_results.py [results_directory]
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / "results"
FIGURES = ROOT / "figures"
FIGURES.mkdir(exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False})


def read(name):
    with (ROOT / name).open() as handle:
        return list(csv.DictReader(handle))


def save(fig, name):
    fig.savefig(FIGURES / f"{name}.png", dpi=170, bbox_inches="tight")
    fig.savefig(FIGURES / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


catalog = [row for row in read("datasets.csv") if row["seed"] == "101"]
fig, axes = plt.subplots(7, 4, figsize=(15, 20), layout="constrained")
for ax, item in zip(axes.flat, catalog):
    points = read(f"datasets/{item['dataset']}.csv")
    x, y = (np.array([float(row[col]) for row in points]) for col in ("f1", "f2"))
    scenario = item["scenario"]
    suffix = ""
    if scenario in ("tiny_units", "huge_units", "mixed_units"):
        x, y = x / np.max(np.abs(x)), y / np.max(np.abs(y))
        suffix = " (axes rescaled)"
    elif scenario == "large_offset":
        x, y = x - 1e14, y - 1e14
        suffix = " (offset removed)"
    color = ["#D97920" if row["planted"].lower() == "true" else "#256A99" for row in points]
    ax.scatter(x, y, c=color, s=19, alpha=.7, linewidths=0)
    ax.set_title(scenario.replace("_", " ") + f"\nn={item['n']}, p={item['features']}" + suffix, fontsize=10)
    ax.set_xlabel("Feature 1"); ax.set_ylabel("Feature 2")
fig.suptitle("28 simulated data families · seed 101\nFirst two features only; duplicates overlap. Orange marks planted diagnostic rows, not a correctness oracle.", fontsize=16)
save(fig, "dataset_atlas")

checks = read("checks.csv")
methods = ["lof", "knn", "kde", "depth", "gaussian", "gaussian_regularized", "beta", "graph"]
names = ["LOF", "kNN", "KDE", "Depth", "Gaussian copula", "Regularized Gaussian", "Beta copula", "Graph"]
properties = ["coordinate_reflection", "feature_permutation", "row_permutation",
              "append_zero_weight_uniform", "append_zero_weight_nonuniform",
              "weight_rescale_1.0e-150", "weight_rescale_1.0e150",
              "common_units_2^-200", "common_units_2^200"]
labels = ["Reflect\ncoordinates", "Reverse\nfeatures", "Permute\nrows", "Add zero weights\nuniform reference",
          "Add zero weights\nweighted reference", "Weights\n× 10⁻¹⁵⁰", "Weights\n× 10¹⁵⁰",
          "Units\n× 2⁻²⁰⁰", "Units\n× 2²⁰⁰"]
groups = defaultdict(list)
for row in checks:
    groups[row["method"], row["check"]].append(row)
rate = np.full((len(methods), len(properties)), np.nan)
annotation = {}
for i, method in enumerate(methods):
    for j, prop in enumerate(properties):
        rows = groups[method, prop]
        if rows:
            failures = sum(row["status"] == "violation" for row in rows)
            rate[i, j] = 100 * failures / len(rows)
            annotation[i, j] = f"{failures}/{len(rows)}"
fig, ax = plt.subplots(figsize=(16, 6), layout="constrained")
cmap = plt.get_cmap("YlOrRd").copy(); cmap.set_bad("#eeeeee")
im = ax.imshow(rate, cmap=cmap, vmin=0, vmax=100, aspect="auto")
ax.set_xticks(range(len(properties)), labels, fontsize=9)
ax.set_yticks(range(len(methods)), names)
for i in range(len(methods)):
    for j in range(len(properties)):
        ax.text(j, i, annotation.get((i, j), "—"), ha="center", va="center",
                color="white" if rate[i, j] > 65 else "#222222", fontsize=10)
fig.colorbar(im, ax=ax, label="Flagged comparisons / attempted comparisons (%)", shrink=.8)
ax.set_title("Consistency flags by method and transformation\nCells show flagged/attempted; pooled transforms and seeds. These are dependent stress tests, not population error rates.", pad=18)
save(fig, "consistency_flags")

examples = read("minimal_examples.csv")
fig, axes = plt.subplots(1, 3, figsize=(15, 4.3), layout="constrained")
for ax, example, title in zip(axes[:2], ["triangle_depth_rank", "weighted_constant_kde"],
                             ["Depth: three noncollinear points", "KDE: three identical weighted points"]):
    rows = [row for row in examples if row["example"] == example]
    ax.bar([int(row["entity_id"]) for row in rows], [float(row["score_observed"]) for row in rows],
           color="#256A99", width=.6, label="Observed percentile")
    ax.axhline(.5, color="#D97920", linewidth=2, linestyle="--", label="Mathematical expectation")
    ax.set(title=title, xlabel="Entity", ylabel="Typicality percentile", ylim=(0, 1), xticks=[1, 2, 3])
    ax.legend(fontsize=8, loc="upper right")
original = [row for row in examples if row["example"] == "lof_original"]
reflected = [row for row in examples if row["example"] == "lof_reflected"]
x = np.arange(1, 6)
ax = axes[2]
for shift, values, color, label in [(-.25, [float(r["raw_observed"]) for r in original], "#256A99", "Original"),
                                  (0, [float(r["raw_observed"]) for r in reflected], "#D97920", "Reflected"),
                                  (.25, [float(r["raw_expected"]) for r in original], "#288071", "Independent oracle")]:
    ax.bar(x + shift, values, width=.24, color=color, label=label)
ax.set(title="LOF: reflection and independent oracle", xlabel="Entity", ylabel="Raw typicality (inverse LOF)", xticks=x)
ax.legend(fontsize=8)
fig.suptitle("Small fixtures with independently established expectations", fontsize=15)
save(fig, "minimal_counterexamples")
print(f"Wrote PNG and PDF figures to {FIGURES}")
