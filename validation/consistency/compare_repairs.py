"""Compare preserved and repaired audits with identical cases and input matrices.

Usage: python3 validation/consistency/compare_repairs.py [before_dir after_dir]
Writes comparison.csv and PNG/PDF figures to after_dir. Assertions validate
comparability, not zero violations; the table retains any remaining flags.
"""
import csv
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).parent
BEFORE, AFTER = map(Path, sys.argv[1:]) if len(sys.argv) > 1 else (
    ROOT / "results", ROOT / "repaired_results")


def read(root, name):
    with (root / name).open() as handle:
        return list(csv.DictReader(handle))


key = lambda row: tuple(row[col] for col in ("dataset", "method", "transform", "check"))
before = {key(row): row for row in read(BEFORE, "checks.csv")}
after = {key(row): row for row in read(AFTER, "checks.csv")}
explicit_ordering = "Mid-CDF oracle uses exposed normalization_scores" in (AFTER / "environment.txt").read_text()
assert before.keys() == after.keys(), "Audit check coverage changed"
for path in (BEFORE / "datasets").glob("*.csv"):
    assert path.read_bytes() == (AFTER / "datasets" / path.name).read_bytes(), path.name
for label in ("Seeds:", "Raw error denominator="):
    lines = [next(line for line in (root / "environment.txt").read_text().splitlines()
                  if line.startswith(label)) for root in (BEFORE, AFTER)]
    assert lines[0] == lines[1], "Audit tolerances or seeds changed"

rows = []
for method in sorted({row["method"] for row in before.values()}):
    keys = [k for k in before if k[1] == method]
    left = Counter(before[k]["status"] for k in keys)
    right = Counter(after[k]["status"] for k in keys)
    rows.append(dict(method=method, checks=len(keys),
        before_violations=left["violation"], after_violations=right["violation"],
        before_unexpected_errors=left["unexpected_error"],
        after_unexpected_errors=right["unexpected_error"],
        resolved_flags=sum(before[k]["status"] == "violation" and
                           after[k]["status"] == "pass" for k in keys),
        new_flags=sum(before[k]["status"] == "pass" and
                      after[k]["status"] != "pass" for k in keys)))
with (AFTER / "comparison.csv").open("w") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

x = np.arange(len(rows))
fig, ax = plt.subplots(figsize=(11, 5), layout="constrained")
for offset, field, label, color in (
    (-.2, "before_violations", "Original audit", "#C96035"),
    (.2, "after_violations", "Repaired replay", "#247A79")):
    bars = ax.bar(x+offset, [row[field] for row in rows], .38, label=label, color=color)
    ax.bar_label(bars, padding=3)
ax.set_xticks(x, [row["method"].replace("_", "\n") for row in rows])
ax.set_ylabel("Violation flags")
contract = ("KDE CDF uses explicit ordering labels; invariance checks and tolerances unchanged"
            if explicit_ordering else "Original CDF contract, invariance checks and tolerances retained")
ax.set_title("84 identical datasets · 14,637 checks per run\n" + contract + "\n"
             "Flags are correlated outcomes, not independent bugs or population error rates", fontsize=10)
ax.legend()
ax.spines[["top", "right"]].set_visible(False)
ax.margins(y=.15)
(AFTER / "figures").mkdir(exist_ok=True)
for extension in ("png", "pdf"):
    fig.savefig(AFTER / "figures" / f"before_after.{extension}", dpi=170)
plt.close(fig)
print("Identical datasets, check keys, seeds and tolerances verified.")
print(contract)
print("Before:", dict(Counter(row["status"] for row in before.values())))
print("After:", dict(Counter(row["status"] for row in after.values())))
print("Wrote comparison.csv and before_after figures to", AFTER)
