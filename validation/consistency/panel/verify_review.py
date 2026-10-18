"""Check the review's evidence links and reported counts, not estimator correctness.

Run from any directory: python3 validation/consistency/panel/verify_review.py
"""
import csv
import hashlib
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PANEL = ROOT / "validation/consistency/panel"
REPORTS = ROOT / "docs/research/panel_followup"


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


roles = ("statistics", "mathematics", "hpc", "julia", "computational_social_science",
         "sociology", "communication")
for role in roles:
    assert (REPORTS / f"{role}.md").is_file(), role
for path in REPORTS.glob("*.md"):
    for link in re.findall(r"\]\(([^)]+)\)", path.read_text()):
        if not link.startswith(("http://", "https://", "#")):
            assert (path.parent / link.split("#")[0]).exists(), (path, link)

source = "\n".join(path.read_text() for path in sorted((ROOT / "src").iterdir()))
digest = hashlib.sha256(source.encode()).hexdigest()
assert digest in (ROOT / "validation/consistency/results/environment.txt").read_text()
flags = [r for r in rows(ROOT / "validation/consistency/results/checks.csv")
         if r["status"] == "violation"]
assert len(flags) == 630
assert sum(float(r["raw_max_relative"]) <= 1e-10 for r in flags) == 525
assert sum(float(r["summary_max_delta"]) > 1e-10 for r in flags) == 603

grouping = rows(PANEL / "statistics/feature_exchange_attribution.csv")
assert len(grouping) == 91
assert all(float(r["baseline_score"]) > 1e-10 for r in grouping)
assert all(float(r["prototype_raw"]) == float(r["prototype_score"]) == 0 for r in grouping)
suite = (PANEL / "integration/existing_suite.txt").read_text()
counts = re.findall(r"\|\s+(\d+)\s+(\d+)\s+[\d.]+s", suite)
assert sum(int(passed) for passed, total in counts) == 911
assert all(passed == total for passed, total in counts)
assert "Error" not in suite and "Fail" not in suite

message = ("Seven specialist reports and their local links verified.\n"
           "Source fingerprint matches original audit: " + digest + "\n"
           "630 flags: 525 raw-within-tolerance, 105 raw-above-tolerance; 603 affect a summary.\n"
           "91 feature-order interventions: original violations retained; prototype differences zero.\n"
           "Existing suite: 911/911 passed in retained clean log.\n"
           "This verifies review artifacts and counts, not resolution of package defects.\n")
(PANEL / "integration/review_integrity.txt").write_text(message)
print(message, end="")
