#!/usr/bin/env python3

"""Compare the benchmark timings of a base and a head revision.

Both revisions are installed in their own environment on the same machine.
The benchmark suite of the *head* revision is run against both installations,
alternating base and head for several rounds, so that slow drift of the
machine hits both sides alike. A benchmark is reported as slower only if the
median change of the per-round head/base ratios exceeds the threshold (+10 %
by default) and every single round shows at least half of it; benchmarks that are too fast to be timed
reliably are reported but never gate.

Usage:
    perf_compare.py --base-python BASE/bin/python --base-src BASE_SRC \
        --head-python HEAD/bin/python --head-src HEAD_SRC [--threshold-percent 10]
"""

import argparse
import json
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

DEFAULT_THRESHOLD_PERCENT = 10.0
DEFAULT_ROUNDS = 3
DEFAULT_MIN_TIME = 1e-4

PYTEST_ARGS = [
    "--benchmark-only",
    "--benchmark-min-rounds=5",
    "--benchmark-max-time=0.3",
    "--benchmark-warmup=off",
    "-q",
    "-p",
    "no:cacheprovider",
]



def load_medians(path):
    """Return a mapping of benchmark name to median seconds from a JSON."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return {
        bench["fullname"]: bench["stats"]["median"]
        for bench in data.get("benchmarks", [])
    }



def classify(ratios, base_median, head_median, threshold, min_time):
    """Return the verdict for the per-round head/base ratios of a benchmark."""
    if max(base_median, head_median) < min_time:
        return "too fast"

    middle = statistics.median(ratios)
    slower_round = 1.0 + (threshold - 1.0) / 2.0

    if middle > threshold and min(ratios) > slower_round:
        return "slower"
    if middle < 1.0 / threshold and max(ratios) < 1.0 / slower_round:
        return "faster"
    return "same"



def compare(base_rounds, head_rounds, threshold, min_time):
    """Compare lists of per-round median mappings of base and head.

    Returns a list of result dictionaries sorted by name.
    """
    names = set()
    for rounds in (base_rounds, head_rounds):
        for medians in rounds:
            names.update(medians)

    results = []
    for name in sorted(names):
        base = [m[name] for m in base_rounds if name in m]
        head = [m[name] for m in head_rounds if name in m]

        if not base:
            results.append({"name": name, "verdict": "new"})
            continue
        if not head:
            results.append({"name": name, "verdict": "missing"})
            continue

        ratios = [
            h / b for b_round, h_round in zip(base_rounds, head_rounds)
            if name in b_round and name in h_round
            for b, h in [(b_round[name], h_round[name])]
        ]
        base_median = statistics.median(base)
        head_median = statistics.median(head)
        results.append(
            {
                "name":
                name,
                "base":
                base_median,
                "head":
                head_median,
                "ratio":
                statistics.median(ratios),
                "verdict":
                classify(
                    ratios, base_median, head_median, threshold, min_time
                ),
            }
        )

    return results



def format_change(ratio):
    """Format a head/base ratio as a signed percentage change."""
    change = round((ratio - 1.0) * 100.0, 1) + 0.0  # no "-0.0"
    return f"{change:+.1f}%"



def format_time(seconds):
    """Format seconds with a readable unit."""
    for unit, factor in (("s", 1.0), ("ms", 1e-3), ("us", 1e-6)):
        if seconds >= factor:
            return f"{seconds / factor:.3g} {unit}"
    return f"{seconds / 1e-9:.3g} ns"



def render_markdown(results, threshold, rounds, header=None):
    """Render a short result with expandable timings and comparison rules."""
    slower = [r for r in results if r["verdict"] == "slower"]
    faster = [r for r in results if r["verdict"] == "faster"]
    compared = [r for r in results if "ratio" in r]
    if not compared:
        title = "Performance comparison unavailable"
    elif slower:
        noun = "regression" if len(slower) == 1 else "regressions"
        title = f"{len(slower)} confirmed performance {noun}"
    else:
        title = "No confirmed performance regressions"
    lines = [f"**{title}**", ""]
    if header:
        lines.extend([header, ""])
    lines.extend(
        [f"{len(compared)} benchmarks compared across {rounds} rounds.", ""]
    )

    def entry(result):
        name = result["name"]
        if "ratio" not in result:
            return f"- `{name}`: {result['verdict']}"
        return (
            f"- `{name}`: {format_time(result['base'])} → "
            f"{format_time(result['head'])} "
            f"({format_change(result['ratio'])}; {result['verdict']})"
        )

    changed = slower + faster
    lines.extend(entry(r) for r in changed[:5])
    if len(changed) > 5:
        lines.append(f"{len(changed) - 5} more changes in the full results.")
    if changed:
        lines.append("")
    lines.extend(
        [
            "<details>",
            "<summary>Full results and comparison rules</summary>",
            "",
            f"A regression needs a median paired-round change above "
            f"{format_change(threshold)} and every paired round above "
            f"{format_change(1.0 + (threshold - 1.0) / 2.0)}.",
            "",
            "Base and head are separate medians. Change is the median of the "
            "paired head/base ratios, so it need not equal the ratio of the "
            "displayed times. Tiny benchmarks are reported but not gated.",
            "",
        ]
    )
    rank = {"slower": 0, "faster": 1, "same": 2, "too fast": 3}
    lines.extend(
        entry(r) for r in
        sorted(results, key=lambda r: (rank.get(r["verdict"], 4), r["name"]))
    )
    lines.extend(["", "</details>"])
    return "\n".join(lines) + "\n"



def package_location(python, cwd):
    """Return the directory that ``python`` imports PQAnalysis from."""
    output = subprocess.run(
        [python, "-c", "import PQAnalysis; print(PQAnalysis.__file__)"],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    return Path(output).resolve().parent.parent



def run_benchmarks(python, benchmarks, config, output, cwd):
    """Run the benchmark suite with ``python``; returns the exit code."""
    result = subprocess.run(
        [
            python,
            "-m",
            "pytest",
            "-c",
            str(config),
            str(benchmarks),
            f"--benchmark-json={output}",
            *PYTEST_ARGS,
        ],
        cwd=cwd,
        check=False,
    )
    return result.returncode



def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--base-python", required=True)
    parser.add_argument("--base-src", required=True)
    parser.add_argument("--head-python", required=True)
    parser.add_argument("--head-src", required=True)
    parser.add_argument(
        "--threshold-percent",
        type=float,
        default=DEFAULT_THRESHOLD_PERCENT,
        help="slowdown in percent above which a benchmark counts as slower",
    )
    parser.add_argument("--rounds", type=int, default=DEFAULT_ROUNDS)
    parser.add_argument("--min-time", type=float, default=DEFAULT_MIN_TIME)
    parser.add_argument("--report", default="perf_report.md")
    parser.add_argument(
        "--no-fail",
        action="store_true",
        help="always exit 0 (the report still lists regressions)",
    )
    args = parser.parse_args()

    threshold = 1.0 + args.threshold_percent / 100.0

    base_src = Path(args.base_src).resolve()
    head_src = Path(args.head_src).resolve()
    benchmarks = head_src / "benchmarks"
    config = benchmarks / "pytest.ini"

    # Run from a neutral directory: the current directory would otherwise
    # shadow the installed package of the other revision.
    with tempfile.TemporaryDirectory() as work:
        for name, python, src in (
            ("base", args.base_python, base_src),
            ("head", args.head_python, head_src),
        ):
            found = package_location(python, work)
            if found != src:
                sys.exit(f"{name} imports PQAnalysis from {found}, not {src}")

        base_rounds, head_rounds = [], []
        head_failed = False
        for index in range(args.rounds):
            order = [("base", args.base_python), ("head", args.head_python)]
            if index % 2:
                order.reverse()

            for side, python in order:
                output = Path(work) / f"{side}-{index}.json"
                code = run_benchmarks(python, benchmarks, config, output, work)
                if side == "head" and code != 0:
                    head_failed = True
                medians = load_medians(output) if output.is_file() else {}
                (base_rounds
                 if side == "base" else head_rounds).append(medians)

    results = compare(base_rounds, head_rounds, threshold, args.min_time)
    report = render_markdown(results, threshold, args.rounds)
    Path(args.report).write_text(report, encoding="utf-8")
    print(report)

    slower = [r for r in results if r["verdict"] == "slower"]
    if head_failed:
        print("benchmarks of the head revision failed", file=sys.stderr)
    if not any("ratio" in r for r in results):
        print("no benchmark could be compared", file=sys.stderr)
        return 1
    if head_failed or (slower and not args.no_fail):
        return 1
    return 0



if __name__ == "__main__":
    sys.exit(main())
