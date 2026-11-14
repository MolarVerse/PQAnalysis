import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "perf_compare.py"
sys.path.insert(0, str(SCRIPT.parent))
SPEC = importlib.util.spec_from_file_location("perf_compare", SCRIPT)
PERF = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PERF)

THRESHOLD = 1.5
MIN_TIME = 1e-4



def verdict(ratios, base=0.01, head=None, threshold=THRESHOLD):
    head = base * ratios[0] if head is None else head
    return PERF.classify(ratios, base, head, threshold, MIN_TIME)



class ClassifyTests(unittest.TestCase):

    def test_equal_timings_are_the_same(self):
        self.assertEqual("same", verdict([1.0, 1.02, 0.98]))

    def test_a_clear_slowdown_in_every_round_is_slower(self):
        self.assertEqual("slower", verdict([1.8, 2.0, 1.7]))

    def test_one_noisy_round_is_not_a_regression(self):
        # median is above the threshold but one round shows no slowdown
        self.assertEqual("same", verdict([2.5, 2.2, 1.05]))

    def test_a_slowdown_below_the_threshold_is_the_same(self):
        self.assertEqual("same", verdict([1.3, 1.4, 1.35]))

    def test_a_clear_speedup_in_every_round_is_faster(self):
        self.assertEqual("faster", verdict([0.4, 0.5, 0.45]))

    def test_benchmarks_that_are_tiny_on_both_sides_are_too_fast(self):
        self.assertEqual(
            "too fast", verdict([5.0, 5.0, 5.0], base=5e-6, head=2.5e-5)
        )

    def test_a_tiny_base_that_became_slow_is_still_judged(self):
        self.assertEqual(
            "slower", verdict([2000.0, 2000.0, 2000.0], base=5e-6, head=1e-2)
        )



class TenPercentLimitTests(unittest.TestCase):

    def verdict(self, ratios):
        return PERF.classify(ratios, 0.01, 0.01 * ratios[0], 1.10, MIN_TIME)

    def test_a_consistent_slowdown_above_ten_percent_is_slower(self):
        self.assertEqual("slower", self.verdict([1.12, 1.15, 1.11]))

    def test_a_slowdown_below_ten_percent_is_the_same(self):
        self.assertEqual("same", self.verdict([1.08, 1.09, 1.07]))

    def test_one_round_below_half_the_limit_prevents_a_regression(self):
        self.assertEqual("same", self.verdict([1.15, 1.14, 1.04]))



class CompareTests(unittest.TestCase):

    def test_rounds_are_paired_and_summarized(self):
        base = [{"a": 0.010}, {"a": 0.012}, {"a": 0.011}]
        head = [{"a": 0.020}, {"a": 0.025}, {"a": 0.022}]

        (result, ) = PERF.compare(base, head, THRESHOLD, MIN_TIME)

        self.assertEqual("a", result["name"])
        self.assertEqual("slower", result["verdict"])
        self.assertAlmostEqual(0.011, result["base"])
        self.assertAlmostEqual(0.022, result["head"])
        self.assertGreater(result["ratio"], 1.9)

    def test_new_and_missing_benchmarks_are_reported(self):
        base = [{"old": 0.01, "both": 0.01}]
        head = [{"new": 0.01, "both": 0.01}]

        verdicts = {
            r["name"]: r["verdict"]
            for r in PERF.compare(base, head, THRESHOLD, MIN_TIME)
        }

        self.assertEqual("missing", verdicts["old"])
        self.assertEqual("new", verdicts["new"])
        self.assertEqual("same", verdicts["both"])

    def test_a_round_where_the_base_failed_is_skipped(self):
        base = [{}, {"a": 0.010}, {"a": 0.010}]
        head = [{"a": 0.05}, {"a": 0.011}, {"a": 0.010}]

        (result, ) = PERF.compare(base, head, THRESHOLD, MIN_TIME)

        self.assertEqual("same", result["verdict"])



class ReportTests(unittest.TestCase):

    def test_load_medians_reads_pytest_benchmark_json(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "result.json"
            path.write_text(
                json.dumps(
                    {
                        "benchmarks": [
                            {
                                "fullname": "x::a[compiled]",
                                "stats": {
                                    "median": 0.5
                                }
                            },
                            {
                                "fullname": "x::a[fallback]",
                                "stats": {
                                    "median": 2.0
                                }
                            },
                        ]
                    }
                ),
                encoding="utf-8",
            )

            self.assertEqual(
                {
                    "x::a[compiled]": 0.5, "x::a[fallback]": 2.0
                },
                PERF.load_medians(path),
            )

    def test_report_lists_regressions_first(self):
        results = PERF.compare(
            [{
                "fast": 0.01, "slow": 0.01
            }] * 3,
            [{
                "fast": 0.01, "slow": 0.03
            }] * 3,
            THRESHOLD,
            MIN_TIME,
        )

        report = PERF.render_markdown(results, THRESHOLD, 3)

        self.assertIn("1 confirmed performance regression", report)
        self.assertLess(report.index("`slow`"), report.index("`fast`"))

    def test_report_without_regressions_says_so(self):
        results = PERF.compare(
            [{
                "a": 0.01
            }] * 3, [{
                "a": 0.0105
            }] * 3, THRESHOLD, MIN_TIME
        )

        report = PERF.render_markdown(results, THRESHOLD, 3)

        self.assertIn("No confirmed performance regressions", report)

    def test_changes_are_reported_as_signed_percentages(self):
        results = PERF.compare(
            [{
                "up": 0.010, "down": 0.010, "flat": 0.010
            }] * 3,
            [{
                "up": 0.0112, "down": 0.0095, "flat": 0.010
            }] * 3,
            1.10,
            MIN_TIME,
        )

        report = PERF.render_markdown(results, 1.10, 3)

        self.assertIn("(+12.0%; slower)", report)
        self.assertIn("(-5.0%; same)", report)
        self.assertIn("(+0.0%; same)", report)
        self.assertIn("median paired-round change above +10.0%", report)

    def test_unchanged_results_are_collapsed_and_no_tables(self):
        results = PERF.compare(
            [{
                "a": 0.01
            }] * 3, [{
                "a": 0.01
            }] * 3, THRESHOLD, MIN_TIME
        )
        report = PERF.render_markdown(results, THRESHOLD, 3)
        self.assertLess(report.index("<details>"), report.index("`a`"))
        self.assertNotIn("|", report)
        self.assertIn("paired head/base ratios", report)

    def test_empty_report_does_not_claim_success(self):
        report = PERF.render_markdown([], THRESHOLD, 3)
        self.assertIn("comparison unavailable", report)
        self.assertNotIn("No confirmed", report)

    def test_format_change(self):
        self.assertEqual("+3.2%", PERF.format_change(1.032))
        self.assertEqual("-10.0%", PERF.format_change(0.9))
        self.assertEqual("+0.0%", PERF.format_change(1.0))
        self.assertEqual("+0.0%", PERF.format_change(0.9996))

    def test_format_time_picks_a_unit(self):
        self.assertEqual("1.5 s", PERF.format_time(1.5))
        self.assertEqual("2.5 ms", PERF.format_time(2.5e-3))
        self.assertEqual("7 us", PERF.format_time(7e-6))



if __name__ == "__main__":
    unittest.main()
