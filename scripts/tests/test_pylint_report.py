"""Exercise real score gates and shifted-line diagnostic matching."""
import importlib.util
import unittest
from pathlib import Path

SCRIPT = Path(__file__
              ).resolve().parents[2] / '.github/scripts/parse_pylint.py'
SPEC = importlib.util.spec_from_file_location('pylint_report', SCRIPT)
REPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPORT)



def result(score, messages=()):
    return {'statistics': {'score': score}, 'messages': list(messages)}



def message(line=2):
    return dict(
        symbol='unused-import',
        module='PQAnalysis.x',
        obj='',
        message='Unused import x',
        path='PQAnalysis/x.py',
        line=line
    )



class PylintReportTests(unittest.TestCase):

    def test_missing_baseline_is_explicit(self):
        report, passed = REPORT.render(result(9.78))
        self.assertTrue(passed)
        self.assertIn('Base comparison unavailable.', report)

    def test_threshold_is_strict(self):
        self.assertFalse(REPORT.render(result(9.75))[1])

    def test_score_drop_is_checked_against_actual_base(self):
        self.assertFalse(REPORT.render(result(9.78), result(9.84))[1])
        self.assertTrue(REPORT.render(result(9.78), result(9.83))[1])

    def test_shifted_lines_are_not_new_findings(self):
        report, _ = REPORT.render(result(9.78, [message(8)]),
                                  result(9.78, [message(2)]))
        self.assertIn('No new findings.', report)
        self.assertIn('unchanged from base', report)

    def test_added_occurrence_is_new_and_linked(self):
        report, _ = REPORT.render(result(9.78, [message(2), message(8)]),
                                  result(9.78, [message(2)]),
                                  'https://github.com/MolarVerse/PQAnalysis/blob/abc')
        self.assertIn('1 new finding.', report)
        self.assertIn('/PQAnalysis/x.py#L8)', report)
        self.assertNotIn('|', report)



if __name__ == '__main__':
    unittest.main()
