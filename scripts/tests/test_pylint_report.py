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

    def test_diagnostics_and_rules_do_not_fill_the_visible_comment(self):
        diagnostic = message()
        diagnostic['message'] = 'Unused import x\nraw code\n<details open>'
        report, _ = REPORT.render(result(9.78, [diagnostic]), result(9.78))
        visible = report.split('<details>')[0]
        self.assertIn('1 new finding.', visible)
        self.assertNotIn('Unused import', visible)
        self.assertNotIn('9.75', visible)
        self.assertIn('- Unused import x', report)
        self.assertNotIn('raw code', report)
        self.assertNotIn('<details open>', report)

    def test_duplicate_code_links_the_actual_files_not_the_reporter(self):
        diagnostic = message()
        diagnostic.update(
            symbol='duplicate-code',
            message='Similar lines in 2 files\n'
            '==PQAnalysis.analysis.msd.msd:[272:279]\n'
            '==PQAnalysis.analysis.vacf.vacf:[247:254]\nraw code',
            path='PQAnalysis/tools/__init__.py',
            line=1
        )
        report, _ = REPORT.render(result(9.78, [diagnostic]), result(9.78),
                                  'https://example.com/blob/abc')
        self.assertIn('/PQAnalysis/analysis/msd/msd.py#L273)', report)
        self.assertIn('/PQAnalysis/analysis/vacf/vacf.py#L248)', report)
        self.assertNotIn('tools/__init__', report)
        self.assertNotIn('raw code', report)

    def test_failure_reason_stays_visible(self):
        report, passed = REPORT.render(result(9.70), result(9.80))
        self.assertFalse(passed)
        visible = report.split('<details>')[0]
        self.assertIn('score must be above 9.75', visible)
        self.assertIn('dropped by more than', visible)

    def test_duplicate_code_sample_changes_are_not_new_findings(self):
        base = message()
        base.update(
            symbol='duplicate-code',
            message=(
                'Similar lines in 2 files\n==PQAnalysis.a:[10:20]\n'
                '==PQAnalysis.b:[30:40]\nexample from a'
            )
        )
        head = dict(
            base,
            message=(
                'Similar lines in 2 files\n==PQAnalysis.b:[30:40]\n'
                '==PQAnalysis.a:[10:20]\nexample from b'
            )
        )
        report, _ = REPORT.render(result(9.78, [head]), result(9.78, [base]))
        self.assertIn('No new findings.', report)

    def test_an_extra_duplicate_code_occurrence_is_still_new(self):
        diagnostic = message()
        diagnostic.update(
            symbol='duplicate-code',
            message=(
                'Similar lines in 2 files\n==PQAnalysis.a:[10:20]\n'
                '==PQAnalysis.b:[30:40]\nexample'
            )
        )
        report, _ = REPORT.render(result(9.78, [diagnostic, diagnostic]),
                                  result(9.78, [diagnostic]))
        self.assertIn('1 new finding.', report)

    def test_duplicate_code_in_other_modules_is_new(self):
        base = message()
        base.update(
            symbol='duplicate-code',
            message=(
                'Similar lines in 2 files\n==PQAnalysis.a:[10:20]\n'
                '==PQAnalysis.b:[30:40]\nexample'
            )
        )
        head = dict(
            base,
            message=base['message'].replace('PQAnalysis.a:', 'PQAnalysis.c:')
        )
        report, _ = REPORT.render(result(9.78, [head]), result(9.78, [base]))
        self.assertIn('1 new finding.', report)

    def test_moved_duplicate_code_block_is_not_new(self):
        base = message()
        base.update(
            symbol='duplicate-code',
            message=(
                'Similar lines in 2 files\n==PQAnalysis.a:[49:63]\n'
                '==PQAnalysis.b:[152:166]\nexample'
            )
        )
        head = dict(
            base,
            message=(
                'Similar lines in 2 files\n==PQAnalysis.a:[58:70]\n'
                '==PQAnalysis.b:[172:184]\nexample'
            )
        )
        report, _ = REPORT.render(result(9.78, [head]), result(9.78, [base]))
        self.assertIn('No new findings.', report)

    def test_artifact_link_and_long_message_summary(self):
        diagnostic = message()
        diagnostic['message'] = 'x' * 500
        report, _ = REPORT.render(result(9.78, [diagnostic]),
                                  artifact_url='https://example.com/artifacts')
        self.assertIn(
            '[Download pylint-report](https://example.com/artifacts)', report
        )
        self.assertIn('x' * 177 + '...', report)
        self.assertNotIn('x' * 181, report)



if __name__ == '__main__':
    unittest.main()
