CI report comments
==================

Pylint reports the score, its change from the PR base, and new diagnostics.
Only the result and finding count are visible by default. Expand New findings
for short, linked summaries; duplicate-code findings link to the compared files.
Raw multiline diagnostics stay in the downloadable artifact. Failed score gates
are explained in the visible comment; passing gate rules stay collapsed.
Duplicate-code comparisons use the reported file ranges, not Pylint's choice of
example excerpt, so a different excerpt does not become a new finding.
Both revisions use the same dependency versions. The score must exceed 9.75
and must not fall by more than 0.05. Push and manual runs have no base and say
so explicitly. Full JSON diagnostics are saved as the pylint-report artifact.

Performance reports lead with confirmed regressions and improvements.
Visible changes use readable benchmark names and state how much slower or faster
they are; exact benchmark identifiers remain in the full results.
Results and comparison rules are expandable; unchanged timings stay there.
The displayed change is the median of paired-round ratios, not necessarily
the ratio of the separate base and head medians.

The Publish performance report workflow runs separately after Performance.
It reads only the report artifact and never checks out or executes PR code.
It ignores closed, stale and external PR runs. The app token is restricted
to pull-request posting and never enters the benchmark job.

Install pq-perf-bot for this repository with pull-request write access and
configure repository or selected-repository organization secrets named
PQ_PERF_BOT_APP_ID and PQ_PERF_BOT_PRIVATE_KEY. Without them, publishing logs
a warning; reports remain in the benchmark summary and artifact. There is
no fallback to posting as github-actions. The publisher must be on the
repository default branch before workflow_run can trigger it.

The bot updates its own marked comment. Historical github-actions comments
are left untouched; GitHub cannot change a comment's author.
