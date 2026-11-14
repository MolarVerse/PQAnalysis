"""Render JSON2 Pylint results and enforce the score gates."""

import argparse
from collections import Counter
import html
import json
from pathlib import Path
from urllib.parse import quote



def finding_key(message):
    """Compare diagnostics independently of shifted line numbers."""
    return tuple(
        message.get(k, "") for k in ("symbol", "module", "obj", "message")
    )



def render(head, base=None, source_url=""):
    """Return a readable report and whether both score gates pass."""
    score = head["statistics"]["score"]
    delta = None if base is None else score - base["statistics"]["score"]
    passed = score > 9.75 and (delta is None or delta >= -0.05 - 1e-9)
    status = "passed" if passed else "failed"
    title = f"**Pylint {status} — {score:.2f}/10"
    if delta is not None:
        title += (
            ", unchanged from base"
            if abs(delta) < 0.005 else f", {delta:+.2f} from base"
        )
    lines = [title + "**", ""]
    lines.append("Minimum required: above 9.75; maximum score drop: 0.05.")
    if base is None:
        lines.append("Base comparison unavailable.")
    messages = head["messages"]
    if base is not None:
        remaining = Counter(finding_key(m) for m in base["messages"])
        new = []
        for message in messages:
            key = finding_key(message)
            if remaining[key]:
                remaining[key] -= 1
            else:
                new.append(message)
        lines.extend(
            [
                "",
                f"{len(new)} new finding{'s' if len(new) != 1 else ''}."
                if new else "No new findings."
            ]
        )
    else:
        new = messages

    def entry(message):
        path = message["path"].replace("\\\\", "/")
        line = message["line"]
        label = f"{path}:{line}"
        if source_url:
            location = f"[{label}]({source_url}/{quote(path)}#L{line})"
        else:
            location = f"`{label}`"
        text = html.escape(message["message"]).replace("`", "'")
        return f"- {location}: {text} (`{message['symbol']}`)"

    if new:
        lines.append("")
        lines.extend(entry(m) for m in new[:5])
        if len(new) > 5:
            lines.append(f"{len(new) - 5} more in the details.")
    lines.extend(
        ["", "<details>", "<summary>Findings and full report</summary>", ""]
    )
    lines.extend(entry(m) for m in messages[:50])
    lines.extend(
        [
            "",
            "All diagnostics are in the pylint-report artifact.",
            "",
            "</details>"
        ]
    )
    return "\n".join(lines) + "\n", passed



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("head")
    parser.add_argument("--base")
    parser.add_argument("--source-url", default="")
    parser.add_argument("--report", default="comment.txt")
    args = parser.parse_args()
    head = json.loads(Path(args.head).read_text(encoding="utf-8"))
    base = (
        json.loads(Path(args.base).read_text(encoding="utf-8"))
        if args.base else None
    )
    report, passed = render(head, base, args.source_url)
    Path(args.report).write_text(report, encoding="utf-8")
    print(report)
    return 0 if passed else 1



if __name__ == "__main__":
    raise SystemExit(main())
