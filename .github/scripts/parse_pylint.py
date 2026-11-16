"""Render JSON2 Pylint results and enforce the score gates."""

import argparse
from collections import Counter
import html
import json
from pathlib import Path
import re
from urllib.parse import quote



def finding_key(message):
    """Match diagnostics without relying on excerpts or line numbers."""
    if message.get("symbol") == "duplicate-code":
        modules = re.findall(
            r"^==([\w.]+):\[\d+:\d+\]", message["message"], re.MULTILINE
        )
        if modules:
            # moved or shrunk blocks keep their key; a new pair or count is new
            return ("duplicate-code", tuple(sorted(modules)))
    return tuple(
        message.get(k, "") for k in ("symbol", "module", "obj", "message")
    )



def render(head, base=None, source_url="", artifact_url=""):
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
    lines = [title + "**"]
    if score <= 9.75:
        lines.extend(["", "The score must be above 9.75."])
    if delta is not None and delta < -0.05 - 1e-9:
        lines.extend(["", "The score dropped by more than the allowed 0.05."])
    if base is None:
        lines.extend(["", "Base comparison unavailable."])
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
        lines.extend(["", f"{len(messages)} findings."])

    def location(path, line):
        path = path.replace("\\", "/")
        label = f"{path}:{line}"
        if source_url:
            return f"[{label}]({source_url}/{quote(path)}#L{line})"
        return f"`{label}`"

    def entry(message):
        text = (message["message"].splitlines() or [message["symbol"]])[0]
        text = text[:177] + "..." if len(text) > 180 else text
        text = html.escape(text).replace("`", "'")
        text = re.sub(r"([\\*_[\]])", r"\\\1", text)
        duplicates = re.findall(
            r"^==([\w.]+):\[(\d+):\d+\]", message["message"], re.MULTILINE
        ) if message["symbol"] == "duplicate-code" else []
        places = [
            location(module.replace(".", "/") + ".py", int(start) + 1)
            for module, start in duplicates
        ] or [location(message["path"], message["line"])]
        return f"- {text} — {', '.join(places)}"

    if new:
        label = "New findings" if base is not None else "Findings"
        lines.extend(
            ["", "<details>", f"<summary>{label} ({len(new)})</summary>", ""]
        )
        lines.extend(entry(m) for m in new[:50])
        if len(new) > 50:
            lines.append(f"\n{len(new) - 50} more in the downloadable report.")
        lines.extend(["", "</details>"])
    lines.extend(
        [
            "",
            "<details>",
            "<summary>Full diagnostics and gate rules</summary>",
            ""
        ]
    )
    diagnostics = (
        f"[Download pylint-report]({artifact_url}) for all diagnostics."
        if artifact_url else
        "All diagnostics are in the pylint-report artifact on the workflow run."
    )
    lines.extend(
        [
            "The score must be above 9.75 and must not drop by more than 0.05.",
            "",
            diagnostics,
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
    parser.add_argument("--artifact-url", default="")
    parser.add_argument("--report", default="comment.txt")
    args = parser.parse_args()
    head = json.loads(Path(args.head).read_text(encoding="utf-8"))
    base = (
        json.loads(Path(args.base).read_text(encoding="utf-8"))
        if args.base else None
    )
    report, passed = render(head, base, args.source_url, args.artifact_url)
    Path(args.report).write_text(report, encoding="utf-8")
    print(report)
    return 0 if passed else 1



if __name__ == "__main__":
    raise SystemExit(main())
