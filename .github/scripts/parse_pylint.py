import sys


def main(file):
    with open(file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    summary = [line for line in lines if line.startswith(
        "Your code has been rated at")][0]

    report_start_index = [i for i, line in enumerate(
        lines) if line.startswith("Raw metrics")][0]
    report_end_index = [i for i, line in enumerate(
        lines) if line.startswith("Your code has been rated at")][0]

    report = lines[report_start_index:report_end_index-2]

    new_report = []
    for line in report:
        if len(line.strip()) == 0:
            new_report.append(line)
        elif len(line.strip().replace("-", "")) == 0:
            line = line.replace("-", "=")
            new_report.append(line)
        elif line.startswith("+="):
            line = line.replace("+", "|")
            line = line.replace("=", "-")
            new_report.append(line)
        elif line.startswith("+-"):
            continue
        else:
            new_report.append(line)

    print("PYLINT REPORT    ")
    print("    ")
    print(summary, "    ")
    print("    ")
    print("<details>")
    print("  <summary>Full report</summary>")
    print("    ")
    for line in new_report:
        print("  ", line, end="")
    print("</details>")


if __name__ == "__main__":
    main(*sys.argv[1:])
