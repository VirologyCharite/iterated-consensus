"""Test helper: append each of argv[2:] as its own line to the file at argv[1].

Used to verify placeholder substitution in commands (e.g. [output].commands)
without depending on shell quoting or any real external tool.
"""

import sys


def main() -> None:
    log_path = sys.argv[1]
    with open(log_path, "a") as f:
        for value in sys.argv[2:]:
            f.write(value + "\n")


if __name__ == "__main__":
    main()
