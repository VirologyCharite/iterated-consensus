"""Test helper: fake 'consensus caller' that writes a fixed sequence, ignoring input.

Used for the BAM-start integration test, where iteration 0 has no {reference}
to copy from.
"""

import sys

FIXED_SEQUENCE = "ACGTACGTACGTACGTACGT"


def main() -> None:
    out_path = sys.argv[1]
    with open(out_path, "w") as f:
        f.write(f">seed\n{FIXED_SEQUENCE}\n")


if __name__ == "__main__":
    main()
