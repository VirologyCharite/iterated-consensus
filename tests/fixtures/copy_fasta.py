"""Test helper: fake 'consensus caller' that just copies a fasta file.

Used in integration tests so the loop converges deterministically without a
real consensus tool installed.
"""

import shutil
import sys


def main() -> None:
    shutil.copyfile(sys.argv[1], sys.argv[2])


if __name__ == "__main__":
    main()
