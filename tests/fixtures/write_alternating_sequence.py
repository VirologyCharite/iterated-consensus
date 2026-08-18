"""Test helper: fake 'consensus caller' that alternates between two different
sequences based on the iteration number embedded in its own output path
(e.g. .../iter_002/consensus.fa) -- a deterministic period-2 cycle, for
testing cycle detection.
"""

import re
import sys

SEQUENCE_A = "AAAACCCC"
SEQUENCE_B = "GGGGTTTT"


def main() -> None:
    out_path = sys.argv[1]
    match = re.search(r"iter_(\d+)", out_path)
    iteration = int(match.group(1)) if match else 0
    sequence = SEQUENCE_A if iteration % 2 == 0 else SEQUENCE_B
    with open(out_path, "w") as f:
        f.write(f">c\n{sequence}\n")


if __name__ == "__main__":
    main()
