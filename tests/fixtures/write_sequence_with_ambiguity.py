"""Test helper: fake 'consensus caller' whose output depends on the iteration
number embedded in its own output path (e.g. .../iter_001/consensus.fa).

Used to exercise ambiguity-transition tracking: the iteration-0 sequence has
an ambiguous site that's resolved by iteration 1, and iteration 1 introduces
a different ambiguous site that iteration 0 didn't have.
"""

import re
import sys

SEQUENCES_BY_ITERATION = {0: "ACGNACGT", 1: "ACGTACGN"}
STABLE_SEQUENCE = "ACGTACGN"  # every iteration from 2 onward emits this


def main() -> None:
    out_path = sys.argv[1]
    match = re.search(r"iter_(\d+)", out_path)
    iteration = int(match.group(1)) if match else 0
    sequence = SEQUENCES_BY_ITERATION.get(iteration, STABLE_SEQUENCE)
    with open(out_path, "w") as f:
        f.write(f">c\n{sequence}\n")


if __name__ == "__main__":
    main()
