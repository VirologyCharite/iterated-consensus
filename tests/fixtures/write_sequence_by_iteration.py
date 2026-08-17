"""Test helper: fake 'consensus caller' whose output depends on the iteration
number embedded in its own output path (e.g. .../iter_002/consensus.fa).

Used to build a run that provably does NOT converge for the first few
iterations, then does -- so resume behavior (continuing past a max_iterations
cutoff) can be tested deterministically.
"""

import re
import sys

SEQUENCES_BY_ITERATION = {0: "AAAA", 1: "CCCC"}
STABLE_SEQUENCE = "GGGG"  # every iteration from 2 onward emits this


def main() -> None:
    out_path = sys.argv[1]
    match = re.search(r"iter_(\d+)", out_path)
    iteration = int(match.group(1)) if match else 0
    sequence = SEQUENCES_BY_ITERATION.get(iteration, STABLE_SEQUENCE)
    with open(out_path, "w") as f:
        f.write(f">c\n{sequence}\n")


if __name__ == "__main__":
    main()
