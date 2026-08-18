"""Test helper: write a tiny valid but *unsorted* BAM to the given path.

Used as a fake mapper's map_cmd, to test that the runner's automatic
`.bam.bai` creation surfaces a clear error for a mapper that forgot to sort
its output, rather than a raw pysam traceback.
"""

import sys

import pysam


def _segment(name: str, pos: int) -> pysam.AlignedSegment:
    s = pysam.AlignedSegment()
    s.query_name = name
    s.query_sequence = "ACGT"
    s.flag = 0
    s.reference_id = 0
    s.reference_start = pos
    s.mapping_quality = 60
    s.cigar = [(0, 4)]
    s.query_qualities = pysam.qualitystring_to_array("IIII")
    return s


def main() -> None:
    out_path = sys.argv[1]
    header = {"HD": {"VN": "1.6"}, "SQ": [{"LN": 100, "SN": "dummy"}]}
    with pysam.AlignmentFile(out_path, "wb", header=header) as f:
        # coordinate order 10, then 0 -- not sorted
        f.write(_segment("r10", 10))
        f.write(_segment("r0", 0))


if __name__ == "__main__":
    main()
