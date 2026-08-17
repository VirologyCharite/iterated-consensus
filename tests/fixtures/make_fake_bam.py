"""Test helper: write a tiny valid single-read BAM to the given path.

Used as a fake mapper's map_cmd in integration tests, so runner tests don't
need a real aligner installed.
"""

import sys

import pysam


def main() -> None:
    out_path = sys.argv[1]
    header = {"HD": {"VN": "1.6"}, "SQ": [{"LN": 100, "SN": "dummy"}]}
    with pysam.AlignmentFile(out_path, "wb", header=header) as f:
        segment = pysam.AlignedSegment()
        segment.query_name = "r1"
        segment.query_sequence = "ACGT"
        segment.flag = 0
        segment.reference_id = 0
        segment.reference_start = 0
        segment.mapping_quality = 60
        segment.cigar = [(0, 4)]
        segment.query_qualities = pysam.qualitystring_to_array("IIII")
        f.write(segment)


if __name__ == "__main__":
    main()
