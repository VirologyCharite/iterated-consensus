# Iterated consensus sequence

Iteratively build a consensus sequence: call a consensus from a BAM
(or an initial reference), build a mapper index from it, remap the
original reads against it, call a new consensus, and repeat until the
consensus stops changing.

Mapping and consensus-calling are both fully user-configurable via a
TOML pipeline config -- you bring your own
`bowtie2`/`bwa`/`minimap2`/... and `ivar`/`samtools`/... commands,
iterated-consensus just drives the loop, tracks convergence, and
writes the results.

The code here was written entirely by Claude Sonnet 5 (model ID
`claude-sonnet-5`), from Anthropic's Claude 5 family. This took about
3.5 hours from giving Claude an initial description, through planning
and writing code, tests, and documentation.

## Install

```
uv add iterated-consensus     # or: pip install iterated-consensus
```

## Quick start

```
iterated-consensus config-template                  # list bundled presets
iterated-consensus config-template bowtie2-ivar > pipelines.toml
# edit pipelines.toml: fill in [input], adjust commands/threads as needed
iterated-consensus run --config pipelines.toml --out-dir results/ --dry-run  # preview
iterated-consensus run --config pipelines.toml --out-dir results/ --progress
```

`--progress` prints a one-line summary after each iteration (reads mapped,
consensus length, identity to the previous consensus, time taken). Every run
also writes `results/index.html` -- open it in a browser for a summary and
full per-iteration detail, no `--progress` needed.

Errors (bad config, a failing command, a mismatched reference, ...) normally
print a short `error: ...` message. Add `--traceback` to instead let them
crash with the full Python traceback, for debugging.

## Config format

A config has one or more `[[mapper]]` tables, one `[consensus]` table, and
optionally `[input]` and `[run]`.

```toml
[[mapper]]
name = "bowtie2"
index_cmd = ["bowtie2-build", "{reference}", "{index_prefix}"]
map_cmd = "bowtie2 -x {index_prefix} -1 {reads_1:,} -2 {reads_2:,} -p {threads} | samtools sort -o {bam}"

# A second [[mapper]] table can be added too -- if more than one mapper is
# configured, every mapper runs each iteration and their BAMs are merged
# before the consensus step sees them.

[consensus]
steps = [
    "samtools mpileup -aa -A -d 0 -Q 0 -f {reference} {bam} | ivar consensus -p {consensus_prefix} -t 0.5",
]
output = "{consensus_prefix}.fa"   # where the harness should find the result

[input]
mate1 = ["a_R1.fastq.gz", "b_R1.fastq.gz"]   # 0 or more paired sets
mate2 = ["a_R2.fastq.gz", "b_R2.fastq.gz"]
unpaired = []                                 # 0 or more single-end files
reference_fasta = "starting_reference.fasta"  # a local file...
# reference_id = "chr2"        # ...and/or a name -- see "Reference resolution" below
# --- or, instead of the FASTQ block above, start from a BAM: ---
# bam = "input.bam"
# reference_id = "chr2"        # only needed if the BAM has >1 reference
# reference_fasta = "chr2.fasta"  # optional -- see "Reference resolution" below
# bam_reads = "ref"            # ref | ref+unal | all -- see below

[run]
threads = 4
max_iterations = 20
convergence_identity = 100.0   # stop once consensus identity to the previous
                                # iteration reaches this percent...
convergence_streak = 1         # ...for this many iterations in a row
```

`[input]` can instead (or partly) be supplied on the command line -- see
`iterated-consensus run --help`. CLI values override the config's `[input]`
field-by-field, so a config can be fully self-contained or left generic and
pointed at different data per invocation.

### Command steps: list or shell string

Every command (`index_cmd`, `map_cmd`, each `[consensus]` step) can be
written as a list of argv tokens (run directly, no shell -- safest, use this
whenever you're just running one program) or as a single string (run via a
shell -- needed for pipes, as in the `samtools mpileup | ivar consensus`
example above).

### Placeholders

- `{reference}` -- the current reference FASTA: the previous iteration's
  consensus, or the starting reference for iteration 0. For a BAM-start run,
  iteration 0's reference is only available if it could be resolved -- see
  "Reference resolution" below; if not, and `[consensus]` uses `{reference}`
  anyway, that's a config error caught before anything runs, not a silent
  guess.
- `{index_prefix}` -- path prefix for this mapper's index this iteration.
- `{bam}` -- path this mapper should write its BAM to.
- `{consensus_prefix}` -- path prefix for the consensus step's output.
- `{threads}` -- from `[run]` threads.
- `{reads_1}`, `{reads_2}`, `{reads_single}` -- the read-file lists (mate1,
  mate2, unpaired). Only present if that category is non-empty for this run
  -- referencing e.g. `{reads_single}` in a run with no unpaired reads is a
  config error, so a mapper template should only reference the categories it
  actually expects.

### Read-list expansion syntax

A read-list placeholder (`reads_1`, `reads_2`, `reads_single`) can be
written plain or with modifiers, to match whatever multi-file syntax your
mapper wants:

| Form | Expands to |
|---|---|
| `{reads_1}` | space-joined: `f1.fq f2.fq f3.fq` |
| `{reads_1:,}` | joined with a literal separator: `f1.fq,f2.fq,f3.fq` |
| `{-1:reads_1}` | prefix before each file: `-1f1.fq -1f2.fq -1f3.fq` |
| `{-1 :reads_1}` | prefix (here with a trailing space) before each file: `-1 f1.fq -1 f2.fq -1 f3.fq` |
| `{-1:reads_1:,}` | prefix + separator together: `-1f1.fq,-1f2.fq,-1f3.fq` |
| `{cat:reads_1}` | concatenates all files into one and substitutes its path -- for mappers (e.g. `bwa mem`, `minimap2`) that only accept exactly one file per mate |

Whether the first colon-separated part is a prefix or the list name itself
is inferred from whether it names a known read list. `{cat:name}` is a
reserved special case in the prefix position -- concatenation only happens
once per run (not per iteration) and only for mappers whose template
actually uses `{cat:...}`.

### `bam_reads`: which reads to use when starting from a BAM

Every iteration remaps the same original read pool (extracted once, reused
throughout) -- it never shrinks to just whatever mapped last time. When that
pool comes from an input BAM rather than FASTQ files, `bam_reads` controls
its scope:

- `ref` (default, strictest) -- only reads aligned to the chosen reference/contig.
- `ref+unal` -- that, plus reads that didn't map anywhere (candidates for
  mapping once the consensus improves).
- `all` -- every read in the BAM, regardless of what it mapped to.

### Reference resolution

Two `[input]` keys between them cover every way of specifying a starting
reference, for both FASTQ-start and BAM-start:

- `reference_id` -- just a *name*, never a file. It can be a record ID
  within `reference_fasta`, a contig name within `bam`, and/or an NCBI
  accession -- the same string can serve more than one of these roles at
  once (see the examples below).
- `reference_fasta` -- a pre-existing local FASTA file. If it has exactly
  one sequence, that sequence is used automatically; if it has more than
  one, `reference_id` must be given to pick which.

Whether you need one, the other, both, or neither depends on the situation:

- **Neither**: BAM-start, the BAM has only one reference, and its name is
  itself an NCBI accession (e.g. `NC_045512.2`) -- it's fetched
  automatically.
- **`reference_id` only**: BAM-start, the BAM has several references, you
  pick one with `reference_id`, and that name is itself an accession.
- **`reference_fasta` only**: FASTQ-start (or BAM-start), and the file has
  just one sequence.
- **Both**: FASTQ-start (or BAM-start) with a multi-sequence
  `reference_fasta` -- `reference_id` picks which record. For BAM-start
  specifically, this is also how you'd give the actual reference the BAM
  was aligned against, rather than relying on auto-fetch.
- **Neither, and unresolvable**: iteration 0 just has no `{reference}`.
  That's fine for a BAM-start run *unless* `[consensus]` actually uses
  `{reference}` -- in which case it's reported as a config error before
  anything runs, since that combination can never succeed (iteration 0
  always runs first, before any consensus this tool computed exists to fall
  back on). A FASTQ-start run always needs *some* reference to build
  iteration 0's mapping index against, so this case is always an error
  there.

For a BAM-start run, whichever way a reference is obtained, it's validated
against the BAM before use: the FASTA record's id must match the resolved
contig name exactly, and its sequence length must match the BAM header's
length for that contig exactly. A mismatch aborts the run immediately --
proceeding would mean calling a pileup-based consensus against a reference
the BAM's own coordinates don't actually match, silently producing garbage.

## Output

Iterations are numbered from 0. `iter_000` is always the bootstrap step: it
produces the first consensus from whatever mapping was already available at
the start, rather than one this tool built by iterating. Its exact shape
depends on how the run started:

- FASTQ-start: `iter_000` builds an index from the reference you gave, maps
  the reads against it, and calls a consensus -- same shape as every later
  iteration, just against a reference you supplied rather than one this tool
  computed.
- BAM-start: `iter_000` calls a consensus directly from the input BAM, with
  no mapping step -- the alignment already exists.

From `iter_001` onward, every iteration has the same shape regardless of how
the run started: build an index from the previous iteration's
`consensus.fasta`, remap the (always-the-same, extracted-once) reads against
it, and call a new consensus. `iter_001` is therefore always the first
iteration with an `identity_to_previous` value, since `iter_000` has nothing
before it to compare against.

```
results/
  reads/                     extracted/concatenated read files (built once)
  reference_initial.fasta    normalized starting reference (FASTQ-start only)
  iter_000/
    <mapper>_index.*         index files (FASTQ-start only -- see above)
    <mapper>.bam              that mapper's mapping output (FASTQ-start only)
    merged.bam                 (only if >1 mapper) merged BAM the consensus step sees
    consensus.fasta             this iteration's consensus
    stats.json                  reads mapped, length, identity to previous, base composition
    logs/                        captured stdout+stderr of every command run
  iter_001/
    <mapper>_index.*         index files, one set per configured mapper
    <mapper>.bam              that mapper's mapping output
    merged.bam                 (only if >1 mapper) merged BAM the consensus step sees
    consensus.fasta             this iteration's consensus
    stats.json                  reads mapped, length, identity to previous, base composition
    logs/                        captured stdout+stderr of every command run
  iter_002/
    ...
  metrics.tsv                 one row per iteration
  summary.json                 iterations run, converged?, total time
  index.html                    human-readable report rendered from summary.json
```

## Resuming a run

If a run stops because it hit `max_iterations` without converging, raise
`max_iterations` in the config and rerun the exact same command (same
`--config`, same `--out-dir`): it picks up from the next iteration rather
than starting over. This is detected automatically from `summary.json` in
`--out-dir` -- there's no separate flag. Already-extracted reads and
already-completed iterations aren't redone.

If a run already converged, rerunning it against the same `--out-dir` is a
no-op: it reports the existing result without redoing anything.

Resume trusts that `--out-dir` corresponds to the same logical run --
pointing it at a config with a different mapper, consensus pipeline, or
input isn't validated or rejected, it'll just continue on top of whatever's
there. Also, resuming is driven entirely by `summary.json`; a run that
crashed before writing it (e.g. killed mid-iteration) can't be resumed and
should be started fresh in a new `--out-dir`.

## Known limitations

- `--dry-run` shows iterations 0 and 1 in full (both always run), but can't
  show iteration 2 onward -- those commands depend on files that don't exist
  yet, and whether the run even reaches them depends on convergence.
- A mapper's command template must match the input categories actually
  present for a given run (e.g. don't reference `{reads_single}` in a config
  meant to also run on paired-only data). There's no conditional templating
  -- write separate configs for meaningfully different input shapes.

## Development

```
uv sync
uv run pytest
```
