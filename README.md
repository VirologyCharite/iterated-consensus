# Iterated consensus sequence

*NOTE* This package is only a day old, so you should expect occasional
breaking changes!

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
iterated-consensus run --config pipelines.toml --output-dir results/ --dry-run  # preview
iterated-consensus run --config pipelines.toml --output-dir results/ --progress
```

`--progress` prints a one-line summary after each iteration (reads mapped,
consensus length, number of ambiguous characters in the consensus, identity
to the previous consensus, time taken, and the consensus sequence's MD5 --
handy for confirming at a glance that two iterations (or two separate runs)
produced byte-for-byte the same sequence).
Every run also writes `results/index.html` -- open it in a browser for a
summary and full per-iteration detail, no `--progress` needed.

`--output-dir` can be left off the command line if the config sets
`[run].output_dir` instead (the flag wins if both are given) -- handy for a
config that's meant to always write to the same place. Giving neither is a
config error.

Errors (bad config, a failing command, a mismatched reference, ...) normally
print a short `error: ...` message. Add `--traceback` to instead let them
crash with the full Python traceback, for debugging.

## Config format

A config has one or more `[[mapper]]` tables, one `[consensus]` table, and
optionally `[input]`, `[output]`, and `[run]`.

```toml
[[mapper]]
name = "bowtie2"
index_cmd = ["bowtie2-build", "{reference}", "{index_prefix}"]
map_cmd = "bowtie2 -x {index_prefix} -1 {reads_1:,} -2 {reads_2:,} -p {threads} | samtools sort -o {bam}"

# A second [[mapper]] table can be added too -- if more than one mapper is
# configured, every mapper runs each iteration and their BAMs are merged
# before the consensus step sees them.

[mapper.tool-versions]
# Optional -- see "Tool versions" below. Attaches to the [[mapper]] table
# immediately above it.
bowtie2 = "bowtie2 --version | head -n 1 | cut -f3 -d' '"

[consensus]
steps = [
    "samtools mpileup -aa -A -d 0 -Q 0 -f {reference} {bam} | ivar consensus -p {consensus_prefix} -t 0.5",
]
output = "{consensus_prefix}.fa"   # where to find the result -- see note below

[consensus.tool-versions]
ivar = "ivar version | head -n 1 | cut -f3 -d' '"

[input]
reads_1 = ["a_R1.fastq.gz", "b_R1.fastq.gz"]  # 0 or more paired sets
reads_2 = ["a_R2.fastq.gz", "b_R2.fastq.gz"]
reads_single = []                             # 0 or more single-end files
reference_fasta = "starting_reference.fasta"  # a local file...
# reference_id = "chr2"        # ...and/or a name -- see "Reference resolution" below
# --- or, instead of the FASTQ block above, start from a BAM: ---
# bam = "input.bam"
# reference_id = "chr2"        # only needed if the BAM has >1 reference
# reference_fasta = "chr2.fasta"  # optional -- see "Reference resolution" below
# bam_reads = "ref"            # ref | ref+unal | all -- see below

[output]
# consensus_fasta = "final_consensus.fasta"   # copy the last iteration's
#   consensus here once the run finishes, converged or not -- a relative
#   path here is relative to the cwd, NOT --output-dir; use
#   "{output_dir}/final_consensus.fasta" to put it under --output-dir. See
#   below.
# consensus_id = "my-sample-name"              # optional new FASTA header
# final_reference_fasta = "{output_dir}/final_reference.fasta"  # symlink to
#   the reference used for the final alignment step -- see below
# final_reference_bam = "{output_dir}/final_reference.bam"      # symlink to
#   the BAM mapped against it and used to call the final consensus
# commands = ["samtools faidx {consensus_fasta}"]   # optional -- run once
#   consensus_fasta/final_reference_fasta/final_reference_bam are written;
#   {consensus_fasta}/{consensus_id}/{final_reference_fasta}/
#   {final_reference_bam} available

[run]
# output_dir = "results/"          # fallback for --output-dir; the CLI flag
                                    # wins if both are given
threads = 4                        # or "auto" -- see "Threads and custom
                                    # [run] variables" below
max_iterations = 20
convergence_identity = 100.0       # stop once consensus identity to the
                                    # previous iteration reaches this percent...
convergence_streak = 1             # ...for this many iterations in a row

# Anything else here becomes a {name} placeholder in every command, e.g.:
# min_depth = 10               # -> {min_depth} in [consensus] steps
# sample_name = "patient-42"   # -> {sample_name} anywhere
```

`[input]` can instead (or partly) be supplied on the command line -- see
`iterated-consensus run --help`. CLI values override the config's `[input]`
field-by-field, so a config can be fully self-contained or left generic and
pointed at different data per invocation.

`[consensus].output` is only used right after the steps run, to find and
read whatever file your tool actually wrote (different tools name it
differently -- `ivar consensus -p PREFIX` writes `PREFIX.fa`, which is why
the example above is `output = "{consensus_prefix}.fa"`, not
`{consensus_prefix}` alone). What gets read there is then copied to this
iteration's own `consensus.fasta` (see "Output" below) -- and it's *that*
fixed-name copy, not the path `output` pointed to, that later becomes
`{reference}` for the next iteration. This is why a mapper's `index_cmd` in
`--dry-run` output references `consensus.fasta` even if your `output`
pattern produces a different filename or extension: `output` only has to
match what your consensus tool actually writes, nothing downstream reads
that path directly.

### Tool versions: `[mapper.tool-versions]` and `[consensus.tool-versions]`

Optional sub-tables -- `[mapper.tool-versions]` attaches to whichever
`[[mapper]]` table comes immediately before it (so with more than one
mapper, give each its own); `[consensus.tool-versions]` attaches to
`[consensus]`. Each entry maps a name you choose to a command whose stdout
is that tool's version:

```toml
[mapper.tool-versions]
bowtie2 = "bowtie2 --version | head -n 1 | cut -f3 -d' '"

[consensus.tool-versions]
ivar = "ivar version | head -n 1 | cut -f3 -d' '"
```

Only stdout is captured (not stderr -- add your own `2>&1` if a tool prints
its version there instead), and it's kept exactly as printed, whitespace
trimmed from each end but otherwise unchanged -- multi-line output is fine.
Each iteration, right before that mapper's `index_cmd`/`map_cmd` run (or,
for `[consensus]`, right before its steps run), every configured
tool-versions command runs and its output is recorded in that iteration's
`stats.json`. A mapper's tool-versions only run in iterations where that
mapper actually runs -- e.g. never for `iter_000` of a BAM-start run, which
has no mapping step at all.

The report (`index.html`) shows each tool's version once, near the top of
its "Logs" section (see "Output" below) -- unless it *changed* partway
through the run, in which case that's flagged prominently instead, listing
every distinct version seen and which iteration(s) reported it. This is the
point of the feature: confirming the same tool binary was used from start
to finish, since a version drifting mid-run (a PATH change, an unpinned
container tag, a background upgrade) can quietly invalidate a comparison
between early and late iterations.

### Command steps: list or shell string

Every command (`index_cmd`, `map_cmd`, each `[consensus]` step) can be
written as a list of argv tokens (run directly, no shell -- safest, use this
whenever you're just running one program) or as a single string (run via a
shell -- needed for pipes, as in the `samtools mpileup | ivar consensus`
example above).

After each mapper's `map_cmd` runs (and after merging, if more than one
mapper is configured), iterated-consensus checks for a `.bam.bai` index next
to the resulting BAM and creates one with `samtools index` if it's missing
-- most consensus tools need one, so you don't have to remember to add an
indexing step to `map_cmd` yourself. Indexing needs the BAM coordinate-sorted
first, so the BAM header's own `SO` tag is checked before doing anything
else: if it's already marked `SO:coordinate`, indexing runs directly (fast,
and avoids a pointless sort pass over a large, already-sorted BAM); if not
(`map_cmd` forgot a `samtools sort`, or a mapper wrote it unsorted), it's
sorted in place automatically first, so `map_cmd` doesn't strictly need its
own sort step either. (A BAM whose header lies about being sorted -- rare,
but not impossible -- is still caught: if indexing the "already-sorted" file
fails anyway, it's sorted for real and indexing is retried.) This only
covers the BAM a mapping step actually produces; it doesn't apply to
`iter_000` of a BAM-start run, which has no mapping step (see "Reference
resolution"
above) -- add your own `samtools index {bam}` step to `[consensus]` if that
BAM needs indexing too (the bundled `bwa-samtools` preset does exactly
this).

### Placeholders

- `{reference}` -- the current reference FASTA: the previous iteration's
  consensus, or the starting reference for iteration 0. For a BAM-start run,
  iteration 0's reference is only available if it could be resolved -- see
  "Reference resolution" below; if not, and `[consensus]` uses `{reference}`
  anyway, that's a config error caught before anything runs, not a silent
  guess.
- `{index_prefix}` -- path prefix for this mapper's index this iteration.
- `{bam}` -- path this mapper should write its BAM to.
- `{consensus_prefix}` -- path prefix for the consensus step's output (see
  the note on `[consensus].output` above -- this is not the same file that
  later becomes `{reference}`).
- `{threads}` -- from `[run]` threads. See "Threads and custom `[run]`
  variables" below for `threads = "auto"` and defining your own placeholders
  alongside it (e.g. for splitting a thread budget across a piped command).
- `{output_dir}` -- the effective output directory: `--output-dir` if given,
  else `[run].output_dir`. Useful in `[output].consensus_fasta` and
  `[output].commands` (see "Final output: `[output]`" below) to put the
  final deliverable, or something derived from it, under the same directory
  as everything else the run produced.
- `{reads_1}`, `{reads_2}`, `{reads_single}` -- the read-file lists, matching
  `[input]`'s `reads_1`/`reads_2`/`reads_single` one-for-one. Only present if
  that category is non-empty for this run -- referencing e.g.
  `{reads_single}` in a run with no unpaired reads is a config error, so a
  mapper template should only reference the categories it actually expects.

### Threads and custom `[run]` variables

Every `[run]` key is available as a `{name}` placeholder in every
`index_cmd`, `map_cmd`, and `[consensus]` step (and, for `{output_dir}`
especially, in `[output].consensus_fasta`/`commands` too -- see "Final
output: `[output]`") -- both the ones with dedicated meaning (`{threads}`,
`{max_iterations}`, `{convergence_identity}`, `{convergence_streak}`,
`{threads_reserve}`, `{output_dir}`) and any custom ones you add. This is
handy for more than just thread counts -- e.g. a logging step that records
what a run was configured with:

```toml
[run]
threads = 8
min_depth = 10
sample_name = "patient-42"
```

```toml
[consensus]
steps = [
    "samtools mpileup -d 0 {bam} | ivar consensus -t 0.5 -m {min_depth} -p {consensus_prefix}",
    'echo "Built {sample_name} consensus at {threads} threads, target {convergence_identity}% over {max_iterations} iterations" >> log.txt',
]
```

A custom variable can't reuse a name iterated-consensus already sets itself
(the dedicated `[run]` fields above, plus the per-iteration placeholders
`reference`, `index_prefix`, `bam`, `consensus_prefix`, `reads_1`, `reads_2`,
`reads_single`, `consensus_fasta`, `consensus_id`, `final_reference_fasta`,
`final_reference_bam`) -- that's rejected at config load time rather than
silently shadowed.

**`threads = "auto"`** resolves to the number of CPUs actually available to
the process (respecting container/cgroup/`taskset` limits on Linux, where
that's exposed; the installed core count elsewhere) at config-load time, once
per run -- not re-detected per iteration. Pair it with `threads_reserve` (an
integer, only valid alongside `threads = "auto"`) to leave some cores free
for other work on the machine:

```toml
[run]
threads = "auto"
threads_reserve = 2   # use (detected CPUs - 2), never less than 1
```

**Splitting a thread budget across a pipe.** `{threads}` is one number, but a
piped command like `bwa mem | samtools sort` or
`samtools mpileup | ivar consensus` runs two programs *concurrently*, each of
which could use its own thread count -- and since they're running at the same
time, those counts add up against your actual core count, they don't each
get to use the whole budget. Using `{threads}` unmodified for more than one
stage of the same pipe oversubscribes the machine. iterated-consensus doesn't
try to auto-split `{threads}` for you -- pipeline stages have wildly
different threading characteristics (some don't support it at all, some
scale linearly, some plateau early), so a generic split would just be a
guess. Instead, treat `{threads}` (or a `threads = "auto"` budget) as the
total, and partition it yourself into named `[run]` variables that add up to
no more than that:

```toml
[run]
threads = "auto"
threads_reserve = 2   # e.g. resolves to 6 on an 8-core machine
map_threads = 5       # give most of the budget to the mapper...
sort_threads = 1       # ...and a thread or two to samtools sort running
                        # alongside it in the same pipe

[[mapper]]
name = "bwa"
map_cmd = "bwa mem -t {map_threads} {index_prefix} {cat:reads_1} {cat:reads_2} | samtools sort -@ {sort_threads} -o {bam}"
```

`{threads}` is still the right placeholder for any command that's just one
program (most `index_cmd`s, or a `[consensus]` step with no pipe) -- reach
for named variables like `map_threads`/`sort_threads` only where a pipe
means two or more programs are genuinely running at once.

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

Your input BAM itself is never rewritten, regardless of any of this. Some of
the above needs it indexed, which needs it coordinate-sorted; if it's
already sorted, an index is created directly beside it if missing (that's
harmless -- purely additive, same as any tool would do), but if it actually
needs sorting, that happens to a separate copy under `--output-dir`, not to
your file. (This is specifically about the BAM you pass in via `bam =`;
BAMs iterated-consensus generates itself, like each iteration's mapping
output, are sorted in place freely -- those are its own working files.)

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

### Convergence: `convergence_identity` and `convergence_streak`

Every iteration from 1 onward computes the identity between its new
consensus and the previous one (iteration 0 has nothing to compare against,
so it's skipped). That identity feeds a streak counter: it increments
whenever `identity >= convergence_identity`, and resets to 0 otherwise. The
run stops, reported as converged, once the streak reaches
`convergence_streak` -- i.e. once `convergence_streak` *consecutive*
iterations have each been at or above `convergence_identity`. If that never
happens, the run stops anyway once `max_iterations` is reached, just
reported as not converged.

With the defaults (`convergence_identity = 100.0`, `convergence_streak =
1`), a run stops as soon as one iteration produces a consensus identical to
the one before it. Raising `convergence_streak` (e.g. to 2 or 3) guards
against declaring convergence on a fluke -- a sequence could hit exactly
100% once by chance (e.g. a low-coverage region that happens to resolve to
the same majority base) and then drift again next iteration; requiring
several consecutive matches is a stronger signal that it's genuinely
settled. Lowering `convergence_identity` below 100 is also legitimate --
useful for a sample that may never perfectly stabilize (e.g. a genuinely
heterogeneous/mixed population), where "close enough" is a more realistic
stopping condition than exact equality.

### Cycle detection

Some inputs never settle on a single consensus: the majority call at one or
more sites flips back and forth (e.g. a near-50/50 heterozygous position),
so the run oscillates among a small set of sequences instead of converging.
Left alone, that would just burn through every iteration up to
`max_iterations`. Instead, every new iteration's consensus MD5 (see
`consensus_md5` in "Output" below) is checked against every *earlier*
iteration's, not just the immediately preceding one -- if it matches, the
run stops immediately, reported as a detected cycle rather than as
converged or as having hit `max_iterations`.

The *first* occurrence of the repeated sequence is treated as the result --
it's what `[output].consensus_fasta`/`commands` run against, and what
`index.html`'s "Final ..." stats reflect -- not the later iteration that
happened to repeat it, since the earlier one is exactly as valid a
representative of the cycle and arrived at with less noise accumulated
along the way. The iterations that actually ran (including the repeat) are
still all recorded normally, so the cycle itself is visible: in
`--progress`/the per-iteration table via the repeating `consensus_md5`
values, and as its own prominent notice at the top of `index.html`.

A cycle of period 1 (immediately repeating the previous iteration exactly)
isn't a special case at all -- that's just ordinary convergence with the
default `convergence_identity = 100.0`, handled the normal way. Cycle
detection specifically covers longer periods (2 or more) that the
adjacent-iteration `identity_to_previous` check can't see. Resuming a run
that stopped this way re-detects the same cycle (from the recorded history)
rather than continuing to iterate.

### Final output: `[output]`

Everything a run produces lives under `--output-dir` regardless, findable as
`iter_NNN/consensus.fasta` for whichever iteration ran last (see "Output"
below) -- `[output]` is an optional convenience on top of that, for when you
want the final result copied somewhere specific rather than having to know
which `iter_NNN` was the last one.

- `consensus_fasta` -- where to copy it to. This is a *template*, not a
  plain path: it's rendered against the same placeholders every other
  command gets (`{threads}`, any custom `[run]` variables, and
  `{output_dir}` -- see "Threads and custom `[run]` variables"), then used
  as-is. A relative result is **not** auto-placed under `--output-dir`: it's
  just a normal relative path, resolved against the current working
  directory like any other file argument. Write
  `consensus_fasta = "{output_dir}/final_consensus.fasta"` to put it
  alongside everything else the run produced; leave `{output_dir}` out to
  get a plain cwd-relative (or absolute) path instead. Either way, its
  parent directory is created if missing.
- `consensus_id` -- the FASTA header to give the copy. Optional; if
  omitted, it keeps whatever id the consensus tool itself assigned.
  Requires `consensus_fasta` to also be given -- renaming with nowhere to
  write doesn't mean anything on its own.
- `final_reference_fasta` -- where to put a copy of the reference used for
  the final alignment step, i.e. the second-last iteration's
  `consensus.fasta` (the last iteration's own `consensus.fasta` is the
  *result*, not something it was aligned against). Same templating rules as
  `consensus_fasta`. Unlike `consensus_fasta`, this is always a *relative
  symlink* into `--output-dir`'s own `iter_NNN/consensus.fasta`, never a
  copy: renaming it would desynchronize it from the reference name embedded
  in `final_reference_bam`'s BAM header (see below), so it's left exactly as
  iteration numbering produced it. If a `.fai` already sits next to that
  `consensus.fasta` (some `[consensus]` pipelines index the reference as a
  side effect, e.g. `samtools mpileup -f`), a matching `.fai` symlink is
  created alongside it too -- unlike `final_reference_bam`'s `.bam.bai`
  (always present, since mapping always indexes its own BAM), a source
  `.fai` isn't guaranteed to exist, so it's skipped rather than left
  dangling when there isn't one. Independent of `consensus_fasta` -- you can
  set either, both, or neither. There's always a second-last iteration to
  point to, since at least two iterations (`iter_000` and `iter_001`) always
  run.
- `final_reference_bam` -- where to put a copy of the BAM that was actually
  mapped against `final_reference_fasta` and used to call the final
  consensus. This lives in the *final* iteration's own directory (one
  iteration dir apart from `final_reference_fasta`, since a BAM sits
  alongside the consensus it produced, not the reference it was mapped
  against). Also always a relative symlink, with a matching `.bam.bai`
  symlink created alongside it automatically. Independent of
  `consensus_fasta`/`final_reference_fasta` too.
- `commands` -- steps (same `index_cmd`/`map_cmd`/`[consensus]` steps
  syntax: a list, each a list of argv tokens or a shell string) to run once
  `consensus_fasta` has been written. Requires `consensus_fasta` to also be
  given. Extra placeholders are available here on top of the usual ones
  (`{threads}` and friends, any custom `[run]` variables -- see "Threads and
  custom `[run]` variables"), resolved to whichever of `consensus_fasta`/
  `final_reference_fasta`/`final_reference_bam` were actually configured:
  `{consensus_fasta}` (the path just written), `{consensus_id}` (the id that
  ended up in it -- the resolved value, so it's set even when `consensus_id`
  itself was left unset in the config), `{final_reference_fasta}`, and
  `{final_reference_bam}` (the symlink paths just created).
  Logs go to `--output-dir/logs/output_command_NN.log`. A failing command
  aborts the run the same way a failing mapper/consensus step would.

This always uses whichever iteration ran *last* -- converged or not, since
even a run that hit `max_iterations` without converging usually still has a
usable "current best" consensus worth having on hand (if a cycle was
detected instead -- see "Cycle detection" below -- the first occurrence of
the repeated consensus is used, not the last iteration actually run). Every
configured deliverable (`consensus_fasta`, `final_reference_fasta`,
`final_reference_bam`) and `commands` happen at the end of every `run()`
call, including a resumed run that turns out to already be converged -- so
re-running the same command is always safe (though note `commands`
therefore also re-runs every time, e.g. on a repeated already-converged
`run()` call -- keep that in mind for a command with side effects
elsewhere, like an upload). Unlike `reference_initial.fasta` (see
"Reference resolution" above), `consensus_fasta`'s copy is always a real
copy, never a symlink: it's meant to be a small, standalone, portable
deliverable, not a pointer back into `--output-dir`'s own working files.
`final_reference_fasta`/`final_reference_bam` are the opposite: always
symlinks, since they're meant to point back at those exact working files
(see above).

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

Every iteration's `consensus.fasta` has `-iteration-N` appended to its FASTA
id before it's written, so each one is uniquely named even when the
underlying sequence repeats across iterations (e.g. after convergence) --
that's what keeps a later iteration's BAM header (which embeds whatever
reference name it was mapped against) matching up with
`[output].final_reference_fasta`'s own id (see "Final output: `[output]`"
above). The one exception is `[output].consensus_fasta` itself, whose
default id has this suffix stripped back off, so it still defaults to
whatever id your consensus tool originally assigned.

```
results/
  reads/                     extracted/concatenated read files (built once)
  reference_initial.fasta    starting reference (FASTQ-start always; BAM-start
                             if one was resolved) -- see note below
  iter_000/
    <mapper>_index.*         index files (FASTQ-start only -- see above)
    <mapper>.bam              that mapper's mapping output (FASTQ-start only)
    merged.bam                 (only if >1 mapper) merged BAM the consensus step sees
    consensus.fasta             this iteration's consensus
    stats.json                  reads mapped, length, identity to previous, consensus MD5,
                                 base composition, tool versions, per-command logs -- see below
    logs/                        captured stdout+stderr of every command run
  iter_001/
    <mapper>_index.*         index files, one set per configured mapper
    <mapper>.bam              that mapper's mapping output
    merged.bam                 (only if >1 mapper) merged BAM the consensus step sees
    consensus.fasta             this iteration's consensus
    stats.json                  reads mapped, length, identity to previous, consensus MD5,
                                 base composition, tool versions, per-command logs -- see below
    logs/                        captured stdout+stderr of every command run
  iter_002/
    ...
  metrics.tsv                 one row per iteration
  summary.json                 iterations run, converged?, total time, [output].commands log
  index.html                    human-readable report -- see "Logs" below
```

`reference_initial.fasta` is a *relative* symlink straight to your original
`reference_fasta` (or the NCBI-fetched cache file) whenever that's safe --
i.e. it already contains exactly the one sequence needed, nothing else --
rather than a copy, so a large reference genome doesn't get needlessly
duplicated. If `reference_id` had to pick one record out of a
multi-sequence `reference_fasta`, symlinking isn't possible (the file has
other sequences in it too), so that case still writes a real single-record
copy. Either way, `{reference}` behaves identically -- every tool that
reads it follows the symlink transparently. The BAM itself is never
symlinked here: getting from a full input BAM to what `iter_000` actually
needs (indexed, coordinate-sorted, and -- for `bam_reads` other than
`all` -- filtered down to reads for the chosen reference) is a real
transformation, not a copy, so `iteration_0_source.bam` under `reads/` is
always a genuine file.

### Consensus composition

`index.html` has a "Consensus composition" card for the final consensus:
the count of every character in it (the four unambiguous bases, plus any
IUPAC ambiguity codes or gap characters actually present), what percent of
the sequence is unambiguous (plain A/C/G/T), and GC content -- computed as
a fraction of just the unambiguous bases, since ambiguity codes and gaps
aren't G or C by definition and would otherwise just dilute the number.
This data lives in each iteration's `stats.json` under `"composition"` too
(every iteration's, not just the final one's), if you want it directly.

### Logs

`index.html` ends with a "Logs" section (only present if there's anything to
show -- a run with no `[tool-versions]` configured and no logged commands
skips it entirely):

- **Tool versions**, once each, if `[mapper.tool-versions]`/
  `[consensus.tool-versions]` are configured -- see "Tool versions" above for
  how they're collected, and how a version that changed mid-run is flagged.
- **One collapsible entry per iteration**, listing every command that
  iteration ran (`index_cmd`, `map_cmd` for each mapper, each `[consensus]`
  step) -- click to expand and see its full standard output/error and how
  long it took. A command that produced nothing just says so, rather than
  showing an empty block.
- **Final output**, if `[output].commands` ran, in the same collapsible
  format.

The same detail lives in each iteration's `stats.json` (`tool_versions`,
`commands`) and, for `[output].commands`, in `summary.json`
(`output_commands`) -- `index.html` is a rendering of that, not a separate
source of truth.

## Resuming a run

If a run stops because it hit `max_iterations` without converging, raise
`max_iterations` in the config and rerun the exact same command (same
`--config`, same `--output-dir`): it picks up from the next iteration rather
than starting over. This is detected automatically from `summary.json` in
`--output-dir` -- there's no separate flag. Already-extracted reads and
already-completed iterations aren't redone.

If a run already converged, rerunning it against the same `--output-dir` is a
no-op: it reports the existing result without redoing anything. The same is
true if it stopped because a cycle was detected (see "Cycle detection"
above) -- rerunning re-reports the same cycle rather than iterating further,
even against a `summary.json` written before cycle detection existed, since
that's checked retroactively too.

Resume trusts that `--output-dir` corresponds to the same logical run --
pointing it at a config with a different mapper, consensus pipeline, or
input isn't validated or rejected, it'll just continue on top of whatever's
there. Also, resuming is driven entirely by `summary.json`; a run that
crashed before writing it (e.g. killed mid-iteration) can't be resumed and
should be started fresh in a new `--output-dir`.

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

Optionally, `uv run pre-commit install` sets up a pre-commit hook that runs
the test suite (and keeps `uv.lock` in sync with `pyproject.toml`
automatically, regenerating and staging it if a commit -- e.g. a version
bump -- leaves it stale) before each commit. This is per-clone setup: it
writes into `.git/hooks/`, which isn't itself tracked by git, so it doesn't
happen automatically just from cloning the repo.
