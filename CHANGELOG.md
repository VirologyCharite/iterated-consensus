# 0.2.6 August 21, 2026

Add symlink to `final-reference.fasta.fai` if the `.fai` file exists.

# 0.2.5 August 21, 2026

Allow `output_dir` in the config. Add `[tool-versions]` so the summary
output can indicate version info. Make all config `[run]` variables
available as placeholders in other parts of the config. Added
iteration cycle detection. Improved `index.html` summary. Added
`final_reference_fasta` and `final_reference_bam` to the config.

# 0.2.4 August 18, 2026

Made the `--progress` output have a header line instead of repeating
the column names on each output line. The final consensus filename is
now relative to the output directory (unless it is given as an
absolute path).

# 0.2.3 August 18, 2026

Added config [output] section. Symlink iteration zero files when
possible. Don't create anything when -n is given. Turn [run] variables
into {name} placeholders that can be used elsewhere. Added better
threads support (including 'auto').

# 0.2.1 August 18, 2026

Added automatic sorting and indexing of BAM files if the mapper
commands happen not to do it.

# 0.2.0 August 18, 2026

Simplified `reference_id` and `reference_fasta` config variables.

# 0.1.2 August 18, 2026

Updated `pyproject.toml`.

# 0.1.1 August 18, 2026

`reference_name` and `reference_id` unified into a single
`reference_id` field/flag, used for both a FASTA record ID start (from
FASTQ) or a BAM start (using a reference from the pre-existing BAM
file).

# 0.1.0 August 18, 2026

Initial release. Code written entirely by Claude Sonnet 5 (model ID
`claude-sonnet-5`), from Anthropic's Claude 5 family.
