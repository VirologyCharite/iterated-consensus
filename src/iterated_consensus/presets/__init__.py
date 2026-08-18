"""Bundled example pipeline configs, for `iterated-consensus config-template`."""

from __future__ import annotations

from importlib import resources

from ..errors import IteratedConsensusError

_DESCRIPTIONS = {
    "bowtie2-ivar": "bowtie2 for mapping, ivar consensus for consensus calling",
    "bwa-samtools": "bwa mem for mapping, `samtools consensus` for consensus calling",
}
_FILENAMES = {
    "bowtie2-ivar": "bowtie2_ivar.toml",
    "bwa-samtools": "bwa_samtools.toml",
}


class PresetError(IteratedConsensusError, ValueError):
    """Raised when an unknown preset name is requested."""


def list_presets() -> list[tuple[str, str]]:
    """Return (name, description) pairs for every bundled preset."""
    return sorted(_DESCRIPTIONS.items())


def get_preset_text(name: str) -> str:
    if name not in _FILENAMES:
        available = ", ".join(sorted(_FILENAMES))
        raise PresetError(f"unknown preset '{name}'; available: {available}")
    return resources.files(__package__).joinpath(_FILENAMES[name]).read_text()
