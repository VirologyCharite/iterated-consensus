import pytest

from iterated_consensus.config import parse_config
from iterated_consensus.presets import PresetError, get_preset_text, list_presets


def test_list_presets_nonempty() -> None:
    names = [name for name, _ in list_presets()]
    assert "bowtie2-ivar" in names
    assert "bwa-samtools" in names


@pytest.mark.parametrize("name", ["bowtie2-ivar", "bwa-samtools"])
def test_preset_text_is_valid_config_once_input_filled_in(name: str) -> None:
    text = get_preset_text(name)
    assert "[[mapper]]" in text
    assert "[consensus]" in text
    assert "[input]" in text
    # Presets ship an [input] header with its contents commented out; replace
    # it with a minimal real one so the bundled mapper/consensus commands
    # themselves parse and validate.
    text_with_input = text.replace(
        "[input]\n", '[input]\nreads_single = ["s.fq"]\nreference_fasta = "ref.fa"\n', 1
    )
    config = parse_config(text_with_input)
    assert config.mappers
    assert config.consensus.steps


def test_unknown_preset_raises() -> None:
    with pytest.raises(PresetError, match="unknown preset"):
        get_preset_text("does-not-exist")
