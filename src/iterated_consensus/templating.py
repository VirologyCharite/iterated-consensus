"""Placeholder templating for mapper/consensus command specs.

Grammar for a placeholder ``{...}``:

    {name}              scalar substitution, or a space-joined file list
    {name:sep}          file list joined with the literal string ``sep``
    {prefix:name}       ``prefix`` immediately before each file, space-joined
    {prefix:name:sep}   ``prefix`` before each file, joined with ``sep``
    {cat:name}          path to a single file made by concatenating all of
                         ``name``'s files (created once and cached by the
                         caller-supplied ``cat_resolver``)

Whether the first colon-separated part is a ``prefix`` or the ``name`` itself
is decided by looking it up in ``values``: if it names a :class:`ReadsList`
already, there is no prefix. ``{{`` and ``}}`` are literal braces, matching
``str.format`` escaping.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path

from .errors import IteratedConsensusError


class TemplateError(IteratedConsensusError, ValueError):
    """Raised for malformed placeholders or invalid template values."""


@dataclass(frozen=True)
class ReadsList:
    """A named list of read files usable in a `{...}` placeholder."""

    name: str
    paths: tuple[Path, ...]

    def __post_init__(self) -> None:
        if not self.paths:
            raise TemplateError(f"reads list '{self.name}' has no files")


CatResolver = Callable[[ReadsList], Path]


def render(
    template: str,
    values: Mapping[str, object],
    *,
    cat_resolver: CatResolver | None = None,
) -> str:
    """Render a single template string, expanding all ``{...}`` placeholders."""
    out: list[str] = []
    i = 0
    n = len(template)
    while i < n:
        c = template[i]
        if c == "{":
            if template[i : i + 2] == "{{":
                out.append("{")
                i += 2
                continue
            j = template.find("}", i)
            if j == -1:
                raise TemplateError(f"unmatched '{{' in template: {template!r}")
            field = template[i + 1 : j]
            out.append(_render_field(field, values, cat_resolver))
            i = j + 1
            continue
        if c == "}":
            if template[i : i + 2] == "}}":
                out.append("}")
                i += 2
                continue
            raise TemplateError(f"unmatched '}}' in template: {template!r}")
        out.append(c)
        i += 1
    return "".join(out)


def _render_field(
    field: str,
    values: Mapping[str, object],
    cat_resolver: CatResolver | None,
) -> str:
    parts = field.split(":")

    prefix: str | None
    sep: str | None
    if len(parts) == 1:
        (name,) = parts
        prefix = sep = None
    elif len(parts) == 2:
        a, b = parts
        if isinstance(values.get(a), ReadsList):
            name, sep = a, b
            prefix = None
        else:
            prefix, name = a, b
            sep = None
    elif len(parts) == 3:
        prefix, name, sep = parts
    else:
        raise TemplateError(f"too many ':' in placeholder {{{field}}}")

    if name not in values:
        raise TemplateError(f"unknown placeholder '{{{field}}}': no value named '{name}'")
    value = values[name]

    if isinstance(value, ReadsList):
        return _render_reads_list(value, prefix, sep, cat_resolver)

    if prefix is not None or sep is not None:
        raise TemplateError(
            f"placeholder '{{{field}}}' uses prefix/separator modifiers, "
            f"but '{name}' is not a file list"
        )
    return str(value)


def _render_reads_list(
    value: ReadsList,
    prefix: str | None,
    sep: str | None,
    cat_resolver: CatResolver | None,
) -> str:
    if prefix == "cat":
        if sep is not None:
            raise TemplateError(
                f"'{{cat:{value.name}:{sep}}}' is invalid: "
                "'cat' cannot be combined with a separator"
            )
        if cat_resolver is None:
            raise TemplateError(
                f"'{{cat:{value.name}}}' used, but no cat_resolver was supplied"
            )
        return str(cat_resolver(value))

    item_prefix = "" if prefix is None else prefix
    item_sep = " " if sep is None else sep
    return item_sep.join(f"{item_prefix}{p}" for p in value.paths)
