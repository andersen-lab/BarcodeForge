#!/usr/bin/env python3
"""Check that the package version is identical across every file that hardcodes it.

Sources checked:
  - pyproject.toml            (``[project]`` -> ``version = "..."``)
  - barcodeforge/__init__.py  (``__version__ = "..."``)

The script prints every version it finds and exits:
  - 0 if all sources agree,
  - 1 if they disagree or a version string cannot be found in a source.

It is dependency-free (regex based, standard library only) so it runs on any
Python >= 3.10 without installing the package or any third-party parser. To add
another source, append an entry to ``SOURCES`` below.

Run locally with:  python scripts/check_version_consistency.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent


def _extract_pyproject_version(text: str) -> str | None:
    """Return the ``version`` declared inside the ``[project]`` table."""
    # Isolate the [project] table so we never match a `version` from another
    # table (e.g. a future [tool.*] section).
    section = re.search(r"(?ms)^\[project\][^\[]*", text)
    if section is None:
        return None
    match = re.search(r"""(?m)^\s*version\s*=\s*["']([^"']+)["']""", section.group(0))
    return match.group(1) if match else None


def _extract_dunder_version(text: str) -> str | None:
    """Return the value of a ``__version__ = "..."`` assignment."""
    match = re.search(r"""(?m)^\s*__version__\s*=\s*["']([^"']+)["']""", text)
    return match.group(1) if match else None


# (label, path, extractor) — extend this list to check additional files.
SOURCES = [
    ("pyproject.toml", REPO_ROOT / "pyproject.toml", _extract_pyproject_version),
    (
        "barcodeforge/__init__.py",
        REPO_ROOT / "barcodeforge" / "__init__.py",
        _extract_dunder_version,
    ),
]


def main() -> int:
    versions: dict[str, str] = {}
    errors: list[str] = []

    for label, path, extractor in SOURCES:
        if not path.exists():
            errors.append(f"{label}: file not found ({path})")
            continue
        version = extractor(path.read_text(encoding="utf-8"))
        if version is None:
            errors.append(f"{label}: no version string found")
            continue
        versions[label] = version
        print(f"  {label}: {version}")

    if errors:
        print("\nERROR: could not read a version from every source:", file=sys.stderr)
        for err in errors:
            print(f"  - {err}", file=sys.stderr)
        return 1

    unique_versions = set(versions.values())
    if len(unique_versions) != 1:
        print("\nERROR: version mismatch across files:", file=sys.stderr)
        for label, version in versions.items():
            print(f"  - {label}: {version}", file=sys.stderr)
        return 1

    print(f"\nOK: all sources report version {unique_versions.pop()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
