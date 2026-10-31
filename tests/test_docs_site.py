"""The documentation website: its pages exist, and the build hook rewrites what it should.

`mkdocs build --strict` in the docs workflow catches a broken site, but only on a
machine with MkDocs installed and only after the fact. The hook in `docs/mkdocs_hooks.py`
is plain Python, so what it promises (links resolved from the file that wrote them, the
CLI reference taken from the parser) is checked here with the rest of the suite.
"""

from __future__ import annotations

import importlib.util
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BLOB = "https://github.com/elmasnuryilmaz/splicescope/blob/main/"
RAW = "https://raw.githubusercontent.com/elmasnuryilmaz/splicescope/main/"


def _hooks():
    path = ROOT / "docs" / "mkdocs_hooks.py"
    spec = importlib.util.spec_from_file_location("mkdocs_hooks", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_every_page_in_the_navigation_exists_and_every_page_is_in_it():
    """A page missing from `nav` is built but unreachable; a `nav` entry with no file
    fails the strict build. Both are cheaper to find here."""
    config = (ROOT / "mkdocs.yml").read_text()
    nav = set(re.findall(r"^\s+- (?:[^:\n]+: )?([\w./-]+\.md)\s*$", config, re.M))
    pages = {p.name for p in (ROOT / "docs").glob("*.md")}
    assert nav, "could not read the navigation from mkdocs.yml"
    assert not nav - pages, f"listed but missing: {sorted(nav - pages)}"
    assert not pages - nav, f"pages the navigation does not list: {sorted(pages - nav)}"


def test_every_included_file_exists_and_maps_back_to_its_page():
    hooks = _hooks()
    for page in sorted((ROOT / "docs").glob("*.md")):
        for included in hooks._INCLUDE.findall(page.read_text()):
            assert (ROOT / included).is_file(), f"{page.name} includes missing {included}"
            assert hooks.INCLUDED_AS_PAGE.get(included) == page.name, (
                f"{included} is included by {page.name}, so links to it should lead there"
            )


def test_links_are_resolved_from_the_file_that_wrote_them():
    """The README is written at the repository root and shown at the site's root; the
    same link has to land in the same place from both."""
    hooks = _hooks()
    text = "\n".join(
        [
            "[m](docs/METHODS.md) [v](validation/README.md#dataset) [s](src/splicescope/io.py)",
            f"![fig]({RAW}docs/showcase.png) [abs]({BLOB}docs/REPRODUCING.md)",
            '<img src="docs/demo.gif"> [ext](https://example.org/x.md) [top](#install)',
            "```bash",
            "[not a link](docs/METHODS.md)",
            "```",
        ]
    )
    out = hooks._rewrite_links(text, "README.md", "index.md", "")
    assert "[m](METHODS.md)" in out
    assert "[v](validation.md#dataset)" in out
    assert f"[s]({BLOB}src/splicescope/io.py)" in out
    assert "![fig](showcase.png)" in out
    assert "[abs](REPRODUCING.md)" in out
    assert '<img src="demo.gif">' in out
    assert "[ext](https://example.org/x.md)" in out and "[top](#install)" in out
    assert "[not a link](docs/METHODS.md)" in out, "fenced code must be left alone"

    # Raw HTML is not rewritten by MkDocs, so it resolves from the page's URL.
    nested = hooks._rewrite_links(
        '<img src="../docs/demo.gif">', "validation/README.md", "validation.md", "validation/"
    )
    assert nested == '<img src="../demo.gif">'


def test_the_cli_reference_is_the_parser_s_own_help():
    hooks = _hooks()
    for command, option in (
        ("", "--version"),
        ("simulate", "--outdir"),
        ("run", "--sj-dir"),
        ("consequence", "--events"),
    ):
        rendered = hooks._cli_help(command)
        assert rendered.startswith("```text\nusage: splicescope"), rendered[:80]
        assert option in rendered, f"`splicescope {command} --help` lost {option}"
