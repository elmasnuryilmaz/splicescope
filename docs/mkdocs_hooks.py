"""MkDocs build hooks for the documentation site (`mkdocs.yml` loads this file).

The site is built from files that also have to read correctly on GitHub: the README,
`validation/README.md`, `CONTRIBUTING.md`, the changelog and `docs/*.md`. Copying them
into `docs/` would give every number two places to drift apart in, so pages include
them instead, and this hook does three things at build time:

``<!-- include: PATH -->``
    Replaced by the contents of PATH (relative to the repository root).

``<!-- cli-help: COMMAND -->``
    Replaced by the live ``--help`` text of ``splicescope COMMAND`` (or of the top-level
    parser when COMMAND is empty), so the CLI reference cannot fall behind the parser.

Relative links
    Every relative link and image, and every absolute link into this repository, in a
    page or in a file it includes, is resolved from the file it was written in. A target
    inside ``docs/`` becomes a path relative to the page; a file the site includes
    becomes a link to that page; anything else in the repository becomes a link to it on
    GitHub. Fenced code blocks are left untouched.

Only the standard library and this package are imported, so the hook adds no
dependency beyond MkDocs itself.
"""

from __future__ import annotations

import contextlib
import io
import os
import posixpath
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REPO_URL = "https://github.com/elmasnuryilmaz/splicescope"
BRANCH = "main"
#: Absolute URLs into this repository. The README uses them so that its links and images
#: also work on PyPI; on the site they are mapped back to pages and local files.
_OWN_URLS = (
    f"{REPO_URL}/blob/{BRANCH}/",
    f"https://raw.githubusercontent.com/elmasnuryilmaz/splicescope/{BRANCH}/",
)

#: Repository files the site shows as a page of its own, and that page.
INCLUDED_AS_PAGE = {
    "README.md": "index.md",
    "validation/README.md": "validation.md",
    "CONTRIBUTING.md": "contributing.md",
    "CODE_OF_CONDUCT.md": "code_of_conduct.md",
    "SUPPORT.md": "support.md",
    "CHANGELOG.md": "changelog.md",
}

_INCLUDE = re.compile(r"^<!-- include: (?P<path>\S+) -->[ \t]*$", re.M)
_CLI_HELP = re.compile(r"^<!-- cli-help:(?P<command>[ \w-]*) -->[ \t]*$", re.M)
_FENCE = re.compile(r"^(```|~~~)", re.M)
_MD_LINK = re.compile(r"(!?\[[^\]\n]*\]\()(?P<target>[^)\s]+)(\s+\"[^\"]*\")?\)")
_HTML_ATTR = re.compile(r"(?P<attr>\b(?:src|href))=\"(?P<target>[^\"]+)\"")


def _is_relative(target: str) -> bool:
    if target.startswith(("#", "/", "mailto:", "data:")):
        return False
    return not re.match(r"^[a-z][a-z0-9+.-]*:", target)


def _resolve(target: str, origin: str, page: str, page_url: str, html: bool) -> str:
    """Rewrite one relative ``target`` written in repository file ``origin``."""
    for prefix in _OWN_URLS:
        if target.startswith(prefix):
            local = _resolve(target[len(prefix) :], "README.md", page, page_url, html)
            return target if local.startswith(REPO_URL) else local
    path, _, anchor = target.partition("#")
    anchor = f"#{anchor}" if anchor else ""
    if not path:
        return target
    repo_path = posixpath.normpath(posixpath.join(posixpath.dirname(origin), path))
    if repo_path in INCLUDED_AS_PAGE:
        repo_path = f"docs/{INCLUDED_AS_PAGE[repo_path]}"
    if repo_path.startswith("docs/") and (ROOT / repo_path).exists():
        in_docs = repo_path[len("docs/") :]
        if html:
            # MkDocs rewrites Markdown links, not raw HTML: resolve from the page's URL.
            start = page_url
            if page_url.endswith(".html"):
                start = posixpath.dirname(page_url)
            return posixpath.relpath(in_docs, start or ".") + anchor
        return posixpath.relpath(in_docs, posixpath.dirname(page) or ".") + anchor
    if (ROOT / repo_path).exists():
        kind = "tree" if (ROOT / repo_path).is_dir() else "blob"
        return f"{REPO_URL}/{kind}/{BRANCH}/{repo_path}{anchor}"
    return target  # unknown: leave it for MkDocs to report


def _rewrite_links(text: str, origin: str, page: str, page_url: str) -> str:
    """Rewrite relative links outside fenced code blocks."""
    parts = re.split(r"(^(?:```|~~~).*?^(?:```|~~~)[^\n]*$)", text, flags=re.M | re.S)
    for i, part in enumerate(parts):
        if _FENCE.match(part):
            continue

        def md(m: re.Match) -> str:
            target = m.group("target")
            if not (_is_relative(target) or target.startswith(_OWN_URLS)):
                return m.group(0)
            new = _resolve(target, origin, page, page_url, html=False)
            return f"{m.group(1)}{new}{m.group(3) or ''})"

        def attr(m: re.Match) -> str:
            target = m.group("target")
            if not (_is_relative(target) or target.startswith(_OWN_URLS)):
                return m.group(0)
            return f'{m.group("attr")}="{_resolve(target, origin, page, page_url, html=True)}"'

        parts[i] = _HTML_ATTR.sub(attr, _MD_LINK.sub(md, part))
    return "".join(parts)


def _cli_help(command: str) -> str:
    from splicescope.cli import main

    argv = [*command.split(), "--help"]
    out = io.StringIO()
    previous = os.environ.get("COLUMNS")
    os.environ["COLUMNS"] = "88"  # argparse wraps to the terminal; pin it
    try:
        with contextlib.redirect_stdout(out), contextlib.suppress(SystemExit):
            main(argv)
    finally:
        if previous is None:
            del os.environ["COLUMNS"]
        else:
            os.environ["COLUMNS"] = previous
    return f"```text\n{out.getvalue().rstrip()}\n```"


def on_page_markdown(markdown: str, page, config, files) -> str:
    page_src = page.file.src_uri
    page_url = page.url
    origin = f"docs/{page_src}"
    markdown = _CLI_HELP.sub(lambda m: _cli_help(m.group("command").strip()), markdown)

    pieces, last = [], 0
    for m in _INCLUDE.finditer(markdown):
        pieces.append(_rewrite_links(markdown[last : m.start()], origin, page_src, page_url))
        included = m.group("path")
        text = (ROOT / included).read_text(encoding="utf-8")
        pieces.append(_rewrite_links(text, included, page_src, page_url))
        last = m.end()
    pieces.append(_rewrite_links(markdown[last:], origin, page_src, page_url))
    return "".join(pieces)
