"""How many tests were collected, for the README to be checked against.

Recorded here rather than measured inside a test: pytest already knows the number, and a
test that shelled out to `pytest --collect-only` would add three seconds to every run —
a hundred and eleven of them to the mutation survey, which runs the whole suite once per
mutation.
"""

from __future__ import annotations

#: Tests collected, when the whole suite was collected. ``None`` for a partial run —
#: a path argument or a ``-k`` selection, where the number means nothing.
COLLECTED: int | None = None


def pytest_collection_modifyitems(session, config, items):
    global COLLECTED
    selected = config.getoption("keyword") or config.getoption("markexpr")
    whole_suite = not selected and list(config.args) == list(config.getini("testpaths") or [])
    COLLECTED = len(items) if whole_suite else None
