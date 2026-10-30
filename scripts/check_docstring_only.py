"""Certify that a working-tree change is docstrings/comments only.

Compares the AST of each given path between HEAD (or --rev) and the
working tree, with every docstring node neutralized: identical dumps
mean no executable line changed. The docstring-pass certificate of
the module walk (AGENTS Verification rule 6); 2026-10-29: caught a
dropped function body that eyeballing the diff missed.

Usage:
    python scripts/check_docstring_only.py [--rev HEAD] path...
"""

from __future__ import annotations

import ast
import subprocess
import sys


def stripped_dump(src: str) -> str:
    """AST dump with every leading docstring REMOVED from its body --
    replaced, added, and deleted docstrings all vanish; any change to
    executable statements survives."""
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef,
                             ast.AsyncFunctionDef, ast.ClassDef)):
            body = node.body
            if (body and isinstance(body[0], ast.Expr)
                    and isinstance(body[0].value, ast.Constant)
                    and isinstance(body[0].value.value, str)):
                del node.body[0]
    return ast.dump(tree)


def main(argv: list[str]) -> int:
    """Compare each path's neutralized AST between the rev (default
    HEAD) and the working tree; exit 0 = every change is
    docstrings/comments only."""
    rev = "HEAD"
    paths = []
    it = iter(argv)
    for arg in it:
        if arg == "--rev":
            rev = next(it)
        else:
            paths.append(arg)
    ok = True
    for path in paths:
        old = subprocess.run(["git", "show", f"{rev}:{path}"],
                             capture_output=True, text=True).stdout
        new = open(path, encoding="utf-8").read()
        same = stripped_dump(old) == stripped_dump(new)
        ok &= same
        print(f"{path}: "
              f"{'AST-IDENTICAL (docstrings/comments only)' if same else 'AST DIFFERS -- CODE CHANGED'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
