"""A no-linter undefined-name check: in every module given, each free
name of a nested scope must be bound in an enclosing scope, at module
level, or be a builtin (name resolution that an import smoke test
cannot reach -- it sits in function bodies that only run later)."""

import builtins
import symtable
import sys

_MODULE_GLOBALS = {"__file__", "__name__", "__doc__", "__spec__",
                   "__package__", "__loader__", "__builtins__", "__path__"}


def _bound(table):
    """Names a scope binds: assigned, imported, parameters, declared."""
    return {s.get_name() for s in table.get_symbols()
            if s.is_assigned() or s.is_imported() or s.is_parameter()}


def check(path: str) -> list[str]:
    top = symtable.symtable(open(path).read(), path, "exec")
    module_names = _bound(top) | _MODULE_GLOBALS
    problems = [f"{path}: module scope references undefined {s.get_name()!r}"
                for s in top.get_symbols()
                if s.is_referenced() and s.get_name() not in module_names
                and not hasattr(builtins, s.get_name())]
    # (scope, the names bound by every scope enclosing it)
    stack = [(child, set()) for child in top.get_children()]
    while stack:
        t, enclosing = stack.pop()
        for s in t.get_symbols():
            name = s.get_name()
            if s.is_referenced() and not (s.is_assigned() or s.is_imported()
                                          or s.is_parameter()) \
                    and name not in enclosing and name not in module_names \
                    and not hasattr(builtins, name):
                problems.append(f"{path}: scope {t.get_name()!r} (line "
                                f"{t.get_lineno()}) references undefined {name!r}")
        stack.extend((child, enclosing | _bound(t)) for child in t.get_children())
    return problems


if __name__ == "__main__":
    bad = [p for f in sys.argv[1:] for p in check(f)]
    print("\n".join(bad) if bad else "no undefined names")
    sys.exit(1 if bad else 0)
