"""The verdict ledger: every check lands one row; the exit code is the
verdict.  ``note`` rows are measured and recorded, never gating (the
new-physics class)."""


class Verdict:
    """A leg's assertion ledger.

    ``check`` gates; ``note`` discloses.  :meth:`report` returns
    whether every check holds -- and refuses an empty ledger: a leg
    with no assertion cannot be green (a notes-only leg dumps
    its ledger with :meth:`dump`).
    """

    def __init__(self):
        self.rows = []

    def check(self, name, ok, detail) -> bool:
        ok = bool(ok)
        self.rows.append((name, ok, "check"))
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}",
              flush=True)
        return ok

    def note(self, name, detail) -> None:
        self.rows.append((name, True, "note"))
        print(f"  [note] {name}: {detail}", flush=True)

    def adopt(self, rows) -> None:
        """Take another ledger's rows as they are -- a worker's, whose
        lines were printed with its case; nothing is printed here."""
        self.rows.extend(tuple(r) for r in rows)

    def _print(self, empty_ok: bool) -> list:
        checks = [(n, ok) for n, ok, kind in self.rows if kind == "check"]
        if not checks and not empty_ok:
            raise ValueError(
                "the verdict carries no assertion -- a leg gates on "
                "something (a notes-only ledger is dumped, not reported)")
        print("\n===== verdict =====", flush=True)
        for name, ok, kind in self.rows:
            tag = ("PASS" if ok else "FAIL") if kind == "check" \
                else "REPORT"
            print(f"  {tag}  {name}", flush=True)
        return checks

    def dump(self) -> None:
        """Print a notes-only ledger (no pass/fail claim)."""
        self._print(empty_ok=True)

    def report(self) -> bool:
        """Print the ledger; True iff every check holds.

        Raises
        ------
        ValueError : no check was registered.
        """
        checks = self._print(empty_ok=False)
        good = sum(ok for _, ok in checks)
        print(f"{good}/{len(checks)} assertions hold", flush=True)
        return good == len(checks)
