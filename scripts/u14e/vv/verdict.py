"""The verdict ledger: every check lands one row; the exit code is the
verdict."""

class Verdict:
    """Every check lands one row; the exit code is the verdict.
    ``note`` rows are measured and recorded, never gating (the
    new-physics class)."""

    def __init__(self):
        self.rows = []

    def check(self, name, ok, detail):
        ok = bool(ok)
        self.rows.append((name, ok, "check"))
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}",
              flush=True)
        return ok

    def note(self, name, detail):
        self.rows.append((name, True, "note"))
        print(f"  [note] {name}: {detail}", flush=True)

    def report(self) -> bool:
        print("\n===== verdict =====", flush=True)
        for name, ok, kind in self.rows:
            tag = ("PASS" if ok else "FAIL") if kind == "check" \
                else "REPORT"
            print(f"  {tag}  {name}", flush=True)
        checks = [(n, ok) for n, ok, k in self.rows if k == "check"]
        good = sum(ok for _, ok in checks)
        print(f"{good}/{len(checks)} assertions hold", flush=True)
        return good == len(checks)
