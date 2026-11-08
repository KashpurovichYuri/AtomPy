#!/bin/sh
# Wait for the route-twin regen to exit, then complete the U14e table:
# warm the UHF anchor record into the store (the cached uhf.chk is
# trusted -- seconds), regenerate TABLE.md + TABLE.csv.  Launched
# detached 2026-11-07; the log is elements/U14e/logs/table_complete.log.
cd /home/kashpur/pylibgrpp/test/AtomTerm || exit 1
while pgrep -f "regen_cells.py U14e" > /dev/null; do sleep 60; done
echo "route batch exited at $(date)" >> elements/U14e/logs/table_complete.log
OMP_NUM_THREADS=8 ~/miniconda3/envs/pyscf/bin/python scripts/uhf_anchor.py U14e \
    >> elements/U14e/logs/table_complete.log 2>&1
OMP_NUM_THREADS=8 ~/miniconda3/envs/pyscf/bin/python scripts/build_methods_tables.py U14e \
    >> elements/U14e/logs/table_complete.log 2>&1
echo "completed at $(date)" >> elements/U14e/logs/table_complete.log
