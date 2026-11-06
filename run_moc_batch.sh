#!/bin/bash
# the moc batch (one serialized runner: concurrent regen invocations
# would race the store writes).  6 CASCI moc cells (the primary +
# the two simplest promoted windows, pol avg + pur) + 8 free optimizer
# moc cells (casscf + sa-casscf x pol/pur on the two promoted windows;
# the anchor-first gate needs --full to materialize them).
set -x
cd /home/kashpur/pylibgrpp/test/AtomTerm
export OMP_NUM_THREADS=64 OPENBLAS_NUM_THREADS=64 PYSCF_NUM_THREADS=64
P=~/miniconda3/envs/pyscf/bin/python
$P scripts/regen_cells.py U14e --start moc --method casci --window 10f13d
$P scripts/regen_cells.py U14e --start moc --method casci --window 10f13d7s
$P scripts/regen_cells.py U14e --start moc --method casci --window 9f10f11f12d13d
$P scripts/regen_cells.py U14e --full --mask - --start moc --method casscf --window 10f13d
$P scripts/regen_cells.py U14e --full --mask - --start moc --method sa-casscf --window 10f13d
$P scripts/regen_cells.py U14e --full --mask - --start moc --method casscf --window 10f13d7s
$P scripts/regen_cells.py U14e --full --mask - --start moc --method sa-casscf --window 10f13d7s
echo BATCH-DONE
