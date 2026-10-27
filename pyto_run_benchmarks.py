"""
One-tap launcher for Pyto (iPad/iPhone) or any Python IDE without a terminal.

Keep this file in the same folder as the other .py files of the repository
(ei_model.py, screening.py, dirac.py, ...), then press Run.
FAST = True skips the two slow screening-onset scans (15 checks instead of 17).
"""
import os
import sys

FAST = True

here = os.path.dirname(os.path.abspath(__file__))
os.chdir(here)
sys.path.insert(0, here)

import run_benchmarks
run_benchmarks.main(fast=FAST)
