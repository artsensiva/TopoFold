"""Held-out decision statistics with the Amendment-2 candidate set (docs/WPD_HELDOUT_MDCATH_AMENDMENT_2.md §1).

Loads the frozen Amendment-1 implementation (heldout_decision.py, unmodified) as a separate module instance and
changes only the best-method candidate set to include A. Estimand, intervals, margins and classes are those of
Amendment 1. Harness code, not production.
"""
import importlib.util
from pathlib import Path

_spec = importlib.util.spec_from_file_location('heldout_decision_a2_impl',
                                               Path(__file__).resolve().parent / 'heldout_decision.py')
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
_mod.BEST_ORDER = ('B2', 'B1', 'B3', 'M1', 'M2', 'A')      # A appended; order breaks exact ties only

BEST_ORDER = _mod.BEST_ORDER
select_best = _mod.select_best
cell_statistics = _mod.cell_statistics
bootstrap = _mod.bootstrap
analyse_cell = _mod.analyse_cell
overall_outcome = _mod.overall_outcome
check_inputs = _mod.check_inputs
conservative_interval = _mod.conservative_interval
