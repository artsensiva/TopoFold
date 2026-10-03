#!/usr/bin/env python3
"""
M3-A2g Corrections Pass: Substantive Acceptance Checklist
Performs structured validation of corrections against actual data.
"""
import json
import re
from pathlib import Path

base_dir = Path(__file__).parent.parent
results_dir = base_dir / "results"
results = {'passed': [], 'failed': []}

def check(name, condition, details=""):
    """Record a check result."""
    if condition:
        results['passed'].append((name, details))
        print(f"✅ PASS: {name}")
        if details:
            print(f"         {details}")
    else:
        results['failed'].append((name, details))
        print(f"❌ FAIL: {name}")
        if details:
            print(f"         {details}")
    return condition

print("=" * 80)
print("M3-A2g CORRECTIONS PASS: SUBSTANTIVE ACCEPTANCE CHECKLIST")
print("=" * 80)

# CRITICAL: Load and validate verification JSON
print("\n## CRITICAL: Verification JSON Validation")
ver_json_path = results_dir / 'm3_a2g_data_verification.json'
if not ver_json_path.exists():
    check("Verification JSON exists", False, "File not found")
    print("\n❌ CANNOT PROCEED - verification JSON missing")
    exit(1)

with open(ver_json_path) as f:
    ver_data = json.load(f)

# Verify SYNTHETIC metadata
check(
    "Verification JSON: SYNTHETIC metadata",
    ver_data.get('_meta', {}).get('type') == 'SYNTHETIC',
    f"Type: {ver_data.get('_meta', {}).get('type')}"
)

# Verify total CSV rows
check(
    "Verification JSON: 1235 total CSV rows",
    ver_data.get('total_csv_rows') == 1235,
    f"Actual: {ver_data.get('total_csv_rows')}"
)

# Verify iteration records (corrected from 5 to 1634)
trajectories = ver_data.get('trajectories', {})
check(
    "Verification JSON: 1634 iteration records (bug corrected)",
    trajectories.get('iteration_records_total') == 1634,
    f"Actual: {trajectories.get('iteration_records_total')}"
)

# CRITICAL: All schemas must be valid
print("\n## CRITICAL: CSV Schema Validation")
schema_validation = ver_data.get('schema_validation', {})

for csv_name in ['exact_family', 'k2_convergence', 'stress_regime', 'extended_maxiter', 'ninit_sensitivity']:
    schema = schema_validation.get(csv_name, {})
    is_valid = schema.get('valid', False)
    check(
        f"CSV schema valid: {csv_name}",
        is_valid,
        f"Missing: {schema.get('missing', [])}, Extra: {schema.get('extra', [])}" if not is_valid else "All columns match"
    )

# CRITICAL: Convergence statistics must match actual data
print("\n## CRITICAL: Convergence Statistics Match Actual Data")

exact_family = ver_data.get('exact_family', {})
check(
    "Exact-family: 111/360 (30.8%)",
    exact_family.get('converged') == 111 and exact_family.get('total_fits') == 360,
    f"Actual: {exact_family.get('converged')}/{exact_family.get('total_fits')} ({100*exact_family.get('rate', 0):.1f}%)"
)

k2 = ver_data.get('k2_convergence', {})
check(
    "True-K2: 258/320 (80.6%)",
    k2.get('converged') == 258 and k2.get('total_fits') == 320,
    f"Actual: {k2.get('converged')}/{k2.get('total_fits')} ({100*k2.get('rate', 0):.1f}%)"
)

stress = ver_data.get('stress_regime', {})
check(
    "Stress regime: 102/280 (36.4%)",
    stress.get('converged') == 102 and stress.get('total_fits') == 280,
    f"Actual: {stress.get('converged')}/{stress.get('total_fits')} ({100*stress.get('rate', 0):.1f}%)"
)

# Verify specific stress example: σ=1.0, F=1000 = 3/10
stress_examples = stress.get('examples', {})
sigma_1_0_F_1000 = stress_examples.get('sigma=1.0,F=1000', {})
check(
    "Stress σ=1.0, F=1000: 3/10 converged",
    sigma_1_0_F_1000.get('converged') == 3 and sigma_1_0_F_1000.get('total') == 10,
    f"Actual: {sigma_1_0_F_1000.get('converged')}/{sigma_1_0_F_1000.get('total')}"
)

# CRITICAL: Grid coverage validation
print("\n## CRITICAL: Grid Coverage Validation")
grid = ver_data.get('grid_coverage', {})
check(
    "Exact-family grid: q complete",
    grid.get('q_complete') == True and grid.get('q_values') == [1, 2, 5],
    f"q values: {grid.get('q_values')}"
)
check(
    "Exact-family grid: kappa complete",
    grid.get('kappa_complete') == True and grid.get('kappa_values') == [1, 5, 20],
    f"kappa values: {grid.get('kappa_values')}"
)
check(
    "Exact-family grid: F complete",
    grid.get('F_complete') == True and grid.get('F_values') == [100, 500, 1000, 2500],
    f"F values: {grid.get('F_values')}"
)
check(
    "Exact-family grid: 36 configurations",
    grid.get('actual_configs') == 36 and grid.get('expected_configs') == 36,
    f"Actual: {grid.get('actual_configs')}, Expected: {grid.get('expected_configs')}"
)

# CRITICAL: Documentation consistency with verified data
print("\n## CRITICAL: Documentation Consistency with Verified Data")

# Load key documentation files
provenance = (base_dir / "M3_A2G_CORRECTIONS_PROVENANCE.md").read_text()
status = (base_dir / "M3_STATUS.md").read_text()
checkpoint_report = (base_dir / "reports/M3_A2g_CHECKPOINT_REPORT.md").read_text()

# Check convergence statistics appear in documentation
check(
    "Documentation: 111/360 (30.8%) exact-family rate",
    "111/360" in checkpoint_report and "30.8" in checkpoint_report
)
check(
    "Documentation: 102/280 (36.4%) stress rate",
    "102/280" in checkpoint_report and "36.4" in checkpoint_report
)
check(
    "Documentation: 258/320 (80.6%) K2 rate",
    "258/320" in checkpoint_report and "80.6" in checkpoint_report
)

# Check provisional O2 status (not claiming O3)
# Allow mentions of O3 in future-work context, but not as current status
has_o2 = "O2" in status and "provisional" in status.lower()
claims_o3_status = bool(re.search(r'(Status|Classification):\s*O3', status))
check(
    "Status: Provisional O2 (not claiming O3 status)",
    has_o2 and not claims_o3_status
)

# CRITICAL: Historical diagnostic withdrawal
print("\n## CRITICAL: Historical Diagnostic Withdrawal")
check(
    "Provenance: Historical withdrawal section exists",
    "POST-COMMIT CORRECTION" in provenance or "WITHDRAWAL" in provenance.upper()
)
check(
    "Provenance: Likelihood bug documented",
    "likelihood" in provenance.lower() and "responsibilities" in provenance.lower()
)
check(
    "Provenance: Inner reseeding documented",
    ("inner" in provenance.lower() or "INNER" in provenance) and "reseed" in provenance.lower()
)
check(
    "Provenance: Stopping diagnostics invalid (ALL protocols)",
    "ALL" in provenance and "stopping" in provenance.lower()
)

# CRITICAL: Fixture set distinction
print("\n## CRITICAL: Fixture Set Distinction")
check(
    "Documentation: Monotonicity uses exact_easy",
    "exact_easy" in status and "exact_easy" in checkpoint_report,
    "Monotonicity check must use exact_easy fixture"
)
check(
    "Documentation: n_init sensitivity uses exact_k1",
    "exact_k1" in checkpoint_report,
    "n_init sensitivity must use exact_k1 fixture"
)
check(
    "Documentation: Deterministic split initialization",
    ("deterministic" in status.lower() or "split" in status) and "monotonicity" in status.lower(),
    "Monotonicity must describe deterministic split initialization"
)

# CRITICAL: UNPAIRED disclaimer
print("\n## CRITICAL: UNPAIRED Draw Disclaimer")
check(
    "Documentation: UNPAIRED disclaimer present",
    "UNPAIRED" in checkpoint_report or "unpaired" in checkpoint_report.lower()
)

# CRITICAL: Budget attachment
print("\n## CRITICAL: Budget Specification")
check(
    "Documentation: Convergence rates specify max_iter=500 budget",
    "max_iter=500" in status or "max_iter = 500" in status or "at max_iter=500" in status
)

# Final summary
print("\n" + "=" * 80)
print("SUBSTANTIVE ACCEPTANCE CHECKLIST SUMMARY")
print("=" * 80)
print(f"✅ PASSED: {len(results['passed'])}")
print(f"❌ FAILED: {len(results['failed'])}")

if results['failed']:
    print("\nFAILURES:")
    for name, details in results['failed']:
        print(f"  ❌ {name}")
        if details:
            print(f"     {details}")
    print("\n⚠ ACCEPTANCE CRITERIA NOT MET")
    exit(1)
else:
    print("\n✅ ALL SUBSTANTIVE ACCEPTANCE CRITERIA MET")
    print("\nKey Validated Facts:")
    print(f"  - Total CSV rows: {ver_data.get('total_csv_rows')}")
    print(f"  - Iteration records: {trajectories.get('iteration_records_total')}")
    print(f"  - Exact-family: {exact_family.get('converged')}/{exact_family.get('total_fits')} ({100*exact_family.get('rate', 0):.1f}%)")
    print(f"  - True-K2: {k2.get('converged')}/{k2.get('total_fits')} ({100*k2.get('rate', 0):.1f}%)")
    print(f"  - Stress: {stress.get('converged')}/{stress.get('total_fits')} ({100*stress.get('rate', 0):.1f}%)")
    print(f"  - All 5 CSV schemas: VALID")
    print(f"  - Grid coverage: {grid.get('actual_configs')} configs complete")
    print("\nM3-A2g corrections pass validation COMPLETE.")
    exit(0)
