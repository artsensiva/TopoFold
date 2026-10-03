#!/usr/bin/env python3
"""
Compute SHA-256 hashes for M3 audit MANIFEST.
Temporary utility for correction pass.
"""
import hashlib
from pathlib import Path

def compute_hash(filepath):
    """Compute SHA-256 hash of a file."""
    sha256 = hashlib.sha256()
    with open(filepath, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            sha256.update(chunk)
    return sha256.hexdigest()

def main():
    base_dir = Path(__file__).parent.parent

    # Define artifact groups
    artifacts = {
        'STATUS': [
            'M3_STATUS.md',
            'README.md',
            'M3_A2G_CHECKPOINT_CLOSURE.md',
            'M3_A2G_COMPLETE_ANALYSIS.md',
            'M3_A2G_CORRECTIONS_PROVENANCE.md',
        ],
        'REPORTS': [
            'reports/M3_A2g_CHECKPOINT_REPORT.md',
            'reports/M3_A2f_CHECKPOINT_REPORT.md',
            'reports/M3_A2f_CLOSURE.md',
            'reports/M3_A2f_LIKELIHOOD_DEGENERACY.md',
            'reports/M3_A2e_CHECKPOINT_REPORT.md',
        ],
        'SCRIPTS_A2G': [
            'scripts/m3_a2g_convergence_audit.py',
            'scripts/m3_a2g_exact_k2_convergence.py',
            'scripts/m3_a2g_monotonicity_check.py',
            'scripts/m3_a2g_ninit_sensitivity.py',
            'scripts/m3_a2g_stopping_criterion_validation.py',
            'scripts/m3_a2g_extended_maxiter.py',
            'scripts/m3_a2g_duplicate_component_investigation.py',
            'scripts/m3_a2g_data_verification.py',
        ],
        'SCRIPTS_FITTERS': [
            'scripts/product_vm_fast_accurate.py',
            'scripts/product_vm_three_estimators.py',
        ],
        'RESULTS_A2G': [
            'results/m3_a2g_exact_family_convergence.csv',
            'results/m3_a2g_k2_convergence.csv',
            'results/m3_a2g_stress_regime.csv',
            'results/m3_a2g_extended_maxiter.csv',
            'results/m3_a2g_ninit_sensitivity.csv',
            'results/m3_a2g_stopping_criterion_validation.csv',
            'results/m3_a2g_duplicate_components.csv',
            'results/m3_a2g_initialization_comparison.csv',
            'results/m3_a2g_max_iter_sensitivity.csv',
            'results/m3_a2g_tolerance_comparison.csv',
            'results/m3_a2g_trajectories.txt',
            'results/m3_a2g_data_verification.json',
        ],
        'RESULTS_A2E': [
            'results/m3_a2e_checkpoint6_exact_family.csv',
            'results/m3_a2e_checkpoint9_k2_recovery.csv',
            'results/m3_a2e_checkpoint10_wrapped_normal.csv',
        ],
        'RESULTS_EARLIER': [
            'results/m3_a2c_exact_null_results.csv',
        ],
    }

    print("=" * 80)
    print("M3 MANIFEST HASH COMPUTATION")
    print("=" * 80)

    for group, files in artifacts.items():
        print(f"\n{group}:")
        for rel_path in files:
            filepath = base_dir / rel_path
            if filepath.exists():
                file_hash = compute_hash(filepath)
                print(f"  {file_hash}  {rel_path}")
            else:
                print(f"  MISSING: {rel_path}")

if __name__ == "__main__":
    main()
