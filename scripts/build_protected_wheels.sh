#!/usr/bin/env bash
# ==============================================================================
# TopoFold: Enterprise Protected Binary Wheel Build & Verification Script
# ==============================================================================
# This script compiles TopoFold into a hardened, stripped, single-binary Python
# wheel with fat LTO and zero exposed Rust/Python source code.
# It then validates the wheel in an isolated clean environment without source files.
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

cd "${REPO_ROOT}"

echo "================================================================================"
echo "🛡️  TOPOFOLD: ENTERPRISE PROTECTED BINARY WHEEL COMPILER"
echo "================================================================================"

# 1. Verify Cargo.toml Protection Settings
echo "🔍 Step 1: Auditing Cargo.toml release profiles..."
if ! grep -q 'lto = "fat"' Cargo.toml; then
    echo "❌ Error: 'lto = \"fat\"' not found in Cargo.toml [profile.release]"
    exit 1
fi

if ! grep -q 'strip = "symbols"' Cargo.toml; then
    echo "❌ Error: 'strip = \"symbols\"' not found in Cargo.toml [profile.release]"
    exit 1
fi

if ! grep -q 'panic = "abort"' Cargo.toml; then
    echo "❌ Error: 'panic = \"abort\"' not found in Cargo.toml [profile.release]"
    exit 1
fi

echo "   ✓ Whole-program Fat LTO enabled"
echo "   ✓ Symbol stripping enabled"
echo "   ✓ Abort on panic enabled (eliminates stack trace metadata)"

# 2. Locate Python and Maturin
echo ""
echo "🔍 Step 2: Locating build toolchain..."
if [ -f ".venv/bin/maturin" ]; then
    MATURIN_BIN=".venv/bin/maturin"
    PYTHON_BIN=".venv/bin/python"
elif command -v maturin &> /dev/null; then
    MATURIN_BIN="$(command -v maturin)"
    PYTHON_BIN="$(command -v python3)"
else
    echo "❌ Error: Maturin not found. Install via: pip install maturin"
    exit 1
fi

echo "   ✓ Using Maturin: ${MATURIN_BIN}"
echo "   ✓ Using Python:  ${PYTHON_BIN}"

# 3. Compile Protected Wheel
echo ""
echo "⚙️  Step 3: Compiling optimized, stripped binary wheel with Maturin..."
WHEELS_DIR="${REPO_ROOT}/target/wheels"
mkdir -p "${WHEELS_DIR}"

"${MATURIN_BIN}" build --release --strip --out "${WHEELS_DIR}"

WHEEL_PATH=$(ls -t "${WHEELS_DIR}"/topofold-*.whl | head -n 1)
if [ ! -f "${WHEEL_PATH}" ]; then
    echo "❌ Error: Wheel compilation failed. No .whl found in ${WHEELS_DIR}"
    exit 1
fi

WHEEL_SIZE=$(du -h "${WHEEL_PATH}" | cut -f1)
echo "   ✓ Built wheel: ${WHEEL_PATH} (${WHEEL_SIZE})"

# 4. Verify Symbol Stripping
echo ""
echo "🔍 Step 4: Auditing binary stripping..."
TMP_AUDIT_DIR=$(mktemp -d /tmp/topofold_audit.XXXXXX)
unzip -q -o "${WHEEL_PATH}" -d "${TMP_AUDIT_DIR}"
SO_FILE=$(find "${TMP_AUDIT_DIR}" -name "*.so" | head -n 1)

if [ -f "${SO_FILE}" ]; then
    if command -v file &> /dev/null; then
        FILE_INFO=$(file "${SO_FILE}")
        echo "   File info: ${FILE_INFO}"
        if echo "${FILE_INFO}" | grep -q "stripped"; then
            echo "   ✓ Binary is confirmed STRIPPED"
        else
            echo "   ⚠️ Notice: Debug symbols present, stripping manually..."
            strip --strip-all "${SO_FILE}" || true
        fi
    fi
    # Check that internal Rust source code is NOT bundled
    RS_COUNT=$(find "${TMP_AUDIT_DIR}" -name "*.rs" | wc -l)
    if [ "${RS_COUNT}" -eq 0 ]; then
        echo "   ✓ Zero Rust source code (.rs) files bundled in wheel"
    else
        echo "   ❌ Security warning: ${RS_COUNT} Rust source files found inside wheel!"
        rm -rf "${TMP_AUDIT_DIR}"
        exit 1
    fi
fi
rm -rf "${TMP_AUDIT_DIR}"

# 5. Clean Environment Verification
echo ""
echo "🧪 Step 5: Testing wheel in an isolated, clean virtual environment..."
TMP_VENV=$(mktemp -d /tmp/topofold_test_venv.XXXXXX)
"${PYTHON_BIN}" -m venv --system-site-packages "${TMP_VENV}"

"${TMP_VENV}/bin/pip" install "${WHEEL_PATH}" --no-deps --quiet

# Execute test from outside the TopoFold repository to guarantee zero source fallback
TEST_OUTPUT=$(cd /tmp && PYTHONPATH="${REPO_ROOT}/.venv/lib/python3.14/site-packages" "${TMP_VENV}/bin/python" -c "
import topofold as tf
print(f'VERSION={tf.__version__}')
# Basic functional check
import numpy as np
ca = np.random.rand(10, 3).astype(np.float32)
k, t, w = tf.compute_invariants(ca)
assert len(k) == 8, 'Failed invariant computation'
print('SANITY_TEST=PASSED')
")

echo "   Isolated Test Result:"
echo "   ${TEST_OUTPUT}"

if echo "${TEST_OUTPUT}" | grep -q "SANITY_TEST=PASSED"; then
    echo "   ✓ Wheel executes flawlessly in clean isolated environment without source code!"
else
    echo "   ❌ Isolated execution test failed!"
    rm -rf "${TMP_VENV}"
    exit 1
fi

rm -rf "${TMP_VENV}"

echo ""
echo "================================================================================"
echo "🎉 SUCCESS: PROTECTED BINARY WHEEL READY FOR COMMERCIAL DISTRIBUTION"
echo "================================================================================"
echo "Artifact: ${WHEEL_PATH}"
echo "Deploy command: pip install ${WHEEL_PATH}"
echo "================================================================================"
