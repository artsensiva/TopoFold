#!/usr/bin/env bash
# ==============================================================================
# TopoFold: Automated Public Showcase Repository Generator
# ==============================================================================
# This script creates an isolated, public-facing portfolio repository in
# `../TopoFold-Showcase` that contains all publication documents, 3D web applications,
# benchmarks, figures, and pre-built stripped binary wheels, while completely
# omitting proprietary Rust source crates (`crates/topofold-*`).
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SHOWCASE_DIR="$(cd "${REPO_ROOT}/.." && pwd)/TopoFold-Showcase"

echo "================================================================================"
echo "🚀  TOPOFOLD: PUBLIC SHOWCASE REPOSITORY BUILDER"
echo "================================================================================"
echo "Source Repository:    ${REPO_ROOT}"
echo "Showcase Destination: ${SHOWCASE_DIR}"
echo ""

# 1. Ensure Pre-built Protected Binary Wheel Exists
echo "🔍 Step 1: Checking for pre-built stripped binary wheel..."
WHEEL_PATH=""
if compgen -G "${REPO_ROOT}/target/wheels/topofold-*.whl" > /dev/null; then
    WHEEL_PATH=$(ls -t "${REPO_ROOT}"/target/wheels/topofold-*.whl | head -n 1)
elif compgen -G "${REPO_ROOT}/dist/topofold-*.whl" > /dev/null; then
    WHEEL_PATH=$(ls -t "${REPO_ROOT}"/dist/topofold-*.whl | head -n 1)
else
    echo "⚙️  Wheel not found. Building protected wheel first..."
    "${SCRIPT_DIR}/build_protected_wheels.sh"
    WHEEL_PATH=$(ls -t "${REPO_ROOT}"/target/wheels/topofold-*.whl | head -n 1)
fi

echo "   ✓ Using binary wheel: ${WHEEL_PATH}"

# 2. Initialize Clean Showcase Directory
echo ""
echo "📂 Step 2: Preparing clean target directory: ${SHOWCASE_DIR}..."
rm -rf "${SHOWCASE_DIR}"
mkdir -p "${SHOWCASE_DIR}"
mkdir -p "${SHOWCASE_DIR}/wheels"
mkdir -p "${SHOWCASE_DIR}/apps"
mkdir -p "${SHOWCASE_DIR}/benchmarks"
mkdir -p "${SHOWCASE_DIR}/docs"
mkdir -p "${SHOWCASE_DIR}/assets"

# 3. Copy Safe Public Assets
echo ""
echo "📋 Step 3: Copying safe public assets (no proprietary Rust source code)..."

# Copy binary wheel
cp "${WHEEL_PATH}" "${SHOWCASE_DIR}/wheels/"
echo "   ✓ Copied pre-built stripped binary wheel to wheels/"

# Copy README and licenses
cp "${REPO_ROOT}/README.md" "${SHOWCASE_DIR}/"
cp "${REPO_ROOT}/LICENSE-APACHE" "${SHOWCASE_DIR}/" 2>/dev/null || true
cp "${REPO_ROOT}/LICENSE-MIT" "${SHOWCASE_DIR}/" 2>/dev/null || true
echo "   ✓ Copied README.md and licenses"

# Copy docs
cp -r "${REPO_ROOT}/docs/"* "${SHOWCASE_DIR}/docs/"
echo "   ✓ Copied docs/ (manuscript drafts, whitepapers, ADRs)"

# Copy assets
cp -r "${REPO_ROOT}/assets/"* "${SHOWCASE_DIR}/assets/"
echo "   ✓ Copied assets/ (all 300 DPI publication figures)"

# Copy apps
cp -r "${REPO_ROOT}/apps/"* "${SHOWCASE_DIR}/apps/"
echo "   ✓ Copied apps/ (Streamlit interactive 3D Web Dashboard)"

# Copy benchmarks and dataset files (excluding large zip archives)
find "${REPO_ROOT}/benchmarks" -maxdepth 1 -name "*.py" -exec cp {} "${SHOWCASE_DIR}/benchmarks/" \;
if [ -d "${REPO_ROOT}/benchmarks/data" ]; then
    mkdir -p "${SHOWCASE_DIR}/benchmarks/data"
    # Copy reference PDBs and DCDs, skip zip archives
    find "${REPO_ROOT}/benchmarks/data" -maxdepth 1 \( -name "*.pdb" -o -name "*.dcd" -o -name "*.npy" \) -exec cp {} "${SHOWCASE_DIR}/benchmarks/data/" \;
fi
echo "   ✓ Copied benchmarks/ and reference datasets"

# 4. Strict Intellectual Property Audit
echo ""
echo "🔒 Step 4: Intellectual property verification..."
CRATES_CHECK=$(find "${SHOWCASE_DIR}" -name "crates" -type d | wc -l)
RS_CHECK=$(find "${SHOWCASE_DIR}" -name "*.rs" -type f | wc -l)

if [ "${CRATES_CHECK}" -gt 0 ] || [ "${RS_CHECK}" -gt 0 ]; then
    echo "❌ SECURITY ALERT: Proprietary Rust source detected in showcase! Aborting."
    rm -rf "${SHOWCASE_DIR}"
    exit 1
fi
echo "   ✓ Verified: 0 Rust source files (.rs) present in showcase."
echo "   ✓ Verified: crates/ folder omitted."

# 5. Generate Showcase pyproject.toml & requirements.txt
echo ""
echo "⚙️  Step 5: Generating packaging and installation files..."
cat << 'EOF' > "${SHOWCASE_DIR}/pyproject.toml"
[build-system]
requires = ["setuptools>=61.0"]
build-backend = "setuptools.build_meta"

[project]
name = "topofold-showcase"
version = "0.8.3"
description = "Public Showcase, Interactive 3D Web Dashboard & Benchmarks for TopoFold"
readme = "README.md"
requires-python = ">=3.9"
dependencies = [
    "numpy>=1.20",
    "scipy>=1.8.0",
    "scikit-learn>=1.0.0",
    "streamlit>=1.28.0",
    "plotly>=5.15.0",
    "matplotlib>=3.5.0",
]
EOF

cat << 'EOF' > "${SHOWCASE_DIR}/requirements.txt"
numpy>=1.20
scipy>=1.8.0
scikit-learn>=1.0.0
streamlit>=1.28.0
plotly>=5.15.0
matplotlib>=3.5.0
EOF

cat << 'EOF' > "${SHOWCASE_DIR}/install.sh"
#!/usr/bin/env bash
set -e
echo "Installing TopoFold stripped native binary wheel..."
pip install wheels/topofold-*.whl
echo "Installing dashboard and benchmark dependencies..."
pip install -r requirements.txt
echo "Installation complete! Run dashboard with: streamlit run apps/streamlit_app.py"
EOF
chmod +x "${SHOWCASE_DIR}/install.sh"

cat << 'EOF' > "${SHOWCASE_DIR}/.gitignore"
__pycache__/
*.py[cod]
.venv/
env/
venv/
*.egg-info/
.DS_Store
.idea/
.vscode/
EOF

# 6. Initialize Fresh Git Repository in Showcase
echo ""
echo "🌱 Step 6: Initializing fresh Git repository in ${SHOWCASE_DIR}..."
cd "${SHOWCASE_DIR}"
git init -b main --quiet
git config user.name "Artem Galukhin" 2>/dev/null || true
git config user.email "artem.galukhin@gmail.com" 2>/dev/null || true
git add .
git commit -m "feat: initial public showcase repository for TopoFold v0.8.3" --quiet

echo ""
echo "================================================================================"
echo "🎉 SUCCESS: PUBLIC SHOWCASE REPOSITORY READY"
echo "================================================================================"
echo "Showcase Location: ${SHOWCASE_DIR}"
echo "Git Branch:        $(git branch --show-current)"
echo "Commit:            $(git rev-parse --short HEAD)"
echo "Total Files:       $(git ls-files | wc -l)"
echo ""
echo "To publish as a public GitHub repository:"
echo "  cd ${SHOWCASE_DIR}"
echo "  git remote add origin https://github.com/artsensiva/TopoFold-Showcase.git"
echo "  git push -u origin main"
echo "================================================================================"
