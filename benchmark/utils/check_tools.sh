#!/bin/bash
#==============================================================================
# check_tools.sh - Verify all tools required by the benchmark suite
#==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

PASS=0
FAIL=0

pass() { echo "  [OK]  $1"; ((PASS++)); }
fail() { echo "  [!!]  $1"; ((FAIL++)); }

echo "=============================================================================="
echo "  freemap benchmark — dependency check"
echo "=============================================================================="
echo ""

# ── 1. freemap binary ────────────────────────────────────────────────────────

echo "=== freemap ==="
FREEMAP_PATH="$PROJECT_DIR/target/release/freemap"
CARGO_TOML="$PROJECT_DIR/Cargo.toml"

if [[ -x "$FREEMAP_PATH" ]]; then
    ver=$(grep '^version' "$CARGO_TOML" | head -1 | cut -d'"' -f2)
    pass "freemap $ver ($FREEMAP_PATH)"
else
    fail "freemap not found — run: cargo build --release"
fi
echo ""

# ── 2. samtools ──────────────────────────────────────────────────────────────

echo "=== samtools ==="
if command -v samtools &>/dev/null; then
    ver=$(samtools version 2>/dev/null | head -1 | awk '{print $2}')
    pass "samtools $ver ($(which samtools))"
else
    fail "samtools not found — needed for depth file generation"
fi
echo ""

# ── 3. Python + packages ────────────────────────────────────────────────────

echo "=== Python ==="
if command -v python3 &>/dev/null; then
    ver=$(python3 --version 2>&1 | awk '{print $2}')
    pass "python3 $ver ($(which python3))"
else
    fail "python3 not found"
fi

for pkg in numpy matplotlib scipy intervaltree; do
    if python3 -c "import $pkg" 2>/dev/null; then
        ver=$(python3 -c "import $pkg; print($pkg.__version__)" 2>/dev/null || echo "?")
        pass "$pkg $ver"
    else
        fail "$pkg not installed — run: pip install $pkg"
    fi
done
echo ""

# ── Summary ──────────────────────────────────────────────────────────────────

echo "=============================================================================="
echo "  Result: $PASS passed, $FAIL failed"
if [[ $FAIL -eq 0 ]]; then
    echo "  All dependencies satisfied. Ready to run benchmarks."
else
    echo "  Fix the above issues before running benchmarks."
fi
echo "=============================================================================="

exit $FAIL
