#!/bin/bash
# config.sh - Shared configuration for freemap benchmark scripts
#
# Source this file from other scripts:
#   source "$(dirname "$0")/config.sh"
#
# Override DATA_DIR via environment variable or first argument to reproduce_paper.sh.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
DATA_DIR="${DATA_DIR:-$PROJECT_DIR/data}"
EVAL_DIR="$SCRIPT_DIR/evaluation"
FIG_DIR="$SCRIPT_DIR/figures"
FREEMAP="$PROJECT_DIR/target/release/freemap"
THREADS="${THREADS:-32}"

export DATA_DIR

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

log_info()  { echo -e "${GREEN}[INFO]${NC} $1"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }
