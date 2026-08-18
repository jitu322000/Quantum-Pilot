#!/bin/bash
# install.sh
#
# Installs gaussbot and/or gamessbot from wherever THIS script
# actually lives -- no hardcoded paths, works no matter where the BOT
# directory was copied to or which user runs it. Each bot gets its
# own venv (matching how they were already set up), plus a portable
# run_gui.sh wrapper and a per-bot .env file holding whatever
# machine-specific settings (Gaussian environment, GAMESS rungms/
# scratch dir, ports) this run collected -- nothing gets written into
# your shell's own ~/.bashrc.
#
# Usage:
#   ./install.sh                 # interactive
#   BOT_INSTALL_CHOICE=3 ./install.sh --yes   # non-interactive, both bots,
#                                               # auto-detected/default settings
#
# Non-interactive overrides (skip prompts when set):
#   BOT_INSTALL_CHOICE   1=gaussbot only, 2=gamessbot only, 3=both (default)
#   BOT_G09ROOT           Gaussian install root, e.g. /home/alice
#   BOT_GAUSS_SCRDIR      Gaussian scratch dir, e.g. /scratch
#   BOT_RUNGMS_PATH       path to GAMESS's rungms script
#   BOT_GAMESS_SCRATCH_DIR  GAMESS scratch dir

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GAUSSBOT_DIR="$SCRIPT_DIR/gaussbot"
GAMESSBOT_DIR="$SCRIPT_DIR/gamessbot"

INTERACTIVE=1
[ -t 0 ] || INTERACTIVE=0
for arg in "$@"; do
  [ "$arg" = "--yes" ] || [ "$arg" = "-y" ] && INTERACTIVE=0
done

echo "=== BOT installer ==="
echo "Installing from: $SCRIPT_DIR"
echo

PYTHON_BIN="${PYTHON_BIN:-python3}"
if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "python3 not found on PATH -- install Python 3.9+ first." >&2
  exit 1
fi

# ------------------------------------------------------- which bot(s)

CHOICE="${BOT_INSTALL_CHOICE:-}"
if [ -z "$CHOICE" ] && [ "$INTERACTIVE" = "1" ]; then
  echo "Which would you like to install?"
  echo "  1) gaussbot only"
  echo "  2) gamessbot only"
  echo "  3) both (default)"
  read -rp "Choice [3]: " CHOICE
fi
CHOICE="${CHOICE:-3}"

INSTALL_GAUSSBOT=0
INSTALL_GAMESSBOT=0
case "$CHOICE" in
  1) INSTALL_GAUSSBOT=1 ;;
  2) INSTALL_GAMESSBOT=1 ;;
  *) INSTALL_GAUSSBOT=1; INSTALL_GAMESSBOT=1 ;;
esac

# ------------------------------------------- Gaussian environment (g09root/GAUSS_SCRDIR)
# Needed by gaussbot itself, and by gamessbot's optional "guess
# geometry, optimize with Gaussian first" intake path.

FOUND_G09ROOT=""
FOUND_GAUSS_SCRDIR=""

detect_or_prompt_gaussian_env() {
  if [ -n "$BOT_G09ROOT" ]; then
    FOUND_G09ROOT="$BOT_G09ROOT"
    FOUND_GAUSS_SCRDIR="${BOT_GAUSS_SCRDIR:-/scratch}"
    echo "Using Gaussian environment from BOT_G09ROOT/BOT_GAUSS_SCRDIR: g09root=$FOUND_G09ROOT GAUSS_SCRDIR=$FOUND_GAUSS_SCRDIR"
    return
  fi

  if [ -n "$g09root" ] && [ -n "$GAUSS_SCRDIR" ]; then
    echo "Using Gaussian environment already set in this shell (g09root=$g09root)."
    FOUND_G09ROOT="$g09root"
    FOUND_GAUSS_SCRDIR="$GAUSS_SCRDIR"
    return
  fi

  local rc root scr
  for rc in "$HOME/.bashrc" "$HOME/.bash_profile" "$HOME/.profile"; do
    [ -f "$rc" ] || continue
    root=$(grep -oP '(?<=export g09root=)\S+' "$rc" 2>/dev/null | tail -1)
    [ -n "$root" ] || continue
    scr=$(grep -oP '(?<=export GAUSS_SCRDIR=)\S+' "$rc" 2>/dev/null | tail -1)
    echo "Found Gaussian settings in $rc: g09root=$root GAUSS_SCRDIR=${scr:-/scratch}"
    if [ "$INTERACTIVE" = "1" ]; then
      read -rp "Use these? [Y/n]: " USE_FOUND
    else
      USE_FOUND="y"
    fi
    if [[ ! "$USE_FOUND" =~ ^[Nn] ]]; then
      FOUND_G09ROOT="$root"
      FOUND_GAUSS_SCRDIR="${scr:-/scratch}"
      return
    fi
    break
  done

  if [ "$INTERACTIVE" = "1" ]; then
    read -rp "Enter g09root (Gaussian install root, e.g. /home/$USER) [skip]: " FOUND_G09ROOT
    if [ -n "$FOUND_G09ROOT" ]; then
      read -rp "Enter GAUSS_SCRDIR (Gaussian scratch dir) [/scratch]: " FOUND_GAUSS_SCRDIR
      FOUND_GAUSS_SCRDIR="${FOUND_GAUSS_SCRDIR:-/scratch}"
    fi
  else
    echo "No Gaussian environment found and none given non-interactively -- skipping" \
         "(Gaussian-dependent features won't work until you edit the generated .env files)." >&2
  fi
}

write_gaussian_env_file() {
  local dir="$1"
  [ -n "$FOUND_G09ROOT" ] || return 0
  local exedir="$FOUND_G09ROOT/g09/bsd:$FOUND_G09ROOT/g09/local:$FOUND_G09ROOT/g09/extras:$FOUND_G09ROOT/g09"
  {
    echo "g09root=$FOUND_G09ROOT"
    echo "GAUSS_EXEDIR=\"$exedir\""
    echo "GAUSS_SCRDIR=$FOUND_GAUSS_SCRDIR"
  } >> "$dir/.env"
}

# ------------------------------------------------------- GAMESS paths

FOUND_RUNGMS=""
FOUND_GAMESS_SCRATCH=""

detect_or_prompt_gamess_paths() {
  FOUND_RUNGMS="${BOT_RUNGMS_PATH:-}"
  FOUND_GAMESS_SCRATCH="${BOT_GAMESS_SCRATCH_DIR:-}"

  if [ -z "$FOUND_RUNGMS" ] && command -v rungms >/dev/null 2>&1; then
    FOUND_RUNGMS="$(command -v rungms)"
    echo "Found rungms on PATH: $FOUND_RUNGMS"
  fi

  if [ -z "$FOUND_RUNGMS" ] && [ "$INTERACTIVE" = "1" ]; then
    read -rp "Path to rungms (GAMESS executable script) [skip, set later]: " FOUND_RUNGMS
  fi

  if [ -z "$FOUND_GAMESS_SCRATCH" ]; then
    local default_scratch="$HOME/gamess_scratch"
    if [ "$INTERACTIVE" = "1" ]; then
      read -rp "GAMESS scratch directory [$default_scratch]: " FOUND_GAMESS_SCRATCH
    fi
    FOUND_GAMESS_SCRATCH="${FOUND_GAMESS_SCRATCH:-$default_scratch}"
  fi
  mkdir -p "$FOUND_GAMESS_SCRATCH"
}

# --------------------------------------------------------- venv/install

setup_venv() {
  local dir="$1"
  echo "--- Setting up venv in $dir/venv ---"
  "$PYTHON_BIN" -m venv "$dir/venv"
  "$dir/venv/bin/pip" install --upgrade pip -q
}

install_gaussbot() {
  setup_venv "$GAUSSBOT_DIR"
  "$GAUSSBOT_DIR/venv/bin/pip" install -q -r "$GAUSSBOT_DIR/requirements.txt"
  "$GAUSSBOT_DIR/venv/bin/pip" install -q -e "$GAUSSBOT_DIR"
  : > "$GAUSSBOT_DIR/.env"
  write_gaussian_env_file "$GAUSSBOT_DIR"
  chmod +x "$GAUSSBOT_DIR/run_gui.sh"
  echo "gaussbot installed."
}

install_gamessbot() {
  setup_venv "$GAMESSBOT_DIR"
  # gaussbot is a sibling local package, not on PyPI -- install it
  # into gamessbot's own venv explicitly, from wherever THIS script
  # found it on THIS machine, rather than a path baked into any
  # config file (see gamessbot/pyproject.toml's comment on this).
  "$GAMESSBOT_DIR/venv/bin/pip" install -q -e "$GAUSSBOT_DIR"
  "$GAMESSBOT_DIR/venv/bin/pip" install -q -r "$GAMESSBOT_DIR/requirements.txt"
  "$GAMESSBOT_DIR/venv/bin/pip" install -q -e "$GAMESSBOT_DIR"
  : > "$GAMESSBOT_DIR/.env"
  write_gaussian_env_file "$GAMESSBOT_DIR"
  {
    [ -n "$FOUND_RUNGMS" ] && echo "GAMESSBOT_RUNGMS_PATH=$FOUND_RUNGMS"
    [ -n "$FOUND_GAMESS_SCRATCH" ] && echo "GAMESSBOT_SCRATCH_DIR=$FOUND_GAMESS_SCRATCH"
  } >> "$GAMESSBOT_DIR/.env"
  chmod +x "$GAMESSBOT_DIR/run_gui.sh"
  echo "gamessbot installed."
}

# ------------------------------------------------------------- run it

if [ "$INSTALL_GAUSSBOT" = "1" ]; then
  echo
  echo "--- Gaussian environment (for gaussbot) ---"
  detect_or_prompt_gaussian_env
fi
if [ "$INSTALL_GAMESSBOT" = "1" ]; then
  echo
  echo "--- GAMESS paths (for gamessbot) ---"
  detect_or_prompt_gamess_paths
  if [ "$INSTALL_GAUSSBOT" != "1" ]; then
    echo
    echo "--- Gaussian environment (for gamessbot's optional guess-geometry intake) ---"
    detect_or_prompt_gaussian_env
  fi
fi

echo
[ "$INSTALL_GAUSSBOT" = "1" ] && install_gaussbot
[ "$INSTALL_GAMESSBOT" = "1" ] && install_gamessbot

chmod +x "$SCRIPT_DIR/home/serve.py" 2>/dev/null || true

echo
echo "=== Done ==="
[ "$INSTALL_GAUSSBOT" = "1" ] && echo "  gaussbot:   $GAUSSBOT_DIR/run_gui.sh   (http://127.0.0.1:8765)"
[ "$INSTALL_GAMESSBOT" = "1" ] && echo "  gamessbot:  $GAMESSBOT_DIR/run_gui.sh   (http://127.0.0.1:8766)"
echo "  home page:  python3 $SCRIPT_DIR/home/serve.py   (http://127.0.0.1:8764, links to both)"
echo
echo "Machine-specific settings (Gaussian env, GAMESS paths, ports) live in each"
echo "bot's own .env file -- edit those directly any time, no need to reinstall."
