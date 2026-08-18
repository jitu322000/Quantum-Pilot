#!/bin/bash
# Launches the gamessbot GUI. Portable: resolves its own location at
# runtime instead of a hardcoded path, so this works wherever the BOT
# directory was copied to. Machine-specific settings (Gaussian env
# vars for the "guess geometry" intake path, GAMESS rungms/scratch
# dir, port) live in ./.env (written by install.sh, or edit it by
# hand) rather than being baked into this script.
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [ -f "$SCRIPT_DIR/.env" ]; then
  set -a
  source "$SCRIPT_DIR/.env"
  set +a
fi

export PYTHONPATH="$SCRIPT_DIR"
exec "$SCRIPT_DIR/venv/bin/python" -m uvicorn gamessbot.webapp:app --host 127.0.0.1 --port "${GAMESSBOT_GUI_PORT:-8766}"
