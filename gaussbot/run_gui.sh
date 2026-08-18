#!/bin/bash
# Launches the gaussbot GUI. Portable: resolves its own location at
# runtime instead of a hardcoded path, so this works wherever the BOT
# directory was copied to. Machine-specific settings (Gaussian env
# vars, port) live in ./.env (written by install.sh, or edit it by
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
exec "$SCRIPT_DIR/venv/bin/python" -m uvicorn gaussbot.webapp:app --host 127.0.0.1 --port "${GAUSSBOT_GUI_PORT:-8765}"
