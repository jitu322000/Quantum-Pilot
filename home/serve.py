#!/usr/bin/env python3
"""
serve.py

Tiny static-file server for the Quantum Pilot landing page -- no
framework needed (just stdlib http.server), since this page is
nothing more than two links out to gaussbot/gamessbot's own
already-running GUIs. Opens a browser tab automatically, same
convenience pattern gaussbot-gui/gamessbot-gui use.

Also serves GET /config -- reads GAUSSBOT_GUI_PORT/GAMESSBOT_GUI_PORT
out of each sibling bot's own .env file (written by install.sh, or
edited by hand), so index.html's links always point at whatever port
each bot is actually configured for instead of a hardcoded default --
per your report that the gamessbot card was pointing at the wrong
port after you moved it to 8767.

Run: python3 serve.py   (opens a browser tab at http://127.0.0.1:8764)
"""

import http.server
import json
import os
import re
import socketserver
import threading
import webbrowser

PORT = int(os.environ.get("BOT_HOME_PORT", "8764"))
HOME_DIR = os.path.dirname(os.path.abspath(__file__))
BOT_ROOT = os.path.dirname(HOME_DIR)


def _port_from_env_file(bot_dir: str, var_name: str, default: int) -> int:
    env_path = os.path.join(BOT_ROOT, bot_dir, ".env")
    try:
        with open(env_path) as f:
            text = f.read()
    except OSError:
        return default
    m = re.search(rf"^{var_name}=(\d+)\s*$", text, re.MULTILINE)
    return int(m.group(1)) if m else default


class Handler(http.server.SimpleHTTPRequestHandler):
    def do_GET(self):
        if self.path == "/config":
            body = json.dumps({
                "gaussbot_port": _port_from_env_file("gaussbot", "GAUSSBOT_GUI_PORT", 8765),
                "gamessbot_port": _port_from_env_file("gamessbot", "GAMESSBOT_GUI_PORT", 8766),
            }).encode()
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return
        super().do_GET()


def main() -> None:
    os.chdir(HOME_DIR)
    with socketserver.TCPServer(("127.0.0.1", PORT), Handler) as httpd:
        url = f"http://127.0.0.1:{PORT}"
        threading.Timer(1.0, lambda: webbrowser.open(url)).start()
        print(f"Quantum Pilot landing page running at {url} (Ctrl+C to stop)")
        httpd.serve_forever()


if __name__ == "__main__":
    main()
