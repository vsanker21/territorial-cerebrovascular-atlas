"""Tiny CORS-open receiver to accept authenticated NITRC file bytes from the browser."""
from __future__ import annotations

import sys
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path
from urllib.parse import unquote

OUT_DIR = Path(__file__).resolve().parents[2] / "data" / "raw" / "arterial_atlas" / "nitrc"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PORT = 8765


class Handler(BaseHTTPRequestHandler):
    def _cors(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "*")

    def do_OPTIONS(self):
        self.send_response(204)
        self._cors()
        self.end_headers()

    def do_POST(self):
        name = unquote(self.path.lstrip("/")) or "download.bin"
        # prevent path traversal
        name = Path(name).name
        length = int(self.headers.get("Content-Length", "0"))
        data = self.rfile.read(length)
        out = OUT_DIR / name
        out.write_bytes(data)
        msg = f"Wrote {out} ({len(data)} bytes)\n"
        print(msg, flush=True)
        body = msg.encode()
        self.send_response(200)
        self._cors()
        self.send_header("Content-Type", "text/plain")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, fmt, *args):
        sys.stderr.write("%s - %s\n" % (self.address_string(), fmt % args))


if __name__ == "__main__":
    print(f"Listening on http://127.0.0.1:{PORT}/  -> {OUT_DIR}", flush=True)
    HTTPServer(("127.0.0.1", PORT), Handler).serve_forever()
