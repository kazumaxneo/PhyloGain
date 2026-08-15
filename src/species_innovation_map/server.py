from __future__ import annotations

import json
import sqlite3
import webbrowser
import zlib
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse


class ProjectHandler(SimpleHTTPRequestHandler):
    project_dir: Path

    def __init__(self, *args, directory=None, **kwargs):
        super().__init__(*args, directory=str(self.project_dir), **kwargs)

    def do_GET(self):  # noqa: N802
        parsed = urlparse(self.path)
        if parsed.path.startswith("/api/"):
            try:
                payload = self._api(parsed.path, parse_qs(parsed.query))
                self._json(payload)
            except (ValueError, sqlite3.Error) as exc:
                self._json({"error": str(exc)}, status=400)
            return
        super().do_GET()

    def _api(self, path: str, query: dict[str, list[str]]):
        if path == "/api/project":
            return json.loads((self.project_dir / "project.json").read_text(encoding="utf-8"))
        with sqlite3.connect(self.project_dir / "species_map.sqlite") as connection:
            connection.row_factory = sqlite3.Row
            if path == "/api/branches":
                rows = connection.execute(
                    "SELECT * FROM branches ORDER BY branch_id"
                ).fetchall()
                return [dict(row) for row in rows]
            if path == "/api/phenotype":
                phenotype = _required(query, "id")
                rows = connection.execute(
                    "SELECT branch_id,event,parent_state,child_state FROM phenotype_events WHERE phenotype_id=?",
                    (phenotype,),
                ).fetchall()
                return [dict(row) for row in rows]
            if path == "/api/branch":
                branch = _required(query, "id")
                event = query.get("event", [""])[0]
                limit = min(max(int(query.get("limit", ["500"])[0]), 1), 5000)
                params: list[object] = [branch]
                where = "branch_id=?"
                if event in {"gain", "loss"}:
                    where += " AND event=?"
                    params.append(event)
                total = connection.execute(
                    f"SELECT COUNT(*) FROM events WHERE {where}", params
                ).fetchone()[0]
                params.append(limit)
                rows = connection.execute(
                    f"SELECT family_id,event,parent_state,child_state FROM events WHERE {where} ORDER BY event,family_id LIMIT ?",
                    params,
                ).fetchall()
                return {"total": total, "events": [dict(row) for row in rows]}
            if path == "/api/family":
                family = _required(query, "id")
                row = connection.execute(
                    "SELECT family_id,members_json FROM families WHERE family_id=?", (family,)
                ).fetchone()
                if row is None:
                    raise ValueError(f"Unknown family: {family}")
                gains = connection.execute(
                    "SELECT branch_id,event FROM events WHERE family_id=? ORDER BY branch_id",
                    (family,),
                ).fetchall()
                stored_members = row["members_json"]
                if isinstance(stored_members, bytes):
                    members = json.loads(zlib.decompress(stored_members).decode("utf-8"))
                else:
                    members = json.loads(stored_members or "{}")
                return {
                    "family_id": row["family_id"],
                    "members": members,
                    "events": [dict(item) for item in gains],
                }
            if path == "/api/candidates":
                phenotype = _required(query, "phenotype")
                limit = min(max(int(query.get("limit", ["200"])[0]), 1), 5000)
                rows = connection.execute(
                    "SELECT family_id,score,coincident_gains,phenotype_gains,family_gains FROM candidates WHERE phenotype_id=? ORDER BY score DESC,family_id LIMIT ?",
                    (phenotype, limit),
                ).fetchall()
                return [dict(row) for row in rows]
        raise ValueError(f"Unknown API endpoint: {path}")

    def _json(self, payload, status: int = 200) -> None:
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)


def _required(query: dict[str, list[str]], name: str) -> str:
    value = query.get(name, [""])[0]
    if not value:
        raise ValueError(f"Missing query parameter: {name}")
    return value


def serve_project(
    project_dir: str | Path,
    host: str = "127.0.0.1",
    port: int = 8000,
    open_browser: bool = True,
) -> None:
    directory = Path(project_dir).resolve()
    required = [directory / "index.html", directory / "project.json", directory / "species_map.sqlite"]
    missing = [path.name for path in required if not path.is_file()]
    if missing:
        raise ValueError(f"Not a Species Innovation Map directory; missing: {', '.join(missing)}")
    handler = type("BoundProjectHandler", (ProjectHandler,), {"project_dir": directory})
    server = ThreadingHTTPServer((host, port), handler)
    url = f"http://{host}:{server.server_port}/"
    print(f"Species Innovation Map: {url}")
    print("Press Ctrl+C to stop.")
    if open_browser:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
