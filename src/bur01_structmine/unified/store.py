"""SQLite metadata, atomic claims and immutable task specifications."""
from __future__ import annotations
import hashlib
import json
import os
import sqlite3
import uuid
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from .contracts import TaskSpec

def now():
    return datetime.now(timezone.utc).isoformat()

def state_root():
    from ..config import project_root
    return Path(os.environ.get("AQUA_STATE_ROOT", project_root() / ".cache/unified")).resolve()

def uid(prefix):
    return f"{prefix}_{uuid.uuid4().hex}"

def dumps(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, allow_nan=False)

class Store:
    def __init__(self, path=None):
        self.path = Path(path or state_root() / "projects.sqlite")
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.connect() as db:
            db.executescript("""
            PRAGMA journal_mode=WAL;
            CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
            INSERT OR IGNORE INTO metadata VALUES ('schema_version','1');
            CREATE TABLE IF NOT EXISTS projects (
                id TEXT PRIMARY KEY, name TEXT NOT NULL, family_id TEXT NOT NULL,
                executor TEXT NOT NULL DEFAULT 'local', created TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS records (
                kind TEXT NOT NULL, id TEXT NOT NULL, project_id TEXT NOT NULL,
                payload TEXT NOT NULL, created TEXT NOT NULL, PRIMARY KEY(kind,id,project_id));
            CREATE TABLE IF NOT EXISTS tasks (
                id TEXT PRIMARY KEY, project_id TEXT NOT NULL, spec TEXT NOT NULL,
                spec_hash TEXT NOT NULL, state TEXT NOT NULL, message TEXT NOT NULL DEFAULT '',
                result TEXT, created TEXT NOT NULL, updated TEXT NOT NULL,
                cancel_requested INTEGER NOT NULL DEFAULT 0);
            CREATE TABLE IF NOT EXISTS events (
                seq INTEGER PRIMARY KEY AUTOINCREMENT, task_id TEXT NOT NULL, time TEXT NOT NULL,
                state TEXT NOT NULL, message TEXT NOT NULL);
            CREATE INDEX IF NOT EXISTS tasks_state ON tasks(state,created);
            """)
            version = db.execute("SELECT value FROM metadata WHERE key='schema_version'").fetchone()[0]
            if version != "1":
                raise ValueError("Unsupported unified database schema")
            db.execute("INSERT OR IGNORE INTO projects VALUES (?,?,?,?,?)",
                       ("default", "默认项目", "bur01", "local", now()))

    @contextmanager
    def connect(self):
        db = sqlite3.connect(self.path, timeout=30)
        db.row_factory = sqlite3.Row
        try:
            with db:
                yield db
        finally:
            db.close()

    def projects(self):
        with self.connect() as db:
            return [dict(r) for r in db.execute("SELECT * FROM projects ORDER BY created")]

    def create_project(self, name, family_id="bur01", executor="local"):
        if not name.strip():
            raise ValueError("Project name is required")
        project_id = uid("project")
        with self.connect() as db:
            db.execute("INSERT INTO projects VALUES (?,?,?,?,?)", (project_id, name.strip(), family_id, executor, now()))
        return project_id

    def record(self, kind, record_id, project_id, payload):
        """Identical imports are idempotent; conflicting scientific records never overwrite."""
        encoded = dumps(payload)
        with self.connect() as db:
            prior = db.execute("SELECT payload FROM records WHERE kind=? AND id=? AND project_id=?",
                               (kind, record_id, project_id)).fetchone()
            if prior and prior[0] != encoded:
                raise ValueError(f"Conflicting {kind} record: {record_id}")
            db.execute("INSERT OR IGNORE INTO records VALUES (?,?,?,?,?)",
                       (kind, record_id, project_id, encoded, now()))

    def records(self, kind, project_id="default"):
        with self.connect() as db:
            return [json.loads(r[0]) for r in db.execute(
                "SELECT payload FROM records WHERE kind=? AND project_id=? ORDER BY created",
                (kind, project_id))]

    def submit(self, spec: TaskSpec):
        encoded = dumps(spec.model_dump(mode="json"))
        digest = hashlib.sha256(encoded.encode()).hexdigest()
        with self.connect() as db:
            prior = db.execute("SELECT spec_hash FROM tasks WHERE id=?", (spec.task_id,)).fetchone()
            if prior and prior[0] != digest:
                raise ValueError("Task ID already belongs to different inputs")
            if not prior:
                db.execute("INSERT INTO tasks(id,project_id,spec,spec_hash,state,created,updated) VALUES(?,?,?,?,?,?,?)",
                           (spec.task_id,spec.project_id,encoded,digest,"queued",now(),now()))
                db.execute("INSERT INTO events(task_id,time,state,message) VALUES(?,?,?,?)",
                           (spec.task_id,now(),"queued","任务已保存"))
        return spec.task_id

    def tasks(self, project_id=None):
        with self.connect() as db:
            sql = "SELECT * FROM tasks" + (" WHERE project_id=?" if project_id else "") + " ORDER BY created DESC"
            rows = db.execute(sql, (project_id,) if project_id else ()).fetchall()
            return [dict(r) for r in rows]

    def task(self, task_id):
        with self.connect() as db:
            row = db.execute("SELECT * FROM tasks WHERE id=?", (task_id,)).fetchone()
            if row is None:
                raise ValueError("Unknown task: " + task_id)
            return dict(row)

    def claim(self):
        with self.connect() as db:
            db.execute("BEGIN IMMEDIATE")
            if db.execute("SELECT 1 FROM tasks WHERE state='running' AND json_extract(spec,'$.executor')='local'").fetchone():
                return None
            row = db.execute("SELECT * FROM tasks WHERE state='queued' AND cancel_requested=0 AND json_extract(spec,'$.executor')='local' ORDER BY created LIMIT 1").fetchone()
            if row is None:
                return None
            db.execute("UPDATE tasks SET state='running',updated=? WHERE id=?", (now(),row["id"]))
            db.execute("INSERT INTO events(task_id,time,state,message) VALUES(?,?,?,?)",
                       (row["id"],now(),"running","开始执行"))
            return TaskSpec.model_validate_json(row["spec"])

    def update(self, task_id, state, message="", result=None):
        if state not in {"queued","running","completed","failed","cancelled","interrupted","blocked"}:
            raise ValueError("Invalid task state")
        with self.connect() as db:
            if isinstance(result,dict) and result.get('remote_updated'):
                # Polling and result downloads can finish concurrently. Merge the
                # download receipt under the same write lock instead of losing it.
                db.execute('BEGIN IMMEDIATE')
                row=db.execute('SELECT result FROM tasks WHERE id=?',(task_id,)).fetchone()
                prior=json.loads(row[0]) if row and row[0] else {}
                previous_time=prior.get('remote_updated','')
                if previous_time>result['remote_updated']:return
                result=dict(result)
                if previous_time==result['remote_updated'] and prior.get('local_output') and not result.get('local_output'):
                    if (Path(prior['local_output'])/'.aqua-download.json').is_file():result['local_output']=prior['local_output']
                if previous_time==result['remote_updated'] and prior.get('download_error') and not result.get('local_output'):
                    result.setdefault('download_error',prior['download_error'])
            db.execute("UPDATE tasks SET state=?,message=?,result=COALESCE(?,result),updated=? WHERE id=?",
                       (state,message,dumps(result) if result is not None else None,now(),task_id))
            db.execute("INSERT INTO events(task_id,time,state,message) VALUES(?,?,?,?)",
                       (task_id,now(),state,message))

    def cancel(self, task_id):
        with self.connect() as db:
            row=db.execute("SELECT state FROM tasks WHERE id=?",(task_id,)).fetchone()
            if row is None:raise ValueError("Unknown task: "+task_id)
            if row[0] in {"completed","failed","interrupted","cancelled"}:return
            db.execute("UPDATE tasks SET cancel_requested=1, state=CASE WHEN state='queued' THEN 'cancelled' ELSE state END, updated=? WHERE id=?",
                       (now(),task_id))
            db.execute("INSERT INTO events(task_id,time,state,message) VALUES(?,?,?,?)",
                       (task_id,now(),"cancelled" if row[0]=="queued" else row[0],"已申请取消此任务"))

    def recover(self):
        """Recover only local jobs whose execution lock is no longer held."""
        from .execution import worker_lock
        with self.connect() as db:
            rows=db.execute("SELECT id FROM tasks WHERE state='running' AND json_extract(spec,'$.executor')='local'").fetchall()
        for row in rows:
            try:
                with worker_lock(self.path.parent/'tasks'/row[0]):
                    self.update(row[0],"interrupted","上次执行进程已结束；可检查输出并重试")
            except RuntimeError:
                pass  # An orphaned but live job must not be submitted twice.

    def retry(self, task_id):
        with self.connect() as db:
            db.execute("BEGIN IMMEDIATE")
            old=db.execute("SELECT * FROM tasks WHERE id=?",(task_id,)).fetchone()
            if old is None:raise ValueError("Unknown task: "+task_id)
            spec=TaskSpec.model_validate_json(old["spec"])
            if spec.kind in {"mechanism","features"}:
                # A repeated retry request reconciles the same checkpoint directory.
                if old["state"] in {"queued","running","completed"}:return task_id
                if old["state"] not in {"failed","interrupted","blocked","cancelled"}:
                    raise ValueError("Only a stopped task can be retried")
                db.execute("UPDATE tasks SET state='queued',cancel_requested=0,message=?,updated=? WHERE id=?",
                           ("从已保存的进度恢复",now(),task_id))
                db.execute("INSERT INTO events(task_id,time,state,message) VALUES(?,?,?,?)",
                           (task_id,now(),"queued","恢复原任务；保留已完成的分段"))
                return task_id
            if old["state"] not in {"failed","interrupted","blocked","cancelled"}:
                raise ValueError("Only a stopped task can be retried")
            prior=db.execute("SELECT id FROM tasks WHERE json_extract(spec,'$.parent_task_id')=? AND state IN ('queued','running','completed') ORDER BY created DESC LIMIT 1",(task_id,)).fetchone()
            if prior:return prior[0]
            spec.parent_task_id=spec.task_id
            spec.task_id=uid("task")
            encoded=dumps(spec.model_dump(mode="json"))
            digest=hashlib.sha256(encoded.encode()).hexdigest()
            db.execute("INSERT INTO tasks(id,project_id,spec,spec_hash,state,created,updated) VALUES(?,?,?,?,?,?,?)",
                       (spec.task_id,spec.project_id,encoded,digest,"queued",now(),now()))
            db.execute("INSERT INTO events(task_id,time,state,message) VALUES(?,?,?,?)",
                       (spec.task_id,now(),"queued","重试已保存的任务"))
            return spec.task_id

def index_history(root, store=None, project_id="default"):
    root=Path(root).resolve()
    store=store or Store()
    count=0
    for area, name in (("runs","run_manifest.json"),("nr_runs","nr_run_manifest.json"),
                       ("predictions","prediction_manifest.json"),("designs","manifest.json"),
                       ("properties","manifest.json")):
        for path in sorted((root/"results"/area).glob("*/"+name)):
            data=path.read_bytes()
            try:
                obj=json.loads(data)
            except (ValueError,UnicodeDecodeError):
                continue
            record={"area":area,"run_id":path.parent.name,"path":str(path),
                    "manifest_sha256":hashlib.sha256(data).hexdigest(),
                    "status":obj.get("status","saved"),"source_root":str(root)}
            store.record("history",hashlib.sha256(str(path).encode()+data).hexdigest(),project_id,record)
            count+=1
    return count
