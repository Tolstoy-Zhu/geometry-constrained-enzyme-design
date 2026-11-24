"""Refresh a small number of remote jobs without blocking Streamlit rendering."""
from __future__ import annotations
import threading,time
from concurrent.futures import ThreadPoolExecutor
from .store import Store
from .execution import RemoteExecutor,server_profiles
_pool=ThreadPoolExecutor(max_workers=2,thread_name_prefix="aqua-status")
_lock=threading.Lock();_pending={};_last={}

def schedule_sync(project_id,store=None):
    store=store or Store()
    profiles=server_profiles()
    import json
    for row in store.tasks(project_id):
        spec=json.loads(row["spec"])
        if spec["executor"]=="local" or spec["executor"] not in profiles:continue
        result=json.loads(row["result"]) if row.get("result") else {}
        needs_download=row["state"] in {"completed","failed","cancelled"} and not result.get("local_output")
        if row["state"] not in {"queued","running","blocked","interrupted"} and not needs_download:continue
        key=(str(store.path),row["id"])
        with _lock:
            if key in _pending or time.monotonic()-_last.get(key,0)<15:continue
            _last[key]=time.monotonic()
            def refresh(task_id=row["id"],profile=profiles[spec["executor"]],key=key):
                try:
                    # A blocked submission is never resent automatically. The explicit reconnect
                    # action reconciles submission identity before retrying a transport failure.
                    remote=RemoteExecutor(profile)
                    remote.sync(task_id,store)
                except (OSError,ValueError,RuntimeError,TimeoutError):
                    pass  # Existing state and last log remain readable during disconnection.
                finally:
                    with _lock:_pending.pop(key,None)
            _pending[key]=_pool.submit(refresh)

