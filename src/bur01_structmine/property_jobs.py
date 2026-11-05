"""One inference queue per web server. Closing a browser tab does not cancel its job."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, replace
from pathlib import Path
from threading import Lock

from .property_screen import profiles_root
from .run_context import make_run_id


@dataclass(frozen=True)
class PropertyJob:
    key: str
    title: str
    state: str
    message: str


class PropertyJobs:
    def __init__(self):
        self._pool = ThreadPoolExecutor(max_workers=1, thread_name_prefix="aqua-property")
        self._jobs = {}
        self._lock = Lock()

    def snapshot(self):
        with self._lock:
            return list(reversed(self._jobs.values()))

    def submit(self, title, predict, **inputs):
        key = make_run_id("web")
        output = profiles_root() / key
        with self._lock:
            self._jobs[key] = PropertyJob(key, title, "queued", "等待本机计算资源")
        self._pool.submit(self._run, key, predict, inputs, output)
        return key

    def _update(self, key, state, message):
        with self._lock:
            self._jobs[key] = replace(self._jobs[key], state=state, message=message)

    def _run(self, key, predict, inputs, output: Path):
        def progress(message):
            self._update(key, "running", message.split(" / ", 1)[0])

        try:
            progress("准备模型输入")
            predict(**inputs, output=output, progress=progress)
        except Exception as exc:
            self._update(key, "failed", str(exc))
        else:
            self._update(key, "complete", "预测结果已保存")
