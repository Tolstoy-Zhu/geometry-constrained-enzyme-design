"""Retrying, cached HTTP access for external biological databases."""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path
from typing import Any

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .io_utils import ensure_dir


class CachedSession:
    """A requests session with retry policy and stale-cache fallback."""

    def __init__(
        self,
        cache_dir: str | Path,
        timeout: int = 30,
        max_age_seconds: int = 30 * 24 * 60 * 60,
    ) -> None:
        self.cache_dir = ensure_dir(cache_dir)
        self.timeout = timeout
        self.max_age_seconds = max_age_seconds
        self.last_error = ""
        retry = Retry(
            total=4,
            connect=4,
            read=4,
            status=4,
            backoff_factor=0.6,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=frozenset({"GET"}),
            respect_retry_after_header=True,
        )
        adapter = HTTPAdapter(max_retries=retry)
        self.session = requests.Session()
        self.session.mount("https://", adapter)
        self.session.mount("http://", adapter)
        self.session.headers.update(
            {"User-Agent": "AquaEsterNet-StructMine/0.4.0 (research workflow)"}
        )

    def _cache_path(self, url: str, params: dict[str, Any] | None, kind: str) -> Path:
        identity = json.dumps(
            {"url": url, "params": params or {}},
            sort_keys=True,
            ensure_ascii=True,
            default=str,
        )
        key = hashlib.sha256(identity.encode("utf-8")).hexdigest()
        return self.cache_dir / f"{key}.{kind}"

    def _fresh(self, path: Path) -> bool:
        return path.exists() and time.time() - path.stat().st_mtime <= self.max_age_seconds

    @staticmethod
    def _write_atomic(path: Path, text: str) -> None:
        tmp = path.with_suffix(path.suffix + ".tmp")
        tmp.write_text(text, encoding="utf-8")
        tmp.replace(path)

    def get_json(self, url: str, params: dict[str, Any] | None = None) -> Any | None:
        cache_path = self._cache_path(url, params, "json")
        if self._fresh(cache_path):
            try:
                return json.loads(cache_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                pass

        self.last_error = ""
        try:
            response = self.session.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            self._write_atomic(cache_path, json.dumps(data, ensure_ascii=False))
            return data
        except (requests.RequestException, ValueError, OSError) as exc:
            self.last_error = f"{type(exc).__name__}: {exc}"
            if cache_path.exists():
                try:
                    return json.loads(cache_path.read_text(encoding="utf-8"))
                except (OSError, json.JSONDecodeError):
                    pass
            return None

    def get_text(self, url: str, params: dict[str, Any] | None = None) -> str | None:
        cache_path = self._cache_path(url, params, "txt")
        if self._fresh(cache_path):
            try:
                return cache_path.read_text(encoding="utf-8")
            except OSError:
                pass

        self.last_error = ""
        try:
            response = self.session.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            text = response.text
            self._write_atomic(cache_path, text)
            return text
        except (requests.RequestException, OSError) as exc:
            self.last_error = f"{type(exc).__name__}: {exc}"
            if cache_path.exists():
                try:
                    return cache_path.read_text(encoding="utf-8")
                except OSError:
                    pass
            return None
