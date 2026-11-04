"""Pinned upstream resources for the offline property predictors."""

from __future__ import annotations

import hashlib
import os
import subprocess
import sys
import tempfile
import time
import venv
from pathlib import Path

import requests

from .config import project_root

SWI_REV = "c28f26d94aa6285ff08270709fc8a95806ddb1db"
TM_REV = "02d4c12b39ffae8a6cde64a401210390ea9cb23f"
ESM_REV = "08e4846e537177426273712802403f7ba8261b6c"
ENZ_REV = "bf9d5388e07d719b94bad3aabc010b78bcd1e2c5"
SOURCES = {
    "swi": f"https://raw.githubusercontent.com/Gardner-BinfLab/SoDoPE_paper_2020/{SWI_REV}/SWI/",
    "adapter": f"https://raw.githubusercontent.com/loschmidt/TmProt/{TM_REV}/tmprot-1.0/src/tmprot/model/",
    "esm": f"https://huggingface.co/facebook/esm2_t33_650M_UR50D/resolve/{ESM_REV}/",
    "enzyme": f"https://raw.githubusercontent.com/maranasgroup/EnzRank/{ENZ_REV}/AICHE/CNN_results_split_final/Final_model.model/",
}
ASSETS = {
    "swi/swi.py": "6e9c0fdb18f65e20c885f5d59b7f5b2cebf557f8c471adb74713be1848e658f3",
    "adapter/adapter_config.json": "7b59b6d27d71cabf861a4d7a400f80f80d3ec3d0ee1b54d913f8d81796d0f3a8",
    "adapter/adapter_model.safetensors": "d7bef27e36de848973fe6b78864882f899c2ad59bd17100d7e9d0eac4c24ff85",
    "esm/config.json": "539095c22efc52a09d6147074ba4ca119f76a890df5901213b2b55f7d2f96b2b",
    "esm/vocab.txt": "0b82cc0a7c7cf9e567b1e5892d793285b9fbae822c964ca48696f7db44598e03",
    "esm/tokenizer_config.json": "7e9161ecdb548ec45a41cbc6b24aa4476fdd418461f491c4207baa99419a29ad",
    "esm/special_tokens_map.json": "3aedcd4211c0d43aec4e607ff60a63255f3174ead795e997350f09a5f8cd9ee1",
    "esm/model.safetensors": "a08adabb949fa67ad3c14b509d04fd60368b35007b0095e3358f81200c4f4db0",
    "enzyme/keras_metadata.pb": "f497a391e0fe22e03bbb3b32e916dafc583ebe24601b8ded19966e739d02415a",
    "enzyme/saved_model.pb": "e3fe10e07fd9b5c93ef42f64307de4488ead986de8aae77c5305c468e040af0a",
    "enzyme/variables/variables.index": "070edd9c1a62a58edd61c6ea5a8812d4d807a7807da0b67e246a483e3c3b7ad4",
    "enzyme/variables/variables.data-00000-of-00001": "ac0f1bfce59ef4063b4a34bc61a496bf1fb6afd9d4a449a48591b18f05e54c62",
}


def asset_root():
    return Path(os.environ.get("AQUA_PROPERTY_MODELS") or project_root() / ".tools/aqua-properties/v1")


def runtime_python():
    configured = os.environ.get("AQUA_PROPERTY_PYTHON")
    if configured:
        return Path(configured)
    return (
        project_root()
        / ".tools/property-runtime"
        / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    )


def file_hash(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_assets():
    directory = asset_root()
    for name, expected in ASSETS.items():
        target = directory / name
        if not target.is_file() or file_hash(target) != expected:
            raise ValueError(f"Missing or altered property model asset: {name}. Run properties-setup.")
    return directory


def _download_backbone(url, path):
    # A multi-GB response may close partway through on a campus proxy. Retry only the current range.
    size = 2609506392
    with path.open("wb") as f:
        for start in range(0, size, 64 * 1024 * 1024):
            end = min(start + 64 * 1024 * 1024, size) - 1
            for attempt in range(3):
                try:
                    with requests.get(
                        url,
                        params={"download": "true", "part": start},
                        headers={"Range": f"bytes={start}-{end}", "Accept-Encoding": "identity"},
                        stream=True,
                        timeout=(20, 90),
                    ) as r:
                        r.raise_for_status()
                        if (
                            r.status_code != 206
                            or r.headers.get("Content-Range") != f"bytes {start}-{end}/{size}"
                        ):
                            raise OSError("Unexpected range response for ESM-2.")
                        f.seek(start)
                        f.truncate()
                        for block in r.iter_content(1024 * 1024):
                            f.write(block)
                        if f.tell() != end + 1:
                            raise OSError("Incomplete ESM-2 download range.")
                    break
                except (requests.RequestException, OSError):
                    if attempt == 2:
                        raise
                    time.sleep(2**attempt)


def install_assets():
    directory = asset_root()
    for name, expected in ASSETS.items():
        target = directory / name
        if target.is_file() and file_hash(target) == expected:
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        group, relative = name.split("/", 1)
        with tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as f:
            tmp = Path(f.name)
        try:
            url = SOURCES[group] + relative
            if name == "esm/model.safetensors":
                _download_backbone(url, tmp)
            else:
                with requests.get(url, stream=True, timeout=(20, 120)) as r:
                    r.raise_for_status()
                    with tmp.open("wb") as out:
                        for block in r.iter_content(1024 * 1024):
                            out.write(block)
            if file_hash(tmp) != expected:
                raise ValueError(f"Downloaded bytes do not match the pinned SHA-256 for {name}.")
            os.replace(tmp, target)
        finally:
            tmp.unlink(missing_ok=True)
    return directory


def install_runtime():
    if os.environ.get("AQUA_PROPERTY_PYTHON"):
        raise ValueError(
            "AQUA_PROPERTY_PYTHON is set. Install requirements-properties.txt in that environment."
        )
    if sys.version_info[:2] not in {(3, 10), (3, 11)}:
        raise ValueError("Use Python 3.10/3.11 to create the isolated property environment.")
    python = runtime_python()
    if not python.is_file():
        venv.EnvBuilder(with_pip=True).create(python.parent.parent)
    requirements = project_root() / "requirements-properties.txt"
    if sys.platform != "darwin":
        subprocess.run(
            [
                str(python),
                "-m",
                "pip",
                "install",
                "torch==2.9.1+cpu",
                "--index-url",
                "https://download.pytorch.org/whl/cpu",
                "--extra-index-url",
                "https://pypi.org/simple",
            ],
            check=True,
        )
    subprocess.run([str(python), "-m", "pip", "install", "-r", str(requirements)], check=True)
    return python
