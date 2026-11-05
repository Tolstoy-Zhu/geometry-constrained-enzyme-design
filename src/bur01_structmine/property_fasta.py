"""Sequence-only inputs for the local property models; coordinates are optional."""

from __future__ import annotations

import re

from .prediction_sequence import FastaRecord
from .property_inference import CANONICAL


def parse_property_fasta(text: str) -> list[FastaRecord]:
    if len(text) > 2 * 1024 * 1024:
        raise ValueError("序列输入超过 2 MiB，请分批提交。")
    lines = [s.strip() for s in text.lstrip("\ufeff").splitlines() if s.strip()]
    if not lines:
        raise ValueError("请粘贴蛋白序列，或上传 FASTA 文件。")
    if not lines[0].startswith(">"):
        if any(s.startswith(">") for s in lines):
            raise ValueError("多条序列的每条记录都需要以 >名称 开头。")
        lines.insert(0, ">sequence_1")
    records, names = [], set()
    for block in "\n".join(lines)[1:].split("\n>"):
        header, _, residues = block.partition("\n")
        name = header.split(maxsplit=1)[0] if header.strip() else ""
        if not re.fullmatch(r"[\w][\w.:|+\-]{0,119}", name):
            raise ValueError("每条 FASTA 记录需要一个有效名称（首个词，最长 120 个字符）。")
        if name in names:
            raise ValueError(f"序列名称重复：{name}。请为每条记录使用不同名称。")
        sequence = "".join(residues.split()).upper()
        if not 1 <= len(sequence) <= 2500:
            raise ValueError(f"{name} 长度为 {len(sequence)}，本次联合预测支持 1–2500 个残基。")
        unsupported = sorted(set(sequence) - CANONICAL)
        if unsupported:
            raise ValueError(f"{name} 含不支持的残基符号 {', '.join(unsupported)}，请使用 20 种标准氨基酸。")
        names.add(name)
        records.append(FastaRecord(name, sequence))
        if len(records) > 100:
            raise ValueError("每次最多接收 100 条序列，请分批提交。")
    return records


def collect_fasta(text: str, first_as_parent: bool):
    records = parse_property_fasta(text)
    reference = records[0]
    if first_as_parent and len(records) < 2:
        raise ValueError("亲本对照需要首条亲本和至少一条候选序列。")
    sequences, samples, members = {}, {}, []
    for index, record in enumerate(records):
        digest = record.sha256
        sequences[digest] = record.sequence
        samples.setdefault(digest, dict(id=record.name, parent=first_as_parent and index == 0, runs=["FASTA"]))
        if first_as_parent and index:
            if len(reference.sequence) != len(record.sequence):
                raise ValueError(f"{record.name} 与首条亲本长度不同。含插入或缺失时请选“分别预测”。")
            changes = [
                f"{a}{i}{b}" for i, (a, b) in enumerate(zip(reference.sequence, record.sequence), 1) if a != b
            ]
            members.append(
                dict(
                    run="FASTA",
                    candidate_id=record.name,
                    mutations_pdb="",
                    mutations_sequence=",".join(changes),
                    mutations_8hfw="",
                    substitutions=len(changes),
                    sequence_sha256=digest,
                    parent_sha256=reference.sha256,
                )
            )
    sources = {
        "FASTA": {
            "input_kind": "fasta",
            "numbering_reference": "input_sequence",
            "comparison": "first_as_parent" if first_as_parent else "absolute",
            "records": [{"id": r.name, "sequence_sha256": r.sha256} for r in records],
        }
    }
    return sequences, samples, members, sources
