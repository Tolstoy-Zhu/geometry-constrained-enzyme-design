"""Validated public data contracts for AquaEsterNet-StructMine."""
from __future__ import annotations

import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


RunStatus = Literal["created", "running", "completed", "failed", "cancelled"]
EventState = Literal["started", "running", "completed", "warning", "failed", "cancelled"]
PredictionPreset = Literal["quick", "standard", "high"]
NrDatabaseScope = Literal["nr_remote_snapshot", "full_nr", "custom_mmseqs"]
NrExecutionMode = Literal["local", "job_bundle", "import"]
RUN_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")


def _validate_run_identifier(value: str | None, label: str) -> str | None:
    if value is None:
        return None
    cleaned = value.strip()
    if not RUN_IDENTIFIER.fullmatch(cleaned):
        raise ValueError(f"{label} contains unsupported path or identifier characters")
    return cleaned


class QueryConfig(BaseModel):
    model_config = ConfigDict(extra="allow")

    pdb_file: str
    fasta_file: str | None = None
    pdb_id: str | None = None
    chain: str | None = None
    residue_range: str | None = None
    source_url: str | None = None


class FilterConfig(BaseModel):
    model_config = ConfigDict(extra="allow")

    min_alntmscore: float = Field(default=0.45, ge=0.0, le=1.0)
    max_rmsd: float = Field(default=4.0, gt=0.0)
    min_qcov: float = Field(default=0.60, ge=0.0, le=1.0)
    min_tcov: float = Field(default=0.40, ge=0.0, le=1.0)
    min_alnlen: int = Field(default=180, ge=1)
    max_sequence_identity: float = Field(default=0.95, ge=0.0, le=1.0)


class AnnotationConfig(BaseModel):
    model_config = ConfigDict(extra="allow")

    enable_rcsb: bool = True
    enable_uniprot: bool = True
    enable_ncbi: bool = True
    enable_alphafold: bool = True
    request_sleep_seconds: float = Field(default=0.15, ge=0.0)
    cache_dir: str = ".cache/annotations"
    max_records: int = Field(default=200, ge=1)


class OutputConfig(BaseModel):
    model_config = ConfigDict(extra="allow")

    out_prefix: str = "results/Bur01_vs_database"
    run_root: str = "results/runs"


class ProjectConfig(BaseModel):
    model_config = ConfigDict(extra="allow")

    query: QueryConfig
    foldseek: dict[str, Any] = Field(default_factory=dict)
    search_databases: list[dict[str, Any]] | dict[str, dict[str, Any]] | None = None
    filter: FilterConfig = Field(default_factory=FilterConfig)
    annotation: AnnotationConfig = Field(default_factory=AnnotationConfig)
    output: OutputConfig = Field(default_factory=OutputConfig)

    @field_validator("search_databases")
    @classmethod
    def validate_search_databases(cls, value: Any) -> Any:
        if value is not None and not value:
            raise ValueError("search_databases cannot be empty")
        return value


class RunRequest(BaseModel):
    config_path: str = "config.yaml"
    run_id: str | None = None
    query_structure: str | None = None
    query_chain: str | None = None
    output_root: str | None = None
    skip_search: bool = False
    input_tsv: str | None = None
    parent_prediction_id: str | None = None
    catalytic_reference_sites: dict[str, int] | None = None


class ProgressEvent(BaseModel):
    timestamp: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))
    stage: str
    state: EventState
    message: str
    progress: float | None = Field(default=None, ge=0.0, le=1.0)
    artifact: str | None = None


class CandidateRecord(BaseModel):
    model_config = ConfigDict(extra="allow")

    target: str
    source_db: str = ""
    search_engine: str = ""
    qtmscore: float | None = None
    ttmscore: float | None = None
    qcov: float | None = None
    tcov: float | None = None
    rmsd: float | None = None
    wase_screen_score: float | None = None
    wase_tier: str = ""


class RunManifest(BaseModel):
    run_id: str
    project_version: str = "0.4.0"
    status: RunStatus = "created"
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))
    started_at: datetime | None = None
    completed_at: datetime | None = None
    config_path: str
    config_sha256: str = ""
    config_snapshot: dict[str, Any] = Field(default_factory=dict)
    query_path: str = ""
    query_sha256: str = ""
    query_pdb_id: str = ""
    query_chain: str = ""
    query_residue_range: str = ""
    query_source_url: str = ""
    parent_prediction_id: str = ""
    catalytic_reference_sites: dict[str, int] = Field(default_factory=dict)
    software_versions: dict[str, str] = Field(default_factory=dict)
    database_versions: dict[str, str] = Field(default_factory=dict)
    artifacts: dict[str, str] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)
    error: str = ""
    raw_hit_count: int = 0
    filtered_hit_count: int = 0
    annotated_hit_count: int = 0
    child_nr_run_ids: list[str] = Field(default_factory=list)


class RunSummary(BaseModel):
    run_id: str
    status: RunStatus
    run_dir: Path
    manifest_path: Path
    artifacts: dict[str, Path] = Field(default_factory=dict)
    raw_hit_count: int = 0
    filtered_hit_count: int = 0
    annotated_hit_count: int = 0
    warnings: list[str] = Field(default_factory=list)

    model_config = ConfigDict(arbitrary_types_allowed=True)


class PredictionRequest(BaseModel):
    """Validated monomer-only AlphaFold2 prediction request."""

    fasta_path: str | None = None
    fasta_text: str | None = None
    candidate_run_id: str | None = None
    candidate_nr_run_id: str | None = None
    accession: str | None = None
    preset: PredictionPreset = "standard"
    use_templates: bool = False
    confirm_public_msa: bool = False
    prediction_id: str | None = None
    output_root: str = "results/predictions"
    msa_cache_dir: str = ".cache/colabfold_msa"
    candidate_run_root: str = "results/runs"
    candidate_nr_run_root: str = "results/nr_runs"
    wsl_distribution: str = "Ubuntu-Foldseek"
    colabfold_executable: str = "/opt/colabfold/bin/colabfold_batch"
    runtime_backend: Literal["wsl", "native"] = "wsl"
    catalytic_override: dict[str, int] | None = None

    @model_validator(mode="after")
    def validate_input_source(self) -> "PredictionRequest":
        direct_sources = sum(bool(value) for value in (self.fasta_path, self.fasta_text))
        structure_candidate = bool(self.candidate_run_id)
        nr_candidate = bool(self.candidate_nr_run_id)
        if structure_candidate and nr_candidate:
            raise ValueError("candidate_run_id and candidate_nr_run_id are mutually exclusive")
        if (structure_candidate or nr_candidate) and not self.accession:
            raise ValueError("accession is required with a candidate run")
        if self.accession and not (structure_candidate or nr_candidate):
            raise ValueError("accession requires candidate_run_id or candidate_nr_run_id")
        if direct_sources + int(structure_candidate or nr_candidate) != 1:
            raise ValueError("provide exactly one FASTA source or one candidate run/accession pair")
        return self

    @field_validator("catalytic_override")
    @classmethod
    def validate_catalytic_override(cls, value: dict[str, int] | None) -> dict[str, int] | None:
        if value is None:
            return None
        expected = {"S111", "D241", "H274"}
        unknown = set(value) - expected
        if unknown:
            raise ValueError(f"unknown catalytic override labels: {', '.join(sorted(unknown))}")
        if any(position < 1 for position in value.values()):
            raise ValueError("catalytic override positions must be positive")
        return value


class PredictionModelRecord(BaseModel):
    rank: int = Field(ge=1)
    model_name: str
    pdb_path: str
    scores_path: str
    mean_plddt: float | None = Field(default=None, ge=0.0, le=100.0)
    ptm: float | None = Field(default=None, ge=0.0, le=1.0)
    max_pae: float | None = Field(default=None, ge=0.0)
    relaxed: bool = False


class PredictionManifest(BaseModel):
    prediction_id: str
    project_version: str = "0.4.0"
    status: RunStatus = "created"
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))
    started_at: datetime | None = None
    completed_at: datetime | None = None
    sequence_name: str = ""
    sequence_length: int = 0
    sequence_sha256: str = ""
    input_fasta_path: str = ""
    source_type: str = "fasta"
    source_run_id: str = ""
    source_nr_run_id: str = ""
    source_accession: str = ""
    preset: PredictionPreset = "standard"
    parameters: dict[str, Any] = Field(default_factory=dict)
    use_templates: bool = False
    public_msa_confirmed: bool = False
    public_msa_service: str = "https://api.colabfold.com"
    msa_cache_key: str = ""
    msa_cache_hit: bool = False
    msa_path: str = ""
    runtime_versions: dict[str, str] = Field(default_factory=dict)
    runtime_patches: list[dict[str, Any]] = Field(default_factory=list)
    model_parameter_source: str = "https://storage.googleapis.com/alphafold/alphafold_params_2021-07-14.tar"
    model_parameter_provenance: dict[str, Any] = Field(default_factory=dict)
    model_parameter_checksums: dict[str, str] = Field(default_factory=dict)
    models: list[PredictionModelRecord] = Field(default_factory=list)
    rank_001_path: str = ""
    mean_plddt: float | None = None
    ptm: float | None = None
    catalytic_mapping: dict[str, Any] = Field(default_factory=dict)
    child_run_ids: list[str] = Field(default_factory=list)
    artifacts: dict[str, str] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)
    error: str = ""


class PredictionSummary(BaseModel):
    prediction_id: str
    status: RunStatus
    prediction_dir: Path
    manifest_path: Path
    rank_001_path: Path | None = None
    model_count: int = 0
    mean_plddt: float | None = None
    ptm: float | None = None
    child_run_ids: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)

    model_config = ConfigDict(arbitrary_types_allowed=True)


class NrExpansionRequest(BaseModel):
    """Validated request for sequence-space expansion from a structure run."""

    source_run_id: str
    config_path: str = "configs/nr_local_demo.yaml"
    nr_run_id: str | None = None
    output_root: str = "results/nr_runs"
    structure_run_root: str = "results/runs"
    database_path: str | None = None
    database_scope: NrDatabaseScope | None = None
    execution_mode: NrExecutionMode = "local"
    scheduler: Literal["shell", "slurm", "pbs"] = "shell"

    @field_validator("source_run_id", "nr_run_id")
    @classmethod
    def validate_run_ids(cls, value: str | None, info: Any) -> str | None:
        return _validate_run_identifier(value, info.field_name)


class NrCandidateRecord(BaseModel):
    """Sequence-only candidate evidence prior to AF2/structure validation."""

    model_config = ConfigDict(extra="allow")

    target: str
    target_accession: str = ""
    source_seed: str = ""
    target_sequence: str = ""
    evalue: float | None = None
    bits: float | None = None
    fident: float | None = None
    qcov: float | None = None
    tcov: float | None = None
    motif_gxsxg_like: bool = False
    mapped_catalytic_identity_count: int = Field(default=0, ge=0, le=3)
    catalytic_mapping_status: str = ""
    sequence_screen_status: str = ""
    af2_ready: bool = False
    wase_tier: str = "Not evaluated - structure required"
    wase_scope: str = "not_evaluated_sequence_only"


class NrRunManifest(BaseModel):
    nr_run_id: str
    project_version: str = "0.4.0"
    status: RunStatus = "created"
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))
    started_at: datetime | None = None
    completed_at: datetime | None = None
    parent_structure_run_id: str
    execution_mode: NrExecutionMode = "local"
    database_scope: NrDatabaseScope
    database_name: str = "NCBI nr"
    database_snapshot: dict[str, Any] = Field(default_factory=dict)
    config_path: str = ""
    config_sha256: str = ""
    config_snapshot: dict[str, Any] = Field(default_factory=dict)
    paper_preset: str = "hong2025"
    resolved_parameters: dict[str, Any] = Field(default_factory=dict)
    software_versions: dict[str, str] = Field(default_factory=dict)
    seed_count: int = 0
    raw_hit_count: int = 0
    deduplicated_hit_count: int = 0
    representative_count: int = 0
    candidate_count: int = 0
    af2_ready_count: int = 0
    child_prediction_ids: list[str] = Field(default_factory=list)
    artifacts: dict[str, str] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)
    error: str = ""


class NrRunSummary(BaseModel):
    nr_run_id: str
    status: RunStatus
    run_dir: Path
    manifest_path: Path
    artifacts: dict[str, Path] = Field(default_factory=dict)
    candidate_count: int = 0
    af2_ready_count: int = 0
    warnings: list[str] = Field(default_factory=list)

    model_config = ConfigDict(arbitrary_types_allowed=True)


class NrNoveltyRequest(BaseModel):
    fasta_path: str | None = None
    fasta_text: str | None = None
    novelty_id: str | None = None
    output_root: str = "results/nr_novelty"
    backend: Literal["ncbi_remote", "mmseqs"] = "ncbi_remote"
    database_path: str | None = None
    database_scope: NrDatabaseScope = "nr_remote_snapshot"
    max_hits: int = Field(default=20, ge=1, le=5000)
    evalue: float = Field(default=1e-5, gt=0.0)
    confirm_ncbi: bool = False
    contact_email: str | None = None
    runtime_backend: Literal["auto", "native", "wsl"] = "auto"
    mmseqs_executable: str = "mmseqs"
    wsl_distribution: str = "Ubuntu-Foldseek"
    wsl_executable: str = "/opt/colabfold/env/bin/mmseqs"

    @field_validator("novelty_id")
    @classmethod
    def validate_novelty_id(cls, value: str | None) -> str | None:
        return _validate_run_identifier(value, "novelty_id")

    @model_validator(mode="after")
    def validate_novelty_input(self) -> "NrNoveltyRequest":
        if sum(bool(value) for value in (self.fasta_path, self.fasta_text)) != 1:
            raise ValueError("provide exactly one FASTA source")
        if self.backend == "ncbi_remote" and not self.confirm_ncbi:
            raise ValueError("confirm_ncbi is required for the NCBI remote backend")
        if self.backend == "mmseqs" and not self.database_path:
            raise ValueError("database_path is required for the MMseqs backend")
        return self


class NrNoveltyManifest(BaseModel):
    novelty_id: str
    project_version: str = "0.4.0"
    status: RunStatus = "created"
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))
    started_at: datetime | None = None
    completed_at: datetime | None = None
    backend: str
    database_scope: NrDatabaseScope
    sequence_name: str = ""
    sequence_sha256: str = ""
    parameters: dict[str, Any] = Field(default_factory=dict)
    database_snapshot: dict[str, Any] = Field(default_factory=dict)
    software_versions: dict[str, str] = Field(default_factory=dict)
    nearest_nr_accession: str = ""
    nearest_nr_identity: float | None = None
    nearest_nr_coverage: float | None = None
    nearest_nr_evalue: float | None = None
    nearest_nr_taxid: str = ""
    nearest_nr_title: str = ""
    nearest_nr_organism: str = ""
    self_hit_warning: str = ""
    artifacts: dict[str, str] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)
    error: str = ""
