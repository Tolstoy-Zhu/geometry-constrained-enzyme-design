"""Validated public data contracts for AquaEsterNet-StructMine."""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


RunStatus = Literal["created", "running", "completed", "failed", "cancelled"]
EventState = Literal["started", "running", "completed", "warning", "failed", "cancelled"]
PredictionPreset = Literal["quick", "standard", "high"]


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
    project_version: str = "0.3.0"
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
    accession: str | None = None
    preset: PredictionPreset = "standard"
    use_templates: bool = False
    confirm_public_msa: bool = False
    prediction_id: str | None = None
    output_root: str = "results/predictions"
    msa_cache_dir: str = ".cache/colabfold_msa"
    candidate_run_root: str = "results/runs"
    wsl_distribution: str = "Ubuntu-Foldseek"
    colabfold_executable: str = "/opt/colabfold/bin/colabfold_batch"
    runtime_backend: Literal["wsl", "native"] = "wsl"
    catalytic_override: dict[str, int] | None = None

    @model_validator(mode="after")
    def validate_input_source(self) -> "PredictionRequest":
        direct_sources = sum(bool(value) for value in (self.fasta_path, self.fasta_text))
        candidate_requested = bool(self.candidate_run_id or self.accession)
        if candidate_requested and not (self.candidate_run_id and self.accession):
            raise ValueError("candidate_run_id and accession must be provided together")
        if direct_sources + int(candidate_requested) != 1:
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
    project_version: str = "0.3.0"
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
