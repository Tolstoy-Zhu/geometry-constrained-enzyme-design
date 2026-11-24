"""Versioned contracts; scientific labels remain distinct from predictions."""
from __future__ import annotations
import hashlib
import re
from typing import Literal, Any
from pydantic import BaseModel, ConfigDict, Field, model_validator, field_validator

AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
IDENTIFIER = r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$"

class Contract(BaseModel):
    model_config = ConfigDict(extra="forbid", allow_inf_nan=False)
    schema_version: Literal['1.0.0'] = "1.0.0"

class CatalyticSite(Contract):
    role: str
    residue_number: int = Field(gt=0)
    amino_acid: str = Field(pattern=r"^[ACDEFGHIKLMNPQRSTVWY]$")
    atoms: list[str] = Field(min_length=1)

class GeometryRule(Contract):
    first_role: str
    second_role: str
    minimum_A: float = Field(gt=0)
    maximum_A: float = Field(gt=0)

    @model_validator(mode="after")
    def bounds(self):
        if self.maximum_A <= self.minimum_A:
            raise ValueError("Geometry upper limit must exceed lower limit")
        return self

class EnzymeFamilyProfile(Contract):
    family_id: str = Field(pattern=IDENTIFIER)
    name: str
    reference_pdb: str
    reference_chain: str = "A"
    reference_sha256: str
    sites: list[CatalyticSite] = Field(min_length=1)
    geometry: list[GeometryRule] = Field(default_factory=list)
    protected_residues: list[int] = Field(default_factory=list)
    protect_serine_motif: bool = False
    min_tm_score: float = Field(default=0.5, ge=0, le=1)
    contact_protection_A: float = Field(default=4.0, gt=0)
    conservation_fraction: float = Field(default=0.85, ge=0, le=1)
    conservation_coverage: float = Field(default=0.7, ge=0, le=1)
    min_homologs: int = Field(default=3, ge=1)
    assembly: str
    cofactors: list[str] = Field(default_factory=list)
    wase_screen: bool = False
    validation_scope: str
    sources: list[str]

    @model_validator(mode="after")
    def distinct_roles(self):
        roles = [s.role for s in self.sites]
        if len(roles) != len(set(roles)):
            raise ValueError("Catalytic roles must be unique")
        for rule in self.geometry:
            if rule.first_role not in roles or rule.second_role not in roles:
                raise ValueError("Geometry refers to an undefined role")
        return self

class EnzymeRecord(Contract):
    enzyme_id: str = Field(pattern=IDENTIFIER)
    sequence: str
    sequence_sha256: str = ""
    name: str
    family_id: str = "bur01"
    source: Literal["natural", "designed", "experimental", "public_reference"]
    parent_id: str | None = None
    construct_id: str
    construct_description: str = ""
    structure_path: str | None = None
    chain: str | None = None
    numbering_map: list[dict[str, Any]] = Field(default_factory=list)
    aliases: list[str] = Field(default_factory=list)
    provenance: dict[str, Any] = Field(default_factory=dict)

    @model_validator(mode="after")
    def identity(self):
        if not self.sequence or set(self.sequence) - AA:
            raise ValueError("Use a nonempty sequence of standard amino acids; no silent truncation")
        digest = hashlib.sha256(self.sequence.encode("ascii")).hexdigest()
        if self.sequence_sha256 and self.sequence_sha256 != digest:
            raise ValueError("Sequence hash mismatch")
        self.sequence_sha256 = digest
        return self

class Molecule(Contract):
    name: str
    smiles: str
    concentration_M: float | None = Field(default=None, ge=0)
    role: Literal["substrate", "product", "solvent", "buffer", "cofactor"] = "substrate"

class ReactionSpec(Contract):
    reaction_id: str = Field(pattern=IDENTIFIER)
    name: str
    molecules: list[Molecule]
    temperature_K: float | None = Field(default=None, gt=0)
    pH: float | None = Field(default=None, ge=0, le=14)
    water_activity: float | None = Field(default=None, ge=0, le=1)
    buffer: str | None = None
    enzyme_state: Literal["free", "acylated"] = "free"
    provenance: dict[str, Any] = Field(default_factory=dict)

class AssayRecord(Contract):
    assay_id: str = Field(pattern=IDENTIFIER)
    enzyme_id: str
    reaction_id: str
    endpoint: Literal["kcat", "Km", "Tm", "solubility", "synthesis_rate", "hydrolysis_rate",
                      "synthesis_hydrolysis_ratio", "relative_activity", "residual_activity"]
    value: float | None = None
    unit: str
    method: str
    replicate_id: str
    measurement_group: str
    qualifier: Literal["=", "<", ">", "not_measured"] = "="
    detection_limit: float | None = None
    data_kind: Literal["experimental", "simulated", "prediction"] = "experimental"
    source: str
    conditions: dict[str, Any] = Field(default_factory=dict)

    @model_validator(mode="after")
    def observed(self):
        if self.qualifier == "not_measured" and self.value is not None:
            raise ValueError("Unmeasured assays must not contain a measured value")
        if self.qualifier == "=" and self.value is None:
            raise ValueError("An exact observation requires a value")
        if self.qualifier in {"<", ">"} and self.detection_limit is None:
            raise ValueError("Censored observations require a detection limit")
        return self

class TaskSpec(Contract):
    task_id: str = Field(pattern=IDENTIFIER)
    project_id: str = Field(pattern=IDENTIFIER)
    kind: Literal["properties", "design", "predict", "mine", "nr", "dataset", "train", "mechanism", "index", "features", "kinetics"]
    family_id: str = "bur01"
    executor: str = "local"
    inputs: dict[str, Any] = Field(default_factory=dict)
    resource: Literal["gpu", "cpu"] = "gpu"
    parent_task_id: str | None = None
    reaction: ReactionSpec | None = None

class PredictionResult(Contract):
    enzyme_id: str
    reaction_id: str | None = None
    endpoint: str
    unit: str
    value: float | None
    model_id: str
    model_version: str
    applicability: dict[str, Any]
    status: Literal["predicted", "not_trained", "unsupported", "out_of_domain", "failed"]
    interval: tuple[float, float] | None = None
    uncertainty_status: str = "not_calibrated"
