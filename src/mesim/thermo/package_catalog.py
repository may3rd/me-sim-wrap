"""Frozen inventory of DWSIM built-in thermodynamic property packages."""

from dataclasses import dataclass
import json
from pathlib import Path

from ..errors import ValidationError


_STATUSES = frozenset({"pending", "partial", "implemented"})


@dataclass(frozen=True, slots=True)
class DWSIMPackageSource:
    product: str
    version: str
    registry_file: str
    registry_method: str


@dataclass(frozen=True, slots=True)
class ThermodynamicPackageRecord:
    id: str
    dwsim_name: str
    dwsim_class: str
    family: str
    source_file: str
    extraction_status: str
    mesim_model_ids: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ThermodynamicPackageCatalog:
    catalog_id: str
    source: DWSIMPackageSource
    packages: tuple[ThermodynamicPackageRecord, ...]

    def package(self, package_id: str) -> ThermodynamicPackageRecord:
        try:
            return next(record for record in self.packages if record.id == package_id)
        except StopIteration as error:
            raise ValidationError(f"unknown thermodynamic package ID: {package_id}") from error


def load_thermodynamic_package_catalog(
    path: str | Path,
) -> ThermodynamicPackageCatalog:
    """Load and strictly validate one immutable DWSIM package inventory."""
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8-sig"))
        if data["schema_version"] != "dwsim-property-package-catalog-1":
            raise ValidationError("unsupported thermodynamic package catalog schema")
        source = DWSIMPackageSource(**data["source"])
        packages = tuple(
            ThermodynamicPackageRecord(
                id=record["id"],
                dwsim_name=record["dwsim_name"],
                dwsim_class=record["dwsim_class"],
                family=record["family"],
                source_file=record["source_file"],
                extraction_status=record["extraction_status"],
                mesim_model_ids=tuple(record["mesim_model_ids"]),
            )
            for record in data["packages"]
        )
        catalog = ThermodynamicPackageCatalog(data["catalog_id"], source, packages)
    except ValidationError:
        raise
    except (KeyError, TypeError, ValueError) as error:
        raise ValidationError(f"invalid thermodynamic package catalog: {error}") from error

    source_values = (
        source.product,
        source.version,
        source.registry_file,
        source.registry_method,
    )
    if not catalog.catalog_id or any(not isinstance(value, str) or not value for value in source_values):
        raise ValidationError("thermodynamic package catalog identity must be non-empty")
    if not packages:
        raise ValidationError("thermodynamic package catalog must not be empty")

    for attribute in ("id", "dwsim_name", "dwsim_class"):
        values = tuple(getattr(record, attribute) for record in packages)
        if any(not isinstance(value, str) or not value for value in values):
            raise ValidationError(f"thermodynamic package {attribute} must be non-empty")
        if len(set(values)) != len(values):
            raise ValidationError(f"thermodynamic package {attribute} values must be unique")
    model_ids: list[str] = []
    for record in packages:
        if (
            not record.family
            or not record.source_file
            or record.extraction_status not in _STATUSES
            or any(not isinstance(value, str) or not value for value in record.mesim_model_ids)
            or len(set(record.mesim_model_ids)) != len(record.mesim_model_ids)
        ):
            raise ValidationError(f"invalid thermodynamic package record: {record.id}")
        if record.extraction_status == "pending" and record.mesim_model_ids:
            raise ValidationError(f"pending package cannot expose model IDs: {record.id}")
        if record.extraction_status != "pending" and not record.mesim_model_ids:
            raise ValidationError(f"extracted package must expose a model ID: {record.id}")
        model_ids.extend(record.mesim_model_ids)
    if len(set(model_ids)) != len(model_ids):
        raise ValidationError("thermodynamic package model IDs must be globally unique")
    return catalog
