from __future__ import annotations

import csv
import math
from pathlib import Path

from .taxonomy import ID_COLUMNS, identifier_candidates


GENOME_SIZE_COLUMNS = {
    "genome_size_bp": 1.0,
    "genome_length_bp": 1.0,
    "total_length": 1.0,
    "genome_size": 1.0,
    "genome_length": 1.0,
    "genome_size_mb": 1_000_000.0,
}

TIP_LABEL_COLUMNS = (
    "strain_name",
    "strain",
    "organism_name",
    "organism",
    "display_name",
    "tip_label",
    "name",
)


def read_genome_sizes(
    path: str | Path,
    species_ids: list[str],
) -> tuple[dict[str, int], dict[str, int | str]]:
    species_lookup: dict[str, set[str]] = {}
    for species_id in species_ids:
        for candidate in identifier_candidates(species_id):
            species_lookup.setdefault(candidate, set()).add(species_id)

    mapped: dict[str, int] = {}
    unmatched_rows = 0
    ambiguous_rows = 0
    invalid_rows = 0
    duplicate_rows = 0
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        id_column = next((column for column in ID_COLUMNS if column in fields), None)
        size_column = next(
            (column for column in GENOME_SIZE_COLUMNS if column in fields), None
        )
        if not id_column or not size_column:
            raise ValueError(
                "Genome metadata TSV needs an ID column "
                f"({', '.join(ID_COLUMNS)}) and a genome-size column "
                f"({', '.join(GENOME_SIZE_COLUMNS)})"
            )
        factor = GENOME_SIZE_COLUMNS[size_column]
        for row in reader:
            row_id = (row.get(id_column) or "").strip()
            raw_size = (row.get(size_column) or "").strip().replace(",", "")
            if not row_id or not raw_size:
                invalid_rows += 1
                continue
            try:
                size_bp = float(raw_size) * factor
            except ValueError:
                invalid_rows += 1
                continue
            if not math.isfinite(size_bp) or size_bp <= 0:
                invalid_rows += 1
                continue
            matches: set[str] = set()
            for candidate in identifier_candidates(row_id):
                matches.update(species_lookup.get(candidate, ()))
            if len(matches) == 1:
                species_id = next(iter(matches))
                if species_id in mapped:
                    duplicate_rows += 1
                    continue
                mapped[species_id] = round(size_bp)
            elif matches:
                ambiguous_rows += 1
            else:
                unmatched_rows += 1

    return mapped, {
        "source_column": size_column,
        "mapped_species": len(mapped),
        "unmapped_species": len(species_ids) - len(mapped),
        "unmatched_rows": unmatched_rows,
        "ambiguous_rows": ambiguous_rows,
        "invalid_rows": invalid_rows,
        "duplicate_rows": duplicate_rows,
    }


def read_tip_labels(
    path: str | Path,
    species_ids: list[str],
) -> tuple[dict[str, str], dict[str, int | str]]:
    species_lookup: dict[str, set[str]] = {}
    for species_id in species_ids:
        for candidate in identifier_candidates(species_id):
            species_lookup.setdefault(candidate, set()).add(species_id)

    mapped: dict[str, str] = {}
    unmatched_rows = 0
    ambiguous_rows = 0
    invalid_rows = 0
    duplicate_rows = 0
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        id_column = next((column for column in ID_COLUMNS if column in fields), None)
        label_column = next(
            (column for column in TIP_LABEL_COLUMNS if column in fields), None
        )
        if not id_column or not label_column:
            raise ValueError(
                "Tip metadata TSV needs an ID column "
                f"({', '.join(ID_COLUMNS)}) and a display-label column "
                f"({', '.join(TIP_LABEL_COLUMNS)})"
            )
        for row in reader:
            row_id = (row.get(id_column) or "").strip()
            label = (row.get(label_column) or "").strip()
            if not row_id or not label:
                invalid_rows += 1
                continue
            matches: set[str] = set()
            for candidate in identifier_candidates(row_id):
                matches.update(species_lookup.get(candidate, ()))
            if len(matches) == 1:
                species_id = next(iter(matches))
                if species_id in mapped:
                    duplicate_rows += 1
                    continue
                mapped[species_id] = label
            elif matches:
                ambiguous_rows += 1
            else:
                unmatched_rows += 1

    return mapped, {
        "source_column": label_column,
        "mapped_species": len(mapped),
        "unmapped_species": len(species_ids) - len(mapped),
        "unmatched_rows": unmatched_rows,
        "ambiguous_rows": ambiguous_rows,
        "invalid_rows": invalid_rows,
        "duplicate_rows": duplicate_rows,
    }
