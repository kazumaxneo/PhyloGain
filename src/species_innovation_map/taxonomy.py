from __future__ import annotations

import csv
import re
from pathlib import Path


RANK_PREFIXES = {
    "d": "domain",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}
RANKS = tuple(RANK_PREFIXES.values())
ID_COLUMNS = ("species_id", "user_genome", "genome_id", "assembly")
TAXONOMY_COLUMNS = ("gtdb_taxonomy", "classification", "taxonomy")


def parse_gtdb_taxonomy(value: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for token in value.split(";"):
        token = token.strip()
        match = re.fullmatch(r"([dpcofgs])__(.*)", token)
        if match and match.group(2).strip():
            result[RANK_PREFIXES[match.group(1)]] = match.group(2).strip()
    return result


def identifier_candidates(value: str) -> list[str]:
    text = value.strip()
    candidates: list[str] = []

    def add(candidate: str) -> None:
        if candidate and candidate not in candidates:
            candidates.append(candidate)

    add(text)
    if "__" in text:
        add(text.split("__", 1)[1])
    for candidate in list(candidates):
        cleaned = re.sub(r"_GTDB$", "", candidate)
        cleaned = re.sub(r"^GTDB_R\d+_", "", cleaned)
        cleaned = re.sub(r"\.(?:fa|faa|fna|fasta|fas)$", "", cleaned, flags=re.I)
        add(cleaned)
        accession = re.search(r"GC[AF]_\d+(?:\.\d+)?", cleaned)
        if accession:
            accession_id = accession.group(0)
            add(accession_id)
            add(accession_id.split(".", 1)[0])
    return candidates


def read_gtdb_taxonomy(
    path: str | Path,
    species_ids: list[str],
) -> tuple[dict[str, dict[str, str]], dict[str, object]]:
    species_lookup: dict[str, set[str]] = {}
    for species_id in species_ids:
        for candidate in identifier_candidates(species_id):
            species_lookup.setdefault(candidate, set()).add(species_id)

    mapped: dict[str, dict[str, str]] = {}
    unmatched_rows = 0
    ambiguous_rows = 0
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        id_column = next((column for column in ID_COLUMNS if column in fields), None)
        taxonomy_column = next(
            (column for column in TAXONOMY_COLUMNS if column in fields), None
        )
        if not id_column or not taxonomy_column:
            raise ValueError(
                "GTDB taxonomy TSV needs an ID column "
                f"({', '.join(ID_COLUMNS)}) and a taxonomy column "
                f"({', '.join(TAXONOMY_COLUMNS)})"
            )
        for row in reader:
            row_id = (row.get(id_column) or "").strip()
            taxonomy = parse_gtdb_taxonomy(row.get(taxonomy_column) or "")
            if not row_id or not taxonomy:
                continue
            matches: set[str] = set()
            for candidate in identifier_candidates(row_id):
                matches.update(species_lookup.get(candidate, ()))
            if len(matches) == 1:
                mapped[next(iter(matches))] = taxonomy
            elif matches:
                ambiguous_rows += 1
            else:
                unmatched_rows += 1

    ranks = [rank for rank in RANKS if any(rank in values for values in mapped.values())]
    return mapped, {
        "mapped_species": len(mapped),
        "unmapped_species": len(species_ids) - len(mapped),
        "unmatched_rows": unmatched_rows,
        "ambiguous_rows": ambiguous_rows,
        "ranks": ranks,
    }
