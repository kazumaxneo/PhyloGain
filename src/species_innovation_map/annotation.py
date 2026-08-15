from __future__ import annotations

import csv
import math
import shutil
import sqlite3
import subprocess
import time
from pathlib import Path
from urllib.request import Request, urlopen


COG_CATEGORIES = {
    "J": "Translation, ribosomal structure and biogenesis",
    "A": "RNA processing and modification",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics",
    "D": "Cell cycle control, cell division and chromosome partitioning",
    "Y": "Nuclear structure",
    "V": "Defense mechanisms",
    "T": "Signal transduction mechanisms",
    "M": "Cell wall, membrane and envelope biogenesis",
    "N": "Cell motility",
    "Z": "Cytoskeleton",
    "W": "Extracellular structures",
    "U": "Intracellular trafficking, secretion and vesicular transport",
    "O": "Posttranslational modification, protein turnover and chaperones",
    "X": "Mobilome: prophages and transposons",
    "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolite biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
}

TERM_COLUMNS = {
    "GOs": "GO",
    "KEGG_ko": "KEGG KO",
    "KEGG_Pathway": "KEGG pathway",
    "KEGG_Module": "KEGG module",
    "KEGG_Reaction": "KEGG reaction",
    "PFAMs": "Pfam",
}

KEGG_LIST_DATABASES = {
    "KEGG KO": "ko",
    "KEGG pathway": "pathway",
    "KEGG module": "module",
    "KEGG reaction": "reaction",
}
KEGG_REST_URL = "https://rest.kegg.jp"


def write_representative_fasta(
    members_path: str | Path,
    proteomes: str | Path,
    output_path: str | Path,
) -> int:
    """Write one protein sequence per orthogroup, choosing its first listed member."""
    wanted: dict[str, str] = {}
    with Path(members_path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        for row in reader:
            if not row:
                continue
            gene = next(
                (
                    cell.split(",", 1)[0].strip()
                    for cell in row[1:]
                    if cell.strip()
                ),
                "",
            )
            if gene:
                wanted[gene] = row[0]

    sequences: dict[str, str] = {}
    proteome_root = Path(proteomes)
    direct_files = [
        path
        for path in proteome_root.iterdir()
        if path.is_file() and path.suffix.lower() in {".fa", ".faa", ".fasta", ".fas"}
    ]
    fasta_files = direct_files or [
        path
        for path in proteome_root.rglob("*")
        if path.is_file() and path.suffix.lower() in {".fa", ".faa", ".fasta", ".fas"}
    ]
    if not fasta_files:
        raise ValueError(f"No protein FASTA files found under: {proteomes}")
    for path in sorted(fasta_files):
        for family, sequence in _read_wanted_fasta(path, wanted):
            if family not in sequences:
                sequences[family] = sequence.rstrip("*")

    missing = len(wanted) - len(sequences)
    if missing:
        raise ValueError(
            f"Could not find representative protein sequences for {missing:,} orthogroups; "
            "check that --proteomes contains the FASTA files used by OrthoFinder"
        )
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="\n") as handle:
        for family in sorted(sequences):
            handle.write(f">{family}\n")
            sequence = sequences[family]
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")
    return len(sequences)


def run_eggnog_mapper(
    fasta_path: str | Path,
    output_directory: str | Path,
    emapper: str = "emapper.py",
    data_dir: str | Path | None = None,
    cpu: int = 1,
) -> Path:
    executable = shutil.which(emapper) if not Path(emapper).is_file() else emapper
    if not executable:
        raise ValueError(
            f"eggNOG-mapper executable was not found: {emapper}; "
            "install eggnog-mapper or pass --eggnog-emapper"
        )
    output_dir = Path(output_directory)
    output_dir.mkdir(parents=True, exist_ok=True)
    command = [
        str(executable),
        "-i",
        str(Path(fasta_path).resolve()),
        "--itype",
        "proteins",
        "--output",
        "eggnog",
        "--output_dir",
        str(output_dir.resolve()),
        "--cpu",
        str(max(cpu, 1)),
        "--override",
    ]
    if data_dir:
        command.extend(["--data_dir", str(Path(data_dir).resolve())])
    subprocess.run(command, check=True)
    result = output_dir / "eggnog.emapper.annotations"
    if not result.is_file():
        raise ValueError("eggNOG-mapper finished without creating eggnog.emapper.annotations")
    return result


def import_eggnog_annotations(
    connection: sqlite3.Connection,
    annotation_path: str | Path,
) -> dict[str, int]:
    valid_families = {row[0] for row in connection.execute("SELECT family_id FROM families")}
    header: list[str] | None = None
    family_rows: list[tuple[str, str, str]] = []
    annotated_families: set[str] = set()
    term_rows: set[tuple[str, str, str]] = set()
    link_rows: set[tuple[str, str, str]] = set()
    with Path(annotation_path).open("r", encoding="utf-8-sig", newline="") as handle:
        for line in handle:
            if line.startswith("##") or not line.strip():
                continue
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\r\n").split("\t")
                continue
            if header is None:
                raise ValueError("eggNOG annotation file has no #query header")
            values = line.rstrip("\r\n").split("\t")
            row = dict(zip(header, values))
            family = row.get("query", "").strip().split(";", 1)[0]
            if family not in valid_families and "_" in family:
                pirate_family = family.rsplit("_", 1)[0]
                if pirate_family in valid_families:
                    family = pirate_family
            if not family or family not in valid_families:
                continue
            preferred = _clean(row.get("Preferred_name", ""))
            description = _clean(row.get("Description", ""))
            family_rows.append((preferred, description, family))
            annotated_families.add(family)
            for source, term_id, term_name in _row_terms(row):
                term_rows.add((source, term_id, term_name))
                link_rows.add((family, source, term_id))

    connection.executemany(
        "UPDATE families SET preferred_name=?,description=? WHERE family_id=?",
        family_rows,
    )
    connection.executemany(
        "INSERT OR IGNORE INTO annotation_terms(source,term_id,term_name) VALUES(?,?,?)",
        sorted(term_rows),
    )
    connection.executemany(
        "INSERT OR IGNORE INTO family_terms(family_id,source,term_id) VALUES(?,?,?)",
        sorted(link_rows),
    )
    _update_kegg_ko_names(connection)
    connection.execute(
        """
        UPDATE annotation_terms
        SET family_count=(
          SELECT COUNT(*) FROM family_terms ft
          WHERE ft.source=annotation_terms.source AND ft.term_id=annotation_terms.term_id
        )
        """
    )
    return {"annotated_families": len(annotated_families), "family_term_links": len(link_rows)}


def import_go_term_names(
    connection: sqlite3.Connection,
    obo_path: str | Path,
) -> int:
    """Replace GO identifiers with official term names from a GO OBO file."""
    names: dict[str, str] = {}
    term_id = ""
    term_name = ""
    alternate_ids: list[str] = []

    def store_term() -> None:
        if term_id and term_name:
            names[term_id] = term_name
            for alternate_id in alternate_ids:
                names[alternate_id] = term_name

    with Path(obo_path).open("r", encoding="utf-8-sig") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line == "[Term]":
                store_term()
                term_id = ""
                term_name = ""
                alternate_ids = []
            elif line.startswith("id: GO:"):
                term_id = line[4:].strip()
            elif line.startswith("name: "):
                term_name = line[6:].strip()
            elif line.startswith("alt_id: GO:"):
                alternate_ids.append(line[8:].strip())
            elif line.startswith("[") and line != "[Term]":
                store_term()
                term_id = ""
                term_name = ""
                alternate_ids = []
        store_term()

    existing = {
        row[0]
        for row in connection.execute(
            "SELECT term_id FROM annotation_terms WHERE source='GO'"
        )
    }
    updates = [(names[term_id], term_id) for term_id in existing if term_id in names]
    connection.executemany(
        "UPDATE annotation_terms SET term_name=? WHERE source='GO' AND term_id=?",
        updates,
    )
    return len(updates)


def fetch_official_kegg_names(
    connection: sqlite3.Connection,
    cache_path: str | Path | None = None,
    base_url: str = KEGG_REST_URL,
) -> dict[str, int]:
    """Fetch official KEGG labels once, cache them, and update stored term names."""
    fetched: dict[str, dict[str, str]] = {}
    for index, (source, database) in enumerate(KEGG_LIST_DATABASES.items()):
        request = Request(
            f"{base_url.rstrip('/')}/list/{database}",
            headers={"User-Agent": "species-innovation-map"},
        )
        with urlopen(request, timeout=120) as response:
            text = response.read().decode("utf-8")
        fetched[source] = _parse_kegg_list(text)
        if index + 1 < len(KEGG_LIST_DATABASES):
            time.sleep(0.4)

    existing = connection.execute(
        "SELECT source,term_id FROM annotation_terms WHERE source LIKE 'KEGG %'"
    ).fetchall()
    updates: list[tuple[str, str, str]] = []
    cached_rows: list[tuple[str, str, str]] = []
    counts = {source: 0 for source in KEGG_LIST_DATABASES}
    for source, term_id in existing:
        lookup_id = term_id
        if source == "KEGG pathway" and term_id.startswith("ko"):
            lookup_id = "map" + term_id[2:]
        term_name = fetched.get(source, {}).get(lookup_id)
        if not term_name:
            continue
        updates.append((term_name, source, term_id))
        cached_rows.append((source, term_id, term_name))
        counts[source] += 1
    unavailable_modules = [
        term_id
        for source, term_id in existing
        if source == "KEGG module" and term_id not in fetched["KEGG module"]
    ]
    retired_label = "Retired KEGG module (not present in the current KEGG release)"
    updates.extend(
        (retired_label, "KEGG module", term_id) for term_id in unavailable_modules
    )
    cached_rows.extend(
        ("KEGG module", term_id, retired_label) for term_id in unavailable_modules
    )
    connection.executemany(
        "UPDATE annotation_terms SET term_name=? WHERE source=? AND term_id=?",
        updates,
    )
    if cache_path:
        destination = Path(cache_path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["source", "term_id", "term_name"])
            writer.writerows(sorted(cached_rows))
    counts["total"] = len(updates)
    counts["retired_modules"] = len(unavailable_modules)
    return counts


def _parse_kegg_list(text: str) -> dict[str, str]:
    names: dict[str, str] = {}
    for line in text.splitlines():
        if "\t" not in line:
            continue
        term_id, term_name = line.split("\t", 1)
        term_id = term_id.removeprefix("ko:")
        if term_id and term_name:
            names[term_id] = term_name
    return names


def _update_kegg_ko_names(connection: sqlite3.Connection) -> int:
    """Use the most common eggNOG family description as a readable KO label."""
    rows = connection.execute(
        """
        SELECT ft.term_id,f.description,COUNT(*) AS frequency
        FROM family_terms ft
        JOIN families f ON f.family_id=ft.family_id
        WHERE ft.source='KEGG KO' AND TRIM(COALESCE(f.description,'')) NOT IN ('','-')
        GROUP BY ft.term_id,f.description
        ORDER BY ft.term_id,frequency DESC,LENGTH(f.description),f.description
        """
    ).fetchall()
    selected: dict[str, str] = {}
    for term_id, description, _frequency in rows:
        selected.setdefault(term_id, description)
    connection.executemany(
        "UPDATE annotation_terms SET term_name=? WHERE source='KEGG KO' AND term_id=?",
        [(description, term_id) for term_id, description in selected.items()],
    )
    return len(selected)


def branch_enrichment(
    connection: sqlite3.Connection,
    branch_id: str,
    event: str,
    limit: int = 20,
    min_overlap: int = 2,
    source: str | None = None,
) -> dict[str, object]:
    if event not in {"gain", "loss"}:
        raise ValueError("event must be gain or loss")
    universe = connection.execute("SELECT COUNT(*) FROM families").fetchone()[0]
    foreground = connection.execute(
        "SELECT COUNT(DISTINCT family_id) FROM events WHERE branch_id=? AND event=?",
        (branch_id, event),
    ).fetchone()[0]
    all_sources = source == "All"
    if source is not None and not all_sources:
        available = {
            row[0]
            for row in connection.execute("SELECT DISTINCT source FROM annotation_terms")
        }
        if source not in available:
            raise ValueError(f"Unknown annotation database: {source}")
    source_filter = " AND ft.source=?" if source is not None and not all_sources else ""
    parameters: list[object] = [branch_id, event]
    if source is not None and not all_sources:
        parameters.append(source)
    rows = connection.execute(
        f"""
        SELECT ft.source,ft.term_id,at.term_name,
               COUNT(DISTINCT ft.family_id) AS overlap,at.family_count
        FROM events e
        JOIN family_terms ft ON ft.family_id=e.family_id
        JOIN annotation_terms at ON at.source=ft.source AND at.term_id=ft.term_id
        WHERE e.branch_id=? AND e.event=?{source_filter}
        GROUP BY ft.source,ft.term_id,at.term_name,at.family_count
        """,
        parameters,
    ).fetchall()
    tested = []
    for row_source, term_id, term_name, overlap, term_total in rows:
        p_value = _hypergeom_sf(overlap, universe, term_total, foreground)
        fold = (overlap / foreground) / (term_total / universe) if foreground else 0.0
        tested.append(
            {
                "source": row_source,
                "term_id": term_id,
                "term_name": term_name,
                "overlap": overlap,
                "foreground": foreground,
                "term_total": term_total,
                "universe": universe,
                "fold_enrichment": fold,
                "p_value": p_value,
            }
        )
    if all_sources:
        sources = {row["source"] for row in tested}
        for row_source in sources:
            _benjamini_hochberg([row for row in tested if row["source"] == row_source])
    else:
        _benjamini_hochberg(tested)
    results = [row for row in tested if row["overlap"] >= min_overlap]
    results.sort(key=lambda row: (row["q_value"], row["p_value"], -row["overlap"]))
    return {
        "branch_id": branch_id,
        "event": event,
        "source": source,
        "foreground": foreground,
        "universe": universe,
        "tested_terms": len(tested),
        "results": results[:limit],
    }


def annotation_sources(connection: sqlite3.Connection) -> list[dict[str, object]]:
    preferred_order = {
        name: index
        for index, name in enumerate(
            [
                "KEGG pathway",
                "KEGG module",
                "KEGG KO",
                "GO",
                "Pfam",
                "COG category",
                "KEGG reaction",
            ]
        )
    }
    rows = connection.execute(
        """
        SELECT at.source,COUNT(DISTINCT at.term_id) AS term_count,
               COUNT(DISTINCT ft.family_id) AS family_count
        FROM annotation_terms at
        LEFT JOIN family_terms ft ON ft.source=at.source AND ft.term_id=at.term_id
        GROUP BY at.source
        """
    ).fetchall()
    results = [
        {"source": row[0], "term_count": row[1], "family_count": row[2]}
        for row in rows
    ]
    results.sort(key=lambda row: (preferred_order.get(row["source"], 999), row["source"]))
    return results


def _read_wanted_fasta(path: Path, wanted: dict[str, str]):
    family = ""
    parts: list[str] = []
    with path.open("r", encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if family:
                    yield family, "".join(parts)
                family = wanted.get(line[1:].split()[0], "")
                parts = []
            elif family:
                parts.append(line)
    if family:
        yield family, "".join(parts)


def _clean(value: str) -> str:
    value = value.strip()
    return "" if value in {"", "-"} else value


def _split_terms(value: str):
    return [item.strip() for item in value.split(",") if item.strip() and item.strip() != "-"]


def _row_terms(row: dict[str, str]):
    for column, source in TERM_COLUMNS.items():
        for term in _split_terms(row.get(column, "")):
            if source == "KEGG KO" and term.startswith("ko:"):
                term = term[3:]
            yield source, term, term
    categories = _clean(row.get("COG_category", ""))
    for category in categories:
        if category not in {"-", " ", "R", "S"}:
            yield "COG category", category, COG_CATEGORIES.get(category, category)


def _hypergeom_sf(observed: int, population: int, successes: int, draws: int) -> float:
    maximum = min(successes, draws)
    if observed <= 0:
        return 1.0
    if observed > maximum or population <= 0:
        return 0.0
    logs = []
    for value in range(observed, maximum + 1):
        if draws - value > population - successes:
            continue
        logs.append(
            _log_comb(successes, value)
            + _log_comb(population - successes, draws - value)
            - _log_comb(population, draws)
        )
    if not logs:
        return 0.0
    peak = max(logs)
    return min(1.0, math.exp(peak) * sum(math.exp(value - peak) for value in logs))


def _log_comb(total: int, selected: int) -> float:
    if selected < 0 or selected > total:
        return float("-inf")
    return math.lgamma(total + 1) - math.lgamma(selected + 1) - math.lgamma(total - selected + 1)


def _benjamini_hochberg(rows: list[dict[str, object]]) -> None:
    ordered = sorted(enumerate(rows), key=lambda item: item[1]["p_value"])
    adjusted = 1.0
    count = len(ordered)
    for rank_from_end in range(count - 1, -1, -1):
        original_index, row = ordered[rank_from_end]
        rank = rank_from_end + 1
        adjusted = min(adjusted, float(row["p_value"]) * count / rank)
        rows[original_index]["q_value"] = min(1.0, adjusted)
