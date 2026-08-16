from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .tree import Node, postorder, preorder


INF = 10**12
PIRATE_METADATA_COLUMNS = 20


@dataclass(frozen=True)
class Transition:
    branch_id: str
    parent_state: int
    child_state: int

    @property
    def event(self) -> str:
        return "gain" if self.parent_state == 0 else "loss"


class Reconstruction:
    def __init__(self, root: Node, gain_cost: float = 1.0, loss_cost: float = 1.0):
        self.root = root
        self.gain_cost = gain_cost
        self.loss_cost = loss_cost
        self.nodes = preorder(root)
        self.post_nodes = postorder(root)
        self.index = {id(node): i for i, node in enumerate(self.nodes)}

    def infer(
        self,
        leaf_states: dict[str, int | None | set[int] | frozenset[int]],
        root_state: str | int = "auto",
    ) -> tuple[list[Transition], int, float]:
        transitions, chosen_root, score, _ = self.infer_states(leaf_states, root_state)
        return transitions, chosen_root, score

    def infer_states(
        self,
        leaf_states: dict[str, int | None | set[int] | frozenset[int]],
        root_state: str | int = "auto",
    ) -> tuple[list[Transition], int, float, dict[str, int]]:
        costs = [[0.0, 0.0] for _ in self.nodes]
        for node in self.post_nodes:
            idx = self.index[id(node)]
            if node.is_leaf:
                state = leaf_states.get(node.label)
                if state is None or state == {0, 1} or state == frozenset({0, 1}):
                    costs[idx] = [0.0, 0.0]
                elif state == 0 or state == {0} or state == frozenset({0}):
                    costs[idx] = [0.0, INF]
                elif state == 1 or state == {1} or state == frozenset({1}):
                    costs[idx] = [INF, 0.0]
                else:
                    raise ValueError(f"Invalid state for {node.label!r}: {state!r}")
                continue
            for parent_state in (0, 1):
                total = 0.0
                for child in node.children:
                    child_cost = costs[self.index[id(child)]]
                    total += min(
                        child_cost[0] + self._transition_cost(parent_state, 0),
                        child_cost[1] + self._transition_cost(parent_state, 1),
                    )
                costs[idx][parent_state] = total

        root_cost = costs[self.index[id(self.root)]]
        if root_state == "auto":
            chosen_root = 0 if root_cost[0] <= root_cost[1] else 1
        else:
            chosen_root = int(root_state)

        transitions: list[Transition] = []
        inferred_states: dict[str, int] = {}

        def choose(node: Node, state: int) -> None:
            inferred_states[node.node_id] = state
            for child in node.children:
                child_cost = costs[self.index[id(child)]]
                choices = [
                    child_cost[0] + self._transition_cost(state, 0),
                    child_cost[1] + self._transition_cost(state, 1),
                ]
                if choices[0] == choices[1]:
                    child_state = state
                else:
                    child_state = 0 if choices[0] < choices[1] else 1
                if child_state != state:
                    transitions.append(
                        Transition(child.branch_id or "", state, child_state)
                    )
                choose(child, child_state)

        choose(self.root, chosen_root)
        return transitions, chosen_root, root_cost[chosen_root], inferred_states

    def _transition_cost(self, parent: int, child: int) -> float:
        if parent == child:
            return 0.0
        return self.gain_cost if parent == 0 else self.loss_cost


TRUE_STATES = {"1", "+", "true", "yes", "present", "positive", "pos"}
FALSE_STATES = {"0", "-", "false", "no", "absent", "negative", "neg"}
UNKNOWN_STATES = {"", "?", "na", "n/a", "nan", "unknown", "missing", "."}


def parse_state(value: str) -> int | None:
    normalized = value.strip().lower()
    if normalized in TRUE_STATES:
        return 1
    if normalized in FALSE_STATES:
        return 0
    if normalized in UNKNOWN_STATES:
        return None
    raise ValueError(
        f"Unsupported phenotype state {value!r}; use +, -, ?, 1, 0, present, absent, or unknown"
    )


def read_phenotypes(
    path: str | Path,
    selected: Iterable[str] | None = None,
) -> tuple[list[str], dict[str, dict[str, int | None]]]:
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "species_id" not in reader.fieldnames:
            raise ValueError("Phenotype TSV must contain a 'species_id' column")
        available = [name for name in reader.fieldnames if name != "species_id"]
        requested = list(selected or available)
        missing = sorted(set(requested) - set(available))
        if missing:
            raise ValueError(f"Phenotype columns not found: {', '.join(missing)}")
        result = {name: {} for name in requested}
        for row_number, row in enumerate(reader, start=2):
            species = (row.get("species_id") or "").strip()
            if not species:
                raise ValueError(f"Missing species_id in phenotype TSV row {row_number}")
            for name in requested:
                result[name][species] = parse_state(row.get(name) or "")
    return requested, result


def find_orthofinder_files(directory: str | Path) -> dict[str, Path | None]:
    root = Path(directory)
    tree_candidates = [
        root / "Species_Tree" / "SpeciesTree_rooted_node_labels.txt",
        root / "Species_Tree" / "SpeciesTree_rooted.txt",
    ]
    tree = next((path for path in tree_candidates if path.is_file()), None)
    counts = root / "Orthogroups" / "Orthogroups.GeneCount.tsv"
    members = root / "Orthogroups" / "Orthogroups.tsv"
    return {
        "root": root,
        "tree": tree,
        "counts": counts if counts.is_file() else None,
        "members": members if members.is_file() else None,
    }


def find_pirate_files(directory: str | Path) -> dict[str, Path | None]:
    root = Path(directory)
    families = root / "PIRATE.gene_families.tsv"
    tree = root / "binary_presence_absence.nwk"
    representatives = root / "representative_sequences.faa"
    return {
        "root": root,
        "tree": tree if tree.is_file() else None,
        "counts": families if families.is_file() else None,
        "members": families if families.is_file() else None,
        "representatives": representatives if representatives.is_file() else None,
    }


def read_count_header(path: str | Path) -> tuple[list[str], str]:
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    if len(header) < 3:
        raise ValueError("Orthogroups.GeneCount.tsv has too few columns")
    first = header[0]
    species = header[1:]
    if species[-1].strip().lower() == "total":
        species = species[:-1]
    return species, first


def iter_gene_counts(path: str | Path) -> Iterable[tuple[str, list[int]]]:
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        has_total = header[-1].strip().lower() == "total"
        expected = len(header)
        for row_number, row in enumerate(reader, start=2):
            if not row or not any(cell.strip() for cell in row):
                continue
            if len(row) != expected:
                raise ValueError(
                    f"Gene-count row {row_number} has {len(row)} columns; expected {expected}"
                )
            values = row[1:-1] if has_total else row[1:]
            try:
                counts = [int(value or 0) for value in values]
            except ValueError as exc:
                raise ValueError(f"Non-integer gene count in row {row_number}") from exc
            yield row[0], counts


def read_pirate_header(path: str | Path) -> tuple[list[str], str]:
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    if len(header) <= PIRATE_METADATA_COLUMNS:
        raise ValueError("PIRATE.gene_families.tsv has too few columns")
    expected = ["allele_name", "gene_family", "consensus_gene_name", "consensus_product"]
    if [value.strip() for value in header[:4]] != expected:
        raise ValueError("PIRATE.gene_families.tsv has an unsupported header")
    return header[PIRATE_METADATA_COLUMNS:], "gene_family"


def split_pirate_loci(cell: str) -> list[str]:
    """Split PIRATE copies while preserving parenthesized fission loci as one entry."""
    return [value.strip() for value in cell.split(";") if value.strip()]


def iter_pirate_counts(path: str | Path) -> Iterable[tuple[str, list[int]]]:
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        expected = len(header)
        for row_number, row in enumerate(reader, start=2):
            if not row or not any(cell.strip() for cell in row):
                continue
            if len(row) != expected:
                raise ValueError(
                    f"PIRATE row {row_number} has {len(row)} columns; expected {expected}"
                )
            yield row[1], [
                len(split_pirate_loci(cell)) for cell in row[PIRATE_METADATA_COLUMNS:]
            ]
