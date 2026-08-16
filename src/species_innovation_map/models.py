from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class GenomeMetadata:
    genome_id: str
    taxonomy: dict[str, str]
    species: str | None = None
    genus: str | None = None
    family: str | None = None
    order: str | None = None
    completeness: float | None = None
    contamination: float | None = None
    assembly_status: str | None = None


@dataclass(frozen=True)
class TaxonGroup:
    taxon_id: str
    rank: str
    member_genomes: tuple[str, ...]

    @property
    def n_genomes(self) -> int:
        return len(self.member_genomes)


@dataclass(frozen=True)
class TaxonGeneState:
    taxon_id: str
    gene_family: str
    present_count: int
    total_count: int
    occupancy: float
    observed_state: str

    @property
    def sankoff_state(self) -> frozenset[int]:
        if self.observed_state == "present":
            return frozenset({1})
        if self.observed_state == "absent":
            return frozenset({0})
        return frozenset({0, 1})
