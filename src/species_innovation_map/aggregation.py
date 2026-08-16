from __future__ import annotations

import math
from collections.abc import Iterable

from .models import TaxonGeneState, TaxonGroup


RANK_PREFIXES = {
    "domain": "d",
    "phylum": "p",
    "class": "c",
    "order": "o",
    "family": "f",
    "genus": "g",
    "species": "s",
}


def taxon_identifier(rank: str, value: str) -> str:
    prefix = RANK_PREFIXES[rank]
    return value if value.startswith(f"{prefix}__") else f"{prefix}__{value}"


def group_genomes_by_taxon(
    genome_ids: Iterable[str],
    taxonomy: dict[str, dict[str, str]],
    rank: str,
) -> tuple[dict[str, TaxonGroup], list[str]]:
    grouped: dict[str, list[str]] = {}
    unmapped: list[str] = []
    for genome_id in genome_ids:
        value = taxonomy.get(genome_id, {}).get(rank, "").strip()
        if not value:
            unmapped.append(genome_id)
            continue
        grouped.setdefault(taxon_identifier(rank, value), []).append(genome_id)
    return (
        {
            taxon_id: TaxonGroup(taxon_id, rank, tuple(members))
            for taxon_id, members in grouped.items()
        },
        unmapped,
    )


def classify_occupancy(
    occupancy: float,
    present_threshold: float = 0.90,
    absent_threshold: float = 0.10,
    *,
    present_count: int | None = None,
    total_count: int | None = None,
    state_method: str = "threshold",
    confidence: float = 0.95,
) -> str:
    validate_thresholds(present_threshold, absent_threshold)
    if state_method not in {"threshold", "confidence"}:
        raise ValueError("state_method must be 'threshold' or 'confidence'")
    if not 0 < confidence < 1:
        raise ValueError("confidence must be between 0 and 1")
    if state_method == "confidence":
        if present_count is None or total_count is None or total_count < 1:
            raise ValueError("confidence state method requires present_count and total_count")
        lower_tail, upper_tail = beta_posterior_tail_probabilities(
            present_count, total_count, absent_threshold, present_threshold
        )
        if upper_tail >= confidence:
            return "present"
        if lower_tail >= confidence:
            return "absent"
        return "polymorphic"
    if occupancy >= present_threshold:
        return "present"
    if occupancy <= absent_threshold:
        return "absent"
    return "polymorphic"


def beta_posterior_tail_probabilities(
    present_count: int,
    total_count: int,
    absent_threshold: float,
    present_threshold: float,
) -> tuple[float, float]:
    """Return P(theta <= absent_threshold), P(theta >= present_threshold).

    For Beta(k+1, n-k+1), the integer-parameter beta CDF is evaluated via
    the equivalent binomial tail, avoiding a SciPy dependency.
    """
    k = max(0, min(int(present_count), int(total_count)))
    n = max(1, int(total_count))
    # Beta(k+1, n-k+1) has integer parameters whose sum is n+2;
    # the equivalent binomial distribution therefore has n+1 trials.
    trials = n + 1

    def pmf_sum(start: int, stop: int, p: float) -> float:
        q = 1.0 - p
        return sum(
            math.comb(trials, j) * (p**j) * (q ** (trials - j))
            for j in range(start, stop + 1)
        )

    # I_x(k+1,n-k+1) = P(Binomial(n+1,x) >= k+1).
    lower = pmf_sum(k + 1, trials, absent_threshold) if k < trials else 0.0
    upper = pmf_sum(0, k, present_threshold)
    return min(1.0, max(0.0, lower)), min(1.0, max(0.0, upper))


def validate_thresholds(present_threshold: float, absent_threshold: float) -> None:
    if not 0 <= absent_threshold < present_threshold <= 1:
        raise ValueError(
            "Taxon thresholds must satisfy 0 <= absent < present <= 1"
        )


def aggregate_family_counts(
    gene_family: str,
    genome_ids: list[str],
    counts: list[int],
    groups: dict[str, TaxonGroup],
    present_threshold: float = 0.90,
    absent_threshold: float = 0.10,
    copy_presence_threshold: int = 1,
    state_method: str = "threshold",
    confidence: float = 0.95,
) -> dict[str, TaxonGeneState]:
    validate_thresholds(present_threshold, absent_threshold)
    if len(genome_ids) != len(counts):
        raise ValueError("Genome IDs and count values have different lengths")
    genome_presence = {
        genome_id: count >= copy_presence_threshold
        for genome_id, count in zip(genome_ids, counts)
    }
    states: dict[str, TaxonGeneState] = {}
    for taxon_id, group in groups.items():
        total = group.n_genomes
        present = sum(bool(genome_presence.get(genome)) for genome in group.member_genomes)
        occupancy = present / total if total else 0.0
        states[taxon_id] = TaxonGeneState(
            taxon_id=taxon_id,
            gene_family=gene_family,
            present_count=present,
            total_count=total,
            occupancy=occupancy,
            observed_state=classify_occupancy(
                occupancy, present_threshold, absent_threshold,
                present_count=present, total_count=total,
                state_method=state_method, confidence=confidence,
            ),
        )
    return states


def aggregate_binary_character(
    genome_states: dict[str, int | None],
    groups: dict[str, TaxonGroup],
    present_threshold: float,
    absent_threshold: float,
    state_method: str = "threshold",
    confidence: float = 0.95,
) -> dict[str, frozenset[int]]:
    result: dict[str, frozenset[int]] = {}
    for taxon_id, group in groups.items():
        known = [genome_states.get(genome) for genome in group.member_genomes]
        known = [state for state in known if state in {0, 1}]
        if not known:
            result[taxon_id] = frozenset({0, 1})
            continue
        occupancy = sum(known) / len(known)
        observed = classify_occupancy(
            occupancy, present_threshold, absent_threshold,
            present_count=sum(known), total_count=len(known),
            state_method=state_method, confidence=confidence,
        )
        result[taxon_id] = (
            frozenset({1})
            if observed == "present"
            else frozenset({0})
            if observed == "absent"
            else frozenset({0, 1})
        )
    return result
