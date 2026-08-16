from __future__ import annotations

from dataclasses import dataclass

from .models import TaxonGroup
from .tree import Node, assign_ids, postorder


@dataclass(frozen=True)
class TaxonCollapseResult:
    root: Node
    groups: dict[str, TaxonGroup]
    excluded_non_monophyletic: tuple[str, ...]
    excluded_small: tuple[str, ...]


def collapse_tree_to_taxa(
    root: Node,
    groups: dict[str, TaxonGroup],
    min_genomes_per_taxon: int = 3,
) -> TaxonCollapseResult:
    if min_genomes_per_taxon < 1:
        raise ValueError("--min-genomes-per-taxon must be at least 1")
    node_leaves: dict[int, frozenset[str]] = {}
    for node in postorder(root):
        node_leaves[id(node)] = (
            frozenset({node.label})
            if node.is_leaf
            else frozenset().union(*(node_leaves[id(child)] for child in node.children))
        )
    clade_by_members: dict[frozenset[str], Node] = {}
    for node in postorder(root):
        leaves = node_leaves[id(node)]
        if leaves:
            clade_by_members[leaves] = node
    accepted: dict[str, TaxonGroup] = {}
    non_monophyletic: list[str] = []
    small: list[str] = []
    replacement: dict[int, TaxonGroup] = {}
    for taxon_id, group in groups.items():
        if group.n_genomes < min_genomes_per_taxon:
            small.append(taxon_id)
            continue
        clade = clade_by_members.get(frozenset(group.member_genomes))
        if clade is None:
            non_monophyletic.append(taxon_id)
            continue
        accepted[taxon_id] = group
        replacement[id(clade)] = group

    def clone(node: Node) -> Node | None:
        group = replacement.get(id(node))
        if group is not None:
            return Node(
                label=group.taxon_id,
                length=node.length,
                support=node.support,
                metadata={
                    "taxon_id": group.taxon_id,
                    "taxon_rank": group.rank,
                    "n_genomes": group.n_genomes,
                    "member_genomes": list(group.member_genomes),
                    "display_label": f"{group.taxon_id} (n={group.n_genomes})",
                },
            )
        children = [child for child in (clone(item) for item in node.children) if child]
        if not children:
            return None
        if len(children) == 1:
            child = children[0]
            if node.length is not None:
                child.length = (child.length or 0.0) + node.length
            return child
        copied = Node(
            label=node.label,
            length=node.length,
            support=node.support,
            children=children,
        )
        for child in children:
            child.parent = copied
        return copied

    collapsed = clone(root)
    if collapsed is None:
        raise ValueError("No taxa remained after taxon-level filtering")
    collapsed.parent = None
    assign_ids(collapsed)
    return TaxonCollapseResult(
        collapsed,
        accepted,
        tuple(sorted(non_monophyletic)),
        tuple(sorted(small)),
    )
