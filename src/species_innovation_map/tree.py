from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Node:
    label: str = ""
    length: float | None = None
    children: list["Node"] = field(default_factory=list)
    parent: "Node | None" = None
    node_id: str = ""
    branch_id: str | None = None
    depth: int = 0

    @property
    def is_leaf(self) -> bool:
        return not self.children


class NewickError(ValueError):
    pass


class _Parser:
    def __init__(self, text: str):
        self.text = text.strip()
        self.pos = 0

    def parse(self) -> Node:
        if not self.text:
            raise NewickError("The species tree is empty")
        root = self._subtree()
        self._space()
        if self.pos < len(self.text) and self.text[self.pos] == ";":
            self.pos += 1
        self._space()
        if self.pos != len(self.text):
            raise NewickError(f"Unexpected Newick content at character {self.pos + 1}")
        return root

    def _space(self) -> None:
        while self.pos < len(self.text) and self.text[self.pos].isspace():
            self.pos += 1

    def _subtree(self) -> Node:
        self._space()
        node = Node()
        if self.pos < len(self.text) and self.text[self.pos] == "(":
            self.pos += 1
            while True:
                child = self._subtree()
                child.parent = node
                node.children.append(child)
                self._space()
                if self.pos >= len(self.text):
                    raise NewickError("Unclosed parenthesis in species tree")
                char = self.text[self.pos]
                if char == ",":
                    self.pos += 1
                    continue
                if char == ")":
                    self.pos += 1
                    break
                raise NewickError(f"Expected ',' or ')' at character {self.pos + 1}")
        node.label = self._label()
        self._space()
        if self.pos < len(self.text) and self.text[self.pos] == ":":
            self.pos += 1
            start = self.pos
            while self.pos < len(self.text) and self.text[self.pos] not in ",();":
                self.pos += 1
            value = self.text[start:self.pos].strip()
            try:
                node.length = float(value)
            except ValueError as exc:
                raise NewickError(f"Invalid branch length: {value!r}") from exc
        return node

    def _label(self) -> str:
        self._space()
        if self.pos >= len(self.text) or self.text[self.pos] in ":,();":
            return ""
        if self.text[self.pos] in "'\"":
            quote = self.text[self.pos]
            self.pos += 1
            result: list[str] = []
            while self.pos < len(self.text):
                char = self.text[self.pos]
                self.pos += 1
                if char == quote:
                    if self.pos < len(self.text) and self.text[self.pos] == quote:
                        result.append(quote)
                        self.pos += 1
                        continue
                    return "".join(result)
                result.append(char)
            raise NewickError("Unclosed quoted label")
        start = self.pos
        while self.pos < len(self.text) and self.text[self.pos] not in ":,();":
            self.pos += 1
        return self.text[start:self.pos].strip()


def parse_newick(path: str | Path) -> Node:
    root = _Parser(Path(path).read_text(encoding="utf-8")).parse()
    assign_ids(root)
    return root


def preorder(root: Node) -> list[Node]:
    result: list[Node] = []

    def visit(node: Node) -> None:
        result.append(node)
        for child in node.children:
            visit(child)

    visit(root)
    return result


def postorder(root: Node) -> list[Node]:
    result: list[Node] = []

    def visit(node: Node) -> None:
        for child in node.children:
            visit(child)
        result.append(node)

    visit(root)
    return result


def assign_ids(root: Node) -> None:
    nodes = preorder(root)
    internal_number = 1
    branch_number = 1
    used: set[str] = set()
    for node in nodes:
        if node.is_leaf:
            node.node_id = node.label
        else:
            candidate = node.label or f"NODE_{internal_number:04d}"
            internal_number += 1
            if candidate in used:
                candidate = f"NODE_{internal_number:04d}"
                internal_number += 1
            node.node_id = candidate
        used.add(node.node_id)
        if node.parent is not None:
            node.branch_id = f"BR{branch_number:06d}"
            branch_number += 1
            node.depth = node.parent.depth + 1


def leaf_labels(root: Node) -> list[str]:
    return [node.label for node in preorder(root) if node.is_leaf]


def as_project_nodes(root: Node) -> list[dict[str, object]]:
    return [
        {
            "id": node.node_id,
            "label": node.label or node.node_id,
            "parent_id": node.parent.node_id if node.parent else None,
            "branch_id": node.branch_id,
            "length": node.length,
            "is_leaf": node.is_leaf,
            "depth": node.depth,
            "children": [child.node_id for child in node.children],
        }
        for node in preorder(root)
    ]
