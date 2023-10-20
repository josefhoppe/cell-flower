"""
Functions to find and handle spanning trees.
"""

from collections import defaultdict
from typing import NamedTuple, Tuple, Dict, Set, List
import numpy as np
from queue import Queue

from .cell_complex import CellComplex

class FlowPotentialSpanningTree(NamedTuple):
    """
    A spanning tree with (multi-dimensional) potentials on nodes
    """
    st_id: int
    parent_node: np.ndarray
    node_level: np.ndarray
    node_potential: np.ndarray

    def least_common_ancestor(self, node_pairs: List[Tuple[int,int]]) -> Dict[Tuple[int,int], int]:
        """
        Implementation of Tarjan's off-line lowest common ancestors algorithm
        (https://en.wikipedia.org/wiki/Tarjan%27s_off-line_lowest_common_ancestors_algorithm, https://doi.org/10.1145%2F322154.322161)

        Returns edge (as tuple (a,b)) -> ancestor

        The time complexity is in O(n + m * α(n + m)) for n = |nodes|, m = |node_pairs|
        (α is the inverse of the Ackermann function f(x) = A(x,x))

        Note that Tarjan later gave an improved algorithm with complexity O(n+m).
        Since the Ackermann function grows so fast, α is almost constant.
        It is especially faster than the sorting we perform afterward.
        """
        node_color = np.zeros_like(self.parent_node, dtype=np.int8)
        ancestor = np.arange(self.parent_node.shape[0])
        partition = UnionFind(node_color.shape[0])

        # need to traverse the spanning tree top-down
        children = defaultdict(list)
        root = -1
        for node, parent in enumerate(self.parent_node):
            if node != parent:
                children[parent].append(node)
            else:
                root = node

        # need to efficiently retrieve queries
        # node -> (other, edge[0], edge[1])
        queries = defaultdict(list)
        for (a,b) in node_pairs:
            queries[a].append((b,a,b))
            queries[b].append((a,a,b))

        result = {}
        
        def lca(node):
            for child in children[node]:
                lca(child)
                partition.unite(node, child)
                ancestor[partition[node]] = node
            node_color[node] = 1
            for other, a,b in queries[node]:
                if node_color[other] == 1:
                    result[(a,b)] = ancestor[partition[other]]
        
        lca(root)

        return result

class UnionFind:
    """
    Starts with `items` different sets and supports two operations:
    1. get unique identifier for set (via __getitem__). The identifier is also part of the set.
    2. unite two sets

    Both operations work in O(α(n)) (with amortisation - individual operations may take longer)
    (α is the inverse of the Ackermann function, so practically constant.)
    """
    parent: np.ndarray

    def __init__(self, items: int) -> None:
        self.parent = np.arange(items)

    def __getitem__(self, idx: int) -> int:
        parent = self.parent[idx]
        # check if parent assigned correctly
        if self.parent[parent] == parent:
            return parent
        # find correct root anestor
        while self.parent[parent] != parent:
            parent = self.parent[parent]
        # update all ancestors to point to root
        p = idx
        while p != parent:
            next_p = self.parent[p]
            self.parent[p] = parent
            p = next_p
        return parent

    def unite(self, a: int, b: int):
        """
        unites the set containing a with the set containing b
        """
        self.parent[self[a]] = self[b]


def max_spanning_tree(cell_compl: CellComplex, flows: np.ndarray,
                      weight: np.ndarray | None = None, st_id = 0) ->\
        Tuple[FlowPotentialSpanningTree, List[Tuple[int,int,int]]]:
    """
    Finds the maximal spanning tree according to `sum(abs(flows))`

    weight: edges are added from highest do lowest weight. default `np.sum(np.abs(flows)`

    Returns: (SpanningTree, List[(edge_idx, node, node)])
    """
    tree_incidences: Dict[int, Set[Tuple[int, int]]] = defaultdict(set)
    non_tree_edges = []

    if weight is None:
        weight = np.sum(np.abs(flows), axis=0)

    edge_list = cell_compl.get_cells(1)
    max_weight_edges = np.argsort(weight)[::-1]
    conn_comps = UnionFind(len(cell_compl.get_cells(0)))

    for edge_idx in max_weight_edges:
        a, b = edge_list[edge_idx]
        if conn_comps[a] != conn_comps[b]:
            # add current edge
            tree_incidences[a].add((b, edge_idx))
            tree_incidences[b].add((a, edge_idx + len(edge_list)))
            conn_comps.unite(a, b)
        else:
            non_tree_edges.append((edge_idx, a, b))

    tree, additional_edges = bfs_spanning_tree(
        rnd=None, node_incidences=tree_incidences, edge_count=len(edge_list), root_node=0, flows=flows, st_id = st_id)
    
    if len(additional_edges) > 0:
        # should never happen, just to catch errors
        raise RuntimeError('Max ST is not a tree')


    return tree, non_tree_edges


def dfs_spanning_tree(rnd: np.random.Generator | None,
                      node_incidences: Dict[int, Set[Tuple[int, int]]], edge_count: int, root_node: int,
                      flows: np.ndarray) -> Tuple[FlowPotentialSpanningTree, List[Tuple[int,int,int]]]:
    """
    Finds a spanning tree via depth-first-search.
    if `rnd` is supplied, it is used to randomize the order in which children are traversed.

    `node_incidences` has the same format as AbstractComplex.node_incidences:
    node -> (neighbor, edge_idx)
    where len(edges) is added to edge_idx if the traversal is from the higher node index to the lower.

    returns (spanning_tree, List[(edge_idx, node, node)])
    """
    permutation = rnd.permutation if rnd is not None else lambda x: x

    nodes = len(node_incidences.keys())
    # spanning tree properties
    parent_node = np.full((nodes,), -1, dtype=np.int32)
    parent_node[root_node] = root_node
    node_level = np.full((nodes,), -1, dtype=np.int32)
    node_level[root_node] = 0
    node_potential = np.zeros((nodes, flows.shape[0]), dtype=np.float64)

    non_tree_edges = []

    def dfs(node):
        for (other, edge) in permutation(np.fromiter(node_incidences[node], dtype=tuple)):
            if other != parent_node[node]:
                edge_flows = flows[:, edge % edge_count] * \
                    (-1) ** (edge // edge_count)
                if parent_node[other] == -1:
                    # undiscovered
                    parent_node[other] = node
                    node_level[other] = node_level[node] + 1
                    node_potential[other] = node_potential[node] - edge_flows
                    dfs(other)
                else:
                    non_tree_edges.append((edge % edge_count, node, other))
    dfs(root_node)

    return FlowPotentialSpanningTree(parent_node, node_level, node_potential), non_tree_edges


def bfs_spanning_tree(rnd: np.random.Generator | None,
                      node_incidences: Dict[int, Set[Tuple[int, int]]], edge_count: int, root_node: int,
                      flows: np.ndarray, st_id = 0) -> Tuple[FlowPotentialSpanningTree, List[Tuple[int,int,int]]]:
    """
    Finds a spanning tree via depth-first-search.
    if `rnd` is supplied, it is used to randomize the order in which children are traversed.

    `node_incidences` has the same format as AbstractComplex.node_incidences:
    node -> (neighbor, edge_idx)
    where len(edges) is added to edge_idx if the traversal is from the higher node index to the lower.

    returns (spanning_tree, List[(edge_idx, node, node)])
    """
    permutation = rnd.permutation if rnd is not None else lambda x: x

    nodes = max(node_incidences.keys()) + 1
    # spanning tree properties
    parent_node = np.full((nodes,), -1, dtype=np.int32)
    parent_node[root_node] = root_node
    node_level = np.full((nodes,), -1, dtype=np.int32)
    node_level[root_node] = 0
    node_potential = np.zeros((nodes, flows.shape[0]), dtype=np.float64)

    non_tree_edges = []

    def bfs(root):
        queue = Queue()
        queue.put(root)
        while not queue.empty():
            node = queue.get()
            # randomize because the BFS tree with root n is not unique:
            #  root
            #  /  \
            # a -- b
            #  \  / <-- two possible BFS trees: in one, c is the child of a, in the other it is the child of b
            #   c
            for (other, edge) in permutation(np.fromiter(node_incidences[node], dtype=tuple)):
                if other != parent_node[node]:
                    edge_flows = flows[:, edge %
                                       edge_count] * (-1) ** (edge // edge_count)
                    if parent_node[other] == -1:
                        # undiscovered
                        parent_node[other] = node
                        node_level[other] = node_level[node] + 1
                        node_potential[other] = node_potential[node] - \
                            edge_flows
                        queue.put(other)
                    else:
                        non_tree_edges.append((edge % edge_count, node, other))
    bfs(root_node)

    return FlowPotentialSpanningTree(st_id, parent_node, node_level, node_potential), non_tree_edges