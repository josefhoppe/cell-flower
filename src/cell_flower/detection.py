import itertools
from typing import Dict, Generator, List, NamedTuple, Set, Tuple
from collections import defaultdict
from heapq import heapify, heappop, heappush
import numpy as np

from numpy.linalg import norm
from sklearn.metrics import pairwise_distances
from sklearn.cluster import KMeans
from scipy.sparse.linalg import lsmr
from scipy.sparse import lil_array, hstack, csc_array
from enum import Enum
from queue import Queue
from .cell_complex import CellComplex, normalize_cell
from .spanning_trees import *

def project_flow(A: np.ndarray, flow: np.ndarray):
    """
    Projects the flow into Im(A).
    """
    flow_star, _, _, _, _, _, _, _ = lsmr(A, flow)
    return A @ flow_star

def score_cell_flows(cell_compl: CellComplex, flows: List[np.array], cell_boundary: csc_array) -> np.ndarray:
    """
    Scores the individual flows, i.e., calculates the element-wise sum of squares of the approximation error for each flow.
    """
    old_boundary = cell_compl.boundary_map(2)
    new_boundary = hstack((old_boundary, cell_boundary))

    scores = np.zeros((len(flows),), dtype=np.float64)
    for idx, flow in enumerate(flows):
        harmonic = flow - project_flow(new_boundary, flow)
        scores[idx] = np.sum(np.square(harmonic))
    return scores

def score_cells_multiple(cell_compl: CellComplex, flows: List[np.array], cells: List[Tuple[tuple, csc_array]], length_penalty: float = 0.01) -> tuple[np.ndarray, list[tuple]]:
    """
    Scores multiple cells on all given flows.
    """
    no_boundary = csc_array((len(cell_compl.get_cells(1)), 0))
    empty_cell = ((), no_boundary)
    cell_idx = [cell[0] for cell in cells + [empty_cell]]
    results = [score_cell_flows(cell_compl, flows, cell[1]) + length_penalty * len(cell)
               for cell in cells + [empty_cell]]
    return np.array(results), cell_idx

def inv_edge(e: Tuple[int, int]):
    return (e[1], e[0])

def get_induced_cycle(edge: Tuple[int, int], st: FlowPotentialSpanningTree) -> tuple:
    """
    Gets the cycle induced by adding edge to the spanning tree modeled by node_level and parent_node
    """
    left = []
    right = []

    a = edge[0]
    b = edge[1]

    if st.node_level[a] < st.node_level[b]:
        a = edge[1]
        b = edge[0]

    while st.node_level[a] > st.node_level[b]:
        left.append(a)
        a = st.parent_node[a]

    while a != b:
        left.append(a)
        right.append(b)
        a = st.parent_node[a]
        b = st.parent_node[b]

    left.append(a)
    return normalize_cell(tuple(left + right[::-1]))

class CellSearchMethod(Enum):
    """
    Methods to find potential cells in cell complexes.
    DFS: depth-first-search to construct a spanning tree, choosing the cell 
         with most circle flow from all cycles induced by adding a single edge to the spanning tree.
    BFS: breadth-first-search (see DFS for details)
    MAX: maximum spanning tree (by `sum(abs(edge_flows))`)
    CLUSTER: cluster by harmonic flow values, then build spanning tree with edges that have maximum similarity to 
    TRIANGLES: consider all triangles
    GROUND_TRUTH: consider all true cells
    """
    DFS = 1
    BFS = 2
    MAX = 3
    CLUSTER = 4
    TRIANGLES = 5


class CellSearchFlowNormalization(Enum):
    """
    How to normalize the integrated cell flows (score = integrated flow / normalization)

    NONE: Don't normalize (i.e., take complete integrated flow)
    LEN: Length of the cell. Efficiently implemented, does not increase asymptotic complexity.
    """
    NONE = 1
    LEN = 2


def find_max_cycles(cell_compl: CellComplex, spanning_tree: FlowPotentialSpanningTree, flows: np.ndarray,
                    edges: List[Tuple[int, int, int]],
                    flow_norm: CellSearchFlowNormalization = CellSearchFlowNormalization.LEN,
                    cycle_flows_heap: list | None = None) -> List[Tuple[float, FlowPotentialSpanningTree, int, int, int]]:
    """
    Checks all cycles induced by `edges` on `spanning_tree` and calculates the total flow (w.r.t `flow_norm`)

    returns a heap of tuples (-normalized_flow, edge_idx, node_a, node_b)  
    Note that the implementation of `heapq` is a min-heap, therefore the negative of the flow is added.  
    However, typically, the flow is only required to `heappop`.
    """
    if cycle_flows_heap is None:
        cycle_flows_heap = []
    # empty or already a heap, so not necessary: heapify(cycle_flows_heap)

    node_potential = spanning_tree.node_potential
    node_level = spanning_tree.node_level
    edge_count = len(cell_compl.get_cells(1))

    # edge -> least common ancestor
    lca: Dict[Tuple[int, int], int]

    if flow_norm == CellSearchFlowNormalization.LEN:
        # O(nodes + edges) operation, so only perform if necessary
        lca = spanning_tree.least_common_ancestor([edge[1:] for edge in edges])

    def add_cycle(edge: int, node: int, other: int, edge_flows: np.ndarray):
        total_flow = np.sum(
            np.abs(node_potential[node] - edge_flows - node_potential[other]))
        if flow_norm == CellSearchFlowNormalization.LEN:
            total_flow /= node_level[node] + node_level[other] + 1 - 2 * node_level[lca[(node, other)]]
        heappush(cycle_flows_heap, (-total_flow, spanning_tree, edge, node, other))

    for edge, node, other in edges:
        edge_flows = flows[:, edge % edge_count] * (-1) ** (edge // edge_count)
        add_cycle(edge, node, other, edge_flows)

    return cycle_flows_heap

def __finalize_candidates(cell_compl: CellComplex, heap: list[tuple[float, FlowPotentialSpanningTree, int, int, int]], 
                        count: int) -> list[tuple[tuple, csc_array, FlowPotentialSpanningTree]]:
    """
    Takes a heap resulting from a candidate search on one or more `FlowPotentialSpanningTree`s and returns a list of candidates for evaluation.

    Also de-duplicates.
    """
    result = []
    found_cycles = set()

    while len(result) < count:
        _, spanning_tree, _, u, v = heappop(heap)
        cycle = get_induced_cycle((u,v), spanning_tree)
        if cycle not in found_cycles:
            result.append((cycle, cell_compl.calc_cell_boundary(cycle), spanning_tree))
            found_cycles.add(cycle)
    
    return result

def cell_candidate_search_st(rnd: np.random.Generator, cell_compl: CellComplex, result_count: int, flows: np.ndarray, method: CellSearchMethod = CellSearchMethod.DFS,
                            max_len=np.inf, n_clusters = 11,
                             flow_norm: CellSearchFlowNormalization = CellSearchFlowNormalization.NONE, random_seed : int | None = None) -> List[Tuple[tuple, csc_array, FlowPotentialSpanningTree]]:
    """
    Searches for Cell candidates using spanning trees.
    The concrete implementation depends on `method`:
    BFS: randomly select `result_count` breadth-first-search spanning trees and get the cell with most integrated flow from each
    DFS: randomly select `result_count` depth-first-search spanning trees and get the cell with most integrated flow from each
    MAX: construct the maximum spanning tree (according to `sum(abs(flows))`) and get the `result_count` cells with the most integrated flow.
    see `cell_candidate_search_rnd_st` and `cell_candidate_search_max_st` respectively for more info.
    """
    max_cycle_heap = []
    if method == CellSearchMethod.MAX:
        spanning_tree, edges = max_spanning_tree(cell_compl, flows)

        max_cycle_heap = find_max_cycles(cell_compl, spanning_tree, flows, edges, flow_norm)
    elif method == CellSearchMethod.CLUSTER:
        # since the orientation of edges is arbitrary, we need to cluster both with and opposite the orientation
        embed = np.concatenate((flows, -flows), axis=1).T
        kmeans = KMeans(n_clusters=n_clusters, n_init='auto', random_state=random_seed).fit(embed)
        dist_center_edges_pos = pairwise_distances(kmeans.cluster_centers_, flows.T)
        dist_center_edges_neg = pairwise_distances(kmeans.cluster_centers_, -flows.T)
        dist_center_edges = np.minimum(dist_center_edges_pos, dist_center_edges_neg)
        for idx, similarity in enumerate(dist_center_edges):
            spanning_tree, edges = max_spanning_tree(cell_compl, flows, -1 * similarity, st_id = idx)

            find_max_cycles(cell_compl, spanning_tree, flows, edges, flow_norm, max_cycle_heap)
    else:
        incidences = cell_compl.node_incidences()
        nodes = np.fromiter(incidences.keys(), dtype=int)
        edge_count = len(cell_compl.get_cells(1))

        for _ in range(result_count):
            root_node = rnd.choice(nodes)
            if method == CellSearchMethod.BFS:
                spanning_tree, edges = bfs_spanning_tree(rnd, incidences, edge_count, root_node, flows)
            else: # CellSearchMethod.DFS
                spanning_tree, edges = dfs_spanning_tree(rnd, incidences, edge_count, root_node, flows)

            find_max_cycles(cell_compl, spanning_tree, flows, edges, flow_norm, max_cycle_heap)
    return __finalize_candidates(cell_compl, max_cycle_heap, result_count)
