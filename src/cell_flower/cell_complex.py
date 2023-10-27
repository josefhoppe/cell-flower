from abc import ABC, abstractmethod
from itertools import chain, combinations
from typing import DefaultDict, Dict, Iterable, List, Set, Tuple

import numpy as np
import numpy.linalg as linalg
from scipy.sparse import csc_array, lil_array, hstack
from scipy.sparse.linalg import lsmr

from .matrix_utils import neg_part, normalize, pos_part

def normalize_cell(cell: tuple) -> tuple:
    """
    Transforms the cell (represented as a sequence of nodes) to normalized form.

    Normalized form: `cell_prime := (m, a, [...], b)` s.t.:
    
    - `m := min(cell)`
    - `a < b`
    - `cell_prime` represents the same cycle as `cell`

    returns `cell_prime`
    """
    min_index = cell.index(min(cell))
    shifted = cell[min_index:] + cell[:min_index]
    if shifted[-1] < shifted[1]:
        shifted = shifted[::-1]
        shifted = shifted[-1:] + shifted[:-1]
    return shifted

def calc_edges(cell: tuple) -> List[tuple]:
    """
    Calculates all edges as 2-tuples of nodes from a 2-dim cell
    """
    res = []
    for i, node in enumerate(cell):
        j = (i + 1) % len(cell)
        res.append(normalize_cell((node, cell[j])))
    return res

class CellComplex():
    """
    Implementation of a cell complex

    Only supports cells of dimension 0, 1, or 2 (i.e., nodes, edges, polygons)
    """
    cell_order_map: Dict[int, List[Tuple[int]]]
    __cell_index: Dict[Tuple[int], int]
    embedding = None
    __cell_boundary_map: csc_array
    __edge_boundary_map: csc_array
    __node_incidences: Dict[int, Dict[int, int]]
    __triangles: np.ndarray[Tuple[tuple,csc_array]] | None = None

    @property
    def triangles(self) -> np.ndarray[Tuple[tuple,csc_array]]:
        """
        All triangles in the complex.
        """
        if self.__triangles is None:
            self.__calc_triangles()
        return self.__triangles
    
    def __calc_triangles(self):
        candidate_set = set()

        edge_set = set(self.get_cells(1))
        incidences = self.node_incidences()

        for node, neighbors in incidences.items():
            for (a, _), (b, _) in combinations(neighbors, 2):
                if normalize_cell((a,b)) in edge_set:
                    candidate_set.add(normalize_cell((node,a,b)))

        self.__triangles = np.ndarray(len(candidate_set), dtype=tuple)
        
        for idx, cycle in enumerate(candidate_set):
            self.__triangles[idx] = (cycle, self.calc_cell_boundary(cycle))

    def skeleton(self, order: int = 1) -> "CellComplex":
        """
        Returns the skeleton of the complex.

        Note: Currently, only the 1-skeleton is implemented efficiently.
        """
        if order == 1:
            cell_order_map = DefaultDict(list)
            cell_index = DefaultDict(lambda: {})
            for i in range(order + 1):
                cell_order_map[i] += self.cell_order_map[i]
                cell_index[i] = self.__cell_index[i].copy()

            empty_boundaries = csc_array((len(self.get_cells(1)), 0), dtype=np.int32)
            return CellComplex(None, cell_order_map=cell_order_map, 
                    cell_boundary_map=empty_boundaries, node_incidences=self.__node_incidences,
                    edge_boundary_map=self.__edge_boundary_map, cell_index=cell_index)
        # very slow
        # todo implement fast version for other orders
        cells = []
        for i in range(order + 1):
            cells += self.cell_order_map[i]
        return CellComplex(cells)

    def add_cell(self, cell: tuple) -> "CellComplex":
        """
        Returns a cell complex with `cell` added.

        Note: the implementation is currently not very efficient, use `add_cell_fast` to add 2-cells efficiently.
        """
        cells = []
        for i in range(3):
            cells += self.cell_order_map[i]
        cells += [cell]
        return CellComplex(cells)
    
    def add_cell_fast(self, cell:tuple, cell_boundary: csc_array) -> "CellComplex":
        """
        Adds a cell efficiently if the boundary is given.

        Returns a new cell complex with the added cell.
        """
        if len(cell) < 3:
            raise NotImplementedError('Can only add 2-cells')
        cell_order_map = DefaultDict(list)
        cell_index = DefaultDict(lambda: {})
        for i in range(3):
            cell_order_map[i] += self.cell_order_map[i]
            cell_index[i] = self.__cell_index[i].copy()
        cell_order_map[2].append(cell)
        cell_index[2][cell] = len(cell_order_map[2]) - 1
        new_compl = CellComplex(None, cell_order_map=cell_order_map, 
                cell_boundary_map=hstack((self.__cell_boundary_map, cell_boundary)),
                edge_boundary_map=self.__edge_boundary_map, node_incidences=self.__node_incidences, cell_index=cell_index)
        if self.__triangles is not None:
            new_compl.__triangles = self.__triangles
        return new_compl

    def __init__(self, cells: List[Tuple[int]] | None, **kwargs):
        """
        Initializes a new Cell complex, either by calculating everything (if cells != None) or by taking all properties from **kwargs
        """
        if cells is None:
            # 'manual' initialization
            self.cell_order_map = kwargs['cell_order_map']
            self.__cell_boundary_map = kwargs['cell_boundary_map']
            self.__edge_boundary_map = kwargs['edge_boundary_map']
            self.__node_incidences = kwargs['node_incidences']
            self.__cell_index = kwargs['cell_index']
            return
        self.cell_order_map = DefaultDict(list)
        self.__cell_index = DefaultDict(lambda: {})

        def add_cell(order, cell):
            if not cell in self.__cell_index[order]:
                self.cell_order_map[order].append(cell)
                self.__cell_index[order][cell] = len(self.cell_order_map[order]) - 1

        for cell in cells:
            if len(cell) == 1:
                add_cell(0, cell)
            elif len(cell) == 2:
                add_cell(0, (cell[0],))
                add_cell(0, (cell[1],))
                add_cell(1, cell)
            else:
                norm_cell = normalize_cell(cell)
                add_cell(2, norm_cell)
                for point in norm_cell:
                    add_cell(0, (point,))
                for edge in calc_edges(norm_cell):
                    add_cell(1, edge)
        
        # calculate cell boundaries
        node_count = len(self.get_cells(0))
        edge_count = len(self.get_cells(1))
        cell_count = len(self.get_cells(2))
        cell_boundary = lil_array((edge_count, cell_count), dtype=np.int32)
        for upper_idx, cell in enumerate(self.get_cells(2)):
            for i, _ in enumerate(cell):
                j = (i + 1) % len(cell)
                lower_idx = self.__cell_index[1][tuple(sorted([cell[i], cell[j]]))]
                orientation = 1 if cell[i] < cell[j] else -1
                cell_boundary[lower_idx, upper_idx] = orientation
        self.__cell_boundary_map = csc_array(cell_boundary)

        edge_boundary = lil_array((node_count, edge_count), dtype=np.int32)
        for upper_idx, edge in enumerate(self.get_cells(1)):
            edge_boundary[self.__cell_index[0][(edge[0],)], upper_idx] = -1
            edge_boundary[self.__cell_index[0][(edge[1],)], upper_idx] = 1
        self.__edge_boundary_map = csc_array(edge_boundary)
        self.__node_incidences = self.__node_incidences()

    def calc_cell_boundary(self, cell: tuple) -> csc_array:
        """
        Get the boundary map for a hypothetical cell. Re-calculates the boundary map, so use the existing boundary map for existing cells.
        """
        edge_count = len(self.get_cells(1))
        cell_boundary = lil_array((edge_count, 1), dtype=np.int32)
        for i, _ in enumerate(cell):
            j = (i + 1) % len(cell)
            lower_idx = self.__cell_index[1][tuple(sorted([cell[i], cell[j]]))]
            orientation = 1 if cell[i] < cell[j] else -1
            cell_boundary[lower_idx, 0] = orientation
        return csc_array(cell_boundary)

    def get_cells(self, k: int = 1) -> List:
        """
        returns the cells $C_k$ of the given dimension `k`.
        """
        return self.cell_order_map[k]

    def boundary_map(self, k: int = 1) -> np.ndarray | csc_array:
        """
        The boundary map $B_k$. Sparse for $k \in {1,2}$.
        """
        if k == 1:
            return self.__edge_boundary_map
        if k == 2:
            return self.__cell_boundary_map
        if k == 0:
            return np.zeros(shape=(0,len(self.get_cells(k))), dtype=np.int32)
        raise NotImplementedError('Only cells up to dimension 2 (nodes, edges, polygons) are supported.')
    
    def node_incidences(self) -> Dict[int, Set[Tuple[int, int]]]:
        """
        Dictionary `node` -> `edge`; edges are represented as tuples of nodes
        """
        return self.__node_incidences
    
    def __node_incidences(self) -> Dict[int, Set[Tuple[int, int]]]:
        """
        For all nodes, get all adjacent nodes and the index of the connecting edge.
        ```
                            ↓ index of the edge. len(edges) is added to the index iff opposite to edge orientation
        per node: (other, edge_index)
                    ↑ other node (as number / name)
        ```
        """
        res = DefaultDict(set)
        edge_count = len(self.get_cells(1))
        for idx, (a, b) in enumerate(self.get_cells(1)):
            res[a].add((b, idx))
            res[b].add((a, idx + edge_count))
        return res
    
def cc_to_nx_graph(CC: CellComplex):
    """
    Converts the 1-skeleton of the given CC to the networkx.Graph format.

    Requires networkx to be on the path.
    """
    import networkx as nx
    G = nx.Graph()
    # add separately to ensure correct order
    G.add_nodes_from([c[0] for c in CC.get_cells(0)])
    G.add_edges_from(CC.get_cells(1))
    return G

def nx_graph_to_cc(G) -> CellComplex:
    """
    Converts the networkx Graph to a Cell Complex
    """
    cells = [(node, ) for node in G.nodes] + [e for e in G.edges]
    return CellComplex(cells)