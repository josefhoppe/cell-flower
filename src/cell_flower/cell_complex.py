from itertools import combinations
from typing import DefaultDict, Any, Literal, Callable

import numpy as np
from scipy.sparse import csc_array, lil_array, hstack

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

def calc_edges(cell: tuple) -> list[tuple]:
    """
    Calculates all edges as 2-tuples of nodes from a 2-dim cell
    """
    res = []
    for i, node in enumerate(cell):
        j = (i + 1) % len(cell)
        res.append(normalize_cell((node, cell[j])))
    return res

def cell_to_index(cell: tuple, node_index_map: dict[Any, int]) -> tuple:
    """Maps the given cell tuple to an int tuple according to the node index map.

    Usage: see `map_cells_to_index`.
    
    Returns: tuple of ints representing the cell.
    """
    return normalize_cell(tuple([node_index_map[n] for n in cell]))

def index_to_cell(cell: tuple, node_list: list[Any]) -> tuple:
    """Maps the given index tuple to a cell tuple according to the node list.

    Usage: see `map_cells_to_index`.
    
    Returns: tuple of original objects representing the cell.
    """
    return tuple([node_list[i] for i in cell])

def map_cells_to_index(cells: list[tuple], sort_key: Literal[None] | Callable[[Any], Any] = None) -> tuple[int, list[tuple], list[Any], dict[Any,int]]:
    """Maps the given cell tuples to int tuples and a dictionary for reverse mapping

    Usage:

    ```python
    import cell_flower as cf
    cells = [('A', 'B'), ('B','C','D')]
    node_count, cells, node_list, node_index_map = cf.map_cells_to_int(cells)
    CC = cf.CellComplex(node_count, cells)
    print(cf.index_to_cell((0,1,2), node_list)) # ('A', 'B', 'C') [ordered like `normalize_cell` for usage in CC]
    print(cf.cell_to_index(('A', 'B', 'C'), node_list)) # (0, 1, 2) [ordered as given for external usage]
    ```

    Parameters:
    - `cells`: List of all cells (as tuples of any type). May include 0-cells as tuples of length 1.
    - `sort_key`: Sort key used to sort nodes by passing it as the `key` parameter to `sorted()`. If `None`, comparable types will still be sorted.
    
    Returns: number of nodes, list of 1-cells and 2-cells as tuples, ordered list of nodes, dictionary to map nodes to int representation
    """
    node_list = list(set().union(*cells))
    if sort_key != None:
        node_list = sorted(node_list, key=sort_key)
    else:
        try:
            node_list = sorted(node_list)
        except:
            pass
    node_index_map = {
        node: idx for node, idx in enumerate(node_list)
    }
    index_cells = [ cell_to_index(cell, node_index_map) for cell in cells if len(cell) > 1 ]
    return len(node_list), index_cells, node_list, node_index_map

class CellComplex():
    """
    Implementation of a cell complex

    Only supports cells of dimension 0, 1, or 2 (i.e., nodes, edges, polygons)
    """
    cell_order_map: dict[int, list[tuple[int]]]
    __cell_index: dict[tuple[int], int]
    embedding = None
    __cell_boundary_map: csc_array
    __edge_boundary_map: csc_array
    __node_incidences: dict[int, dict[int, int]]
    __triangles: np.ndarray[tuple[tuple,csc_array]] | None = None

    @property
    def triangles(self) -> np.ndarray[tuple[tuple,csc_array]]:
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
            return CellComplex(-1, [], cell_order_map=cell_order_map, 
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
        for i in range(1,3):
            cells += self.cell_order_map[i]
        cells += [cell]
        return CellComplex(len(self.cell_order_map[0]), cells)
    
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
        new_compl = CellComplex(-1, [], cell_order_map=cell_order_map, 
                cell_boundary_map=hstack((self.__cell_boundary_map, cell_boundary)),
                edge_boundary_map=self.__edge_boundary_map, node_incidences=self.__node_incidences, cell_index=cell_index)
        if self.__triangles is not None:
            new_compl.__triangles = self.__triangles
        return new_compl

    def __init__(self, node_count: int, cells: list[tuple[int]], **kwargs):
        """
        Initializes a new Cell complex, either by calculating everything (if node_count != -1) or by taking all properties from **kwargs
        """
        if node_count == -1:
            # 'manual' initialization (see also: add_cell_fast)
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
        
        def assert_cell(order, cell):
            if not cell in self.__cell_index[order]:
                raise RuntimeError(f"Unknown but required cell {cell}. Please use `cf.map_cells_to_index()`.")

        for node in range(node_count):
            # TODO: Refactor to replace 0-cells with node_count attribute?
            add_cell(0, (node,))

        for cell in cells:
            if len(cell) == 1:
                assert_cell(0, cell)
            elif len(cell) == 2:
                assert_cell(0, (cell[0],))
                assert_cell(0, (cell[1],))
                add_cell(1, cell)
            else:
                norm_cell = normalize_cell(cell)
                add_cell(2, norm_cell)
                for point in norm_cell:
                    assert_cell(0, (point,))
                for edge in calc_edges(norm_cell):
                    add_cell(1, edge)
        
        # calculate cell boundaries
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

    def get_cells(self, k: int = 1) -> list:
        """
        returns the cells $C_k$ of the given dimension `k`.
        """
        return self.cell_order_map[k]

    def boundary_map(self, k: int = 1) -> np.ndarray | csc_array:
        """The boundary map $B_k$. Sparse for $k \in {1,2}$."""
        if k == 1:
            return self.__edge_boundary_map
        if k == 2:
            return self.__cell_boundary_map
        if k == 0:
            return np.zeros(shape=(0,len(self.get_cells(k))), dtype=np.int32)
        raise NotImplementedError('Only cells up to dimension 2 (nodes, edges, polygons) are supported.')
    
    def node_incidences(self) -> dict[int, set[tuple[int, int]]]:
        """
        dictionary `node` -> `edge`; edges are represented as tuples of nodes
        """
        return self.__node_incidences
    
    def __node_incidences(self) -> dict[int, set[tuple[int, int]]]:
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

def nx_graph_to_cc(G) -> tuple[CellComplex, list[Any], dict[Any,int]]:
    """
    Converts the networkx Graph to a Cell Complex.

    Always index-converts the given nodes using `map_cells_to_index`, consult that function for more info.
    Also returns the node list and node index map from `map_cells_to_index`.

    Returns: Tuple (cell complex, ordered list of nodes, dictionary to map nodes to int representation)
    """
    cells = [(node, ) for node in G.nodes] + [e for e in G.edges]
    node_count, idx_cells, node_list, node_index_map = map_cells_to_index(cells)
    return CellComplex(node_count, idx_cells), node_list, node_index_map