"""
Cell FLOWer. Implements algorithms and requirements from [1].

[1]: Josef Hoppe and Michael T. Schaub: Representing Edge Flows on Graphs via Sparse Cell Complexes. In: The Second Learning on Graphs Conference, 2023. https://openreview.net/pdf?id=qix189lq5D
"""

from .cell_complex import CellComplex, normalize_cell, nx_graph_to_cc, cc_to_nx_graph, map_cells_to_index, cell_to_index, index_to_cell
from .detection import CellCandidateHeuristic, CellSearchFlowNormalization, cell_inference_approximation