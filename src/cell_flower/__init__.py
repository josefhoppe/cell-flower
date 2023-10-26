"""
Cell FLOWer. Implements algorithms and requirements from [1].

[1]: Josef Hoppe and Michael T. Schaub: Representing Edge Flows on Graphs via Sparse Cell Complexes. In: arXiv pre-prints, 2023. https://arxiv.org/abs/2309.01632
"""

from .cell_complex import CellComplex, normalize_cell
from .detection import CellCandidateHeuristic, CellSearchFlowNormalization, cell_inference_approximation