# Cell FLOWer

<img align="right" width="200" style="margin-top:-5px" src="https://raw.githubusercontent.com/josefhoppe/cell-flower/main/readme_src/LOGO_ERC-FLAG_FP.png">

[![Published at Learning on Graphs Conference 2023](https://img.shields.io/badge/Published_at-LoG_2023-003380.svg)](https://openreview.net/forum?id=qix189lq5D)
[![Best Paper at LoG 2023 - LoG Conference on Twitter](https://img.shields.io/badge/LoG_2023:_Best_Paper-003380.svg?logo=data:image/png;base64,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)](https://twitter.com/LogConference/status/1729167588851143023)
[![arXiv:2309.01632](https://img.shields.io/badge/arXiv-2309.01632-b31b1b.svg?logo=arxiv)](https://arxiv.org/abs/2309.01632)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/josefhoppe/cell-flower/blob/main/LICENSE)
![Python version](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Fjosefhoppe%2Fcell-flower%2Fmain%2Fpyproject.toml&logo=python&logoColor=ffd242)
![Build status](https://github.com/josefhoppe/cell-flower/actions/workflows/python-test.yml/badge.svg)
[![Package version on PyPI](https://img.shields.io/pypi/v/cell-flower?logo=pypi&logoColor=ffd242)](https://pypi.org/project/cell-flower/)

Cell FLOWer processes edge flows using cell complexes.
It was developed for the paper [*Representing Edge Flows on Graphs via Sparse Cell Complexes*](https://arxiv.org/abs/2309.01632).
You can find the evaluation workflow [here](https://github.com/josefhoppe/edge-flow-cell-complexes).
Install it using:

```
pip install cell-flower
```

How to use it (Also check out the complete [examples](https://github.com/josefhoppe/cell-flower/tree/main/examples)):

```python
import cell_flower as cf
import networkx as nx

G = ... # nx.Graph
CC, node_list, node_index_map = cf.nx_graph_to_cc(G) # converts nodes in G to int and also returns the mapping
flows = ... # np.ndarray, shape (samples, edges)

CC_prime = cf.cell_inference_approximation(CC, flows, 2, 2, n_clusters=5)
# Check to see the cells that were recovered; map back to node labels
[ cf.index_to_cell(cell, node_list) for cell in CC_prime.get_cells(2) ]
```

If you use Cell FLOWer, please cite the following paper:

```
@inproceedings{
      hoppe2023representing,
      title={Representing Edge Flows on Graphs via Sparse Cell Complexes},
      author={Josef Hoppe and Michael T Schaub},
      booktitle={The Second Learning on Graphs Conference},
      year={2023},
      url={https://openreview.net/forum?id=qix189lq5D}
}
```

## Acknowledgements

Funded by the European Union (ERC, HIGH-HOPeS, 101039827). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them.
