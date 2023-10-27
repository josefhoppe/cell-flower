# Cell FLOWer

<img align="right" width="200" style="margin-top:-5px" src="https://raw.githubusercontent.com/josefhoppe/cell-flower/main/readme_src/LOGO_ERC-FLAG_FP.png">

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

CC = cf.CellComplex([(0,1), (1,2), (2,3), (0,3), (3,4), (0,4)])
flows = ...

CC_prime = cf.cell_inference_approximation(CC, flows, 2, 2, n_clusters=5)
# Check to see the cells that were recovered
CC_prime.get_cells(2)
```

If you use Cell FLOWer, please cite the following paper:

```
@misc{hoppe2023representing,
      title={Representing Edge Flows on Graphs via Sparse Cell Complexes}, 
      author={Josef Hoppe and Michael T. Schaub},
      year={2023},
      eprint={2309.01632},
      archivePrefix={arXiv},
      primaryClass={cs.SI}
}
```

## Acknowledgements

Funded by the European Union (ERC, HIGH-HOPeS, 101039827). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them.
