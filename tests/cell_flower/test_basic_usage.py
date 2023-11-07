import cell_flower as cf
import pandas as pd
import numpy as np
import networkx as nx

def test_basic_run():
    # Initialize underlying graph
    CC = cf.CellComplex(5, [(0,1), (1,2), (2,3), (0,3), (3,4), (0,4)])

    # Cell boundary maps expressed as pd.Series
    # Cell (0,1,2,3)
    square_cell = pd.Series({
            (0,1): 1,
            (1,2): 1,
            (2,3): 1,
            (0,3): -1
        }, index=CC.get_cells(1)).fillna(0).sort_index()

    # Cell (0,3,4)
    triangle_cell = pd.Series({
            (0,3): 1,
            (3,4): 1,
            (0,4): -1
        }, index=CC.get_cells(1)).fillna(0).sort_index()

    # Flows

    flow1 = square_cell * 1 - triangle_cell * .5
    flow2 = square_cell * .3 + triangle_cell * .6

    flows = np.array([flow1.to_numpy(), flow2.to_numpy()])

    CC_prime = cf.cell_inference_approximation(CC, flows, 2, 2, n_clusters=5)
    # Check to see the cells are correct
    assert (0,1,2,3) in CC_prime.get_cells(2)
    assert (0,3,4) in CC_prime.get_cells(2)

def test_nx_compatibility():
    G = nx.karate_club_graph()
    CC, _, _ = cf.nx_graph_to_cc(G)
    G2 = cf.cc_to_nx_graph(CC)
    assert nx.is_isomorphic(G, G2)