
from ..myProteinTree import *


def test_init_empty_ProteinTree():
    prottree = ProteinTree()


def test_init_ProteinTree():
    prottree = ProteinTree(data={1: [(2, 0.01), (3, 0.01)]},
                           info={1: {'taxon': 'Pan'},
                                 2: {'taxon': 'Pan troglodytes'},
                                 3: {'taxon': 'Pan paniscus'}},
                           root=1)
    assert prottree.root is not None
    assert prottree.data
    assert prottree.info
