
import pytest
from ..myPhylTree import *


def test_format_uniq_suffix_neg_raises_error():
    with pytest.raises(ValueError):
        format_uniq_suffix(-1)

def test_format_uniq_suffix_0_is_A():
    assert format_uniq_suffix(0) == 'A'

def test_format_uniq_suffix_25_is_Z():
    assert format_uniq_suffix(25) == 'Z'

def test_format_uniq_suffix_26_is_BA():
    assert format_uniq_suffix(26) == 'BA'

def test_format_uniq_suffix_27_is_BB():
    assert format_uniq_suffix(27) == 'BB'


class Test_PhylogeneticTree:
    items = {'HomoPan': [('Pan', 3), ('Homo sapiens', 6)],
             'Pan': [('Pan troglodytes', 3), ('Pan paniscus', 3)]}
    root = 'HomoPan'
    officialName = {}
    
    def test_load_from_tuple_skipInit(self):
        # SkipInit has no effect
        PhylogeneticTree((self.items, self.root, self.officialName), skipInit=True)

    def test_load_from_tuple(self):
        PhylogeneticTree((self.items, self.root, self.officialName), skipInit=False)

    def test_load_from_tuple_init(self):
        phyltree = PhylogeneticTree((self.items, self.root, self.officialName), skipInit=True)
        phyltree.reinitTree()
        assert hasattr(phyltree, 'parent')
        assert hasattr(phyltree, 'species')
        assert hasattr(phyltree, 'fileName')
        assert hasattr(phyltree, 'dicLinks')
        assert hasattr(phyltree, 'dicParents')
        assert hasattr(phyltree, 'allDescendants')

