
import pytest

import topiary
import topiary.draw.prettytree as prettytree
import ete3

import toytree
import numpy as np


def test_PrettyTree():

    # Read as string
    tree = "(((A:1.0,B:4.0)AB:1.0,((C:1.0,D:1.0)CD:1.0,E:1.0)CDE:1.0)ABCDE:1.0,(F:1.0,G:1.0)FG)ABCDEFG;"
    pt = prettytree.PrettyTree(T=tree)
    assert pt.tT.ntips == 7

    # read as ete3 tree
    tree = "(((A:1.0,B:4.0)AB:1.0,((C:1.0,D:1.0)CD:1.0,E:1.0)CDE:1.0)ABCDE:1.0,(F:1.0,G:1.0)FG)ABCDEFG;"
    T = ete3.Tree(tree,format=1)
    pt = prettytree.PrettyTree(T=tree)
    assert pt.tT.ntips == 7

    pt = prettytree.PrettyTree(T=tree,font_size=20)
    assert pt._font_size == 20

    pt = prettytree.PrettyTree(T=tree,stroke_width=20)
    assert pt._stroke_width == 20

    pt = prettytree.PrettyTree(T=tree,vertical_pixels_per_tip=100)
    assert pt._vertical_pixels_per_tip == 100

    # check artwork parameters
    bad_float = [-1,None,int,float,[],(1,)]
    for b in bad_float:
        with pytest.raises(ValueError):
            pt = prettytree.PrettyTree(T=tree,font_size=b)
        with pytest.raises(ValueError):
            pt = prettytree.PrettyTree(T=tree,stroke_width=b)
        with pytest.raises(ValueError):
            pt = prettytree.PrettyTree(T=tree,vertical_pixels_per_tip=b)

def test_integrated_single():

    T = toytree.rtree.rtree(50)

    tree_data = list(T.idx_dict.keys())
    v1 = dict([(k,i) for i, k in enumerate(tree_data)])
    v2 = dict([(k,0) for i, k in enumerate(tree_data)])

    T = T.set_node_values(feature="test_feature",values=v1)
    T = T.set_node_values(feature="other_feature",values=v2)

    pt = topiary.draw.PrettyTree(T,tip_labels_align=True)

    pt.draw_nodes("test_feature")
    pt.draw_nodes("other_feature",color="pink",size=5)
    pt.draw_node_labels("test_feature")
    pt.draw_scale_bar()
    pt.draw_node_legend()

def test_integrated_gradient():

    T = toytree.rtree.rtree(50)

    # Change in toytree api moved idx_dict to _idx_dict
    if not hasattr(T,"idx_dict"):
        T.idx_dict = T._idx_dict

    for n in T.idx_dict:
        T.idx_dict[n].add_feature("test_feature",n)
        T.idx_dict[n].add_feature("other_feature",len(T.idx_dict)-n)

    pt = topiary.draw.PrettyTree(T,tip_labels_align=True)

    pt.draw_nodes("test_feature",color=("white","red"),size=(5,20))
    pt.draw_nodes("other_feature",color=("white","blue"),size=(7,7))
    pt.draw_node_labels("test_feature")
    pt.draw_scale_bar()
    pt.draw_node_legend()


def test_integrated_categories():

    T = toytree.rtree.rtree(50)

    # Change in toytree api moved idx_dict to _idx_dict
    if not hasattr(T,"idx_dict"):
        T.idx_dict = T._idx_dict

    for n in T.idx_dict:
        T.idx_dict[n].add_feature("test_feature",np.random.choice(["A","B","C","D"]))
        T.idx_dict[n].add_feature("other_feature",len(T.idx_dict)-n)

    pt = topiary.draw.PrettyTree(T,tip_labels_align=True)

    pt.draw_nodes("test_feature",
                  color={"A":"#00FF00","B":"pink","C":np.array((1,1,0)),"D":np.array((0.5,0.5,0.5))},
                  size={"A":5,"B":10,"C":20,"D":30})
    pt.draw_nodes("other_feature",color=("white","blue"),size=(7,7))
    pt.draw_node_labels("test_feature")
    pt.draw_scale_bar()
    pt.draw_node_legend()
    
def test_PrettyTree__get_node_values():

    pass

def test_PrettyTree_canvas():

    pass

def test_PrettyTree_default_size():

    pass

def test_PrettyTree_draw_node_labels():

    pass

def test_PrettyTree_draw_node_legend():

    pass

def test_PrettyTree_draw_nodes():

    pass

def test_PrettyTree_draw_scale_bar():

    pass

def test_PrettyTree_legend_ax():

    pass

def test_PrettyTree_plotted_properties():

    pass

def test_PrettyTree_render():

    pass

def test_PrettyTree_as_html():
    
    pass

def test_PrettyTree_tT():

    pass

def test_PrettyTree_tree_ax():

    pass

def test_PrettyTree_tree_mark():

    pass
