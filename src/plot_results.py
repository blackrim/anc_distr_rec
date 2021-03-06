import sys
import os
from ete3 import Tree, TreeStyle, TextFace,faces, AttrFace

def mylayout(node):
    nf = faces.ImgFace(conf.outdir+node.name+".png")
    nf.margin_bottom = 20
    nf.margin_right = 20
    nf.margin_left = 20
    faces.add_face_to_node(nf,node,0,position = "branch-top")
    if node.is_leaf():
        ff = AttrFace("name",fsize=30)
        ff.margin_left = 20
        faces.add_face_to_node(ff, node, column=1)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("python "+sys.argv[0]+" treefile dir")
        sys.exit(0)
    inf = open(sys.argv[1],'r')
    treestring = inf.readline().strip()
    inf.close()
    t = Tree(treestring,format=1)

    # Basic tree style
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.layout_fn = mylayout
    # Add two text faces to different columns
    t.show(tree_style=ts)
