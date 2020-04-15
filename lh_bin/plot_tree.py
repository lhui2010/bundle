#!/usr/bin/env python
from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle
import sys
import re   
#t=Tree('/lustre/home/liuhui/project/buzzo/interproscan/orthofinder/total_tps.aln.format.nwk')
t=Tree(sys.argv[1])
ts = TreeStyle()
ts.show_leaf_name=True
ts.show_branch_support=True
#ts.mode = "c"

ts.scale=2000

#species_list = ['B73', 'W22', 'A188','Mo17', 'PH207', 'sorghum',]
#color_list = ['#ff3f7c', '#ffa7fc', '#ffe2f4', '#7fff3f', '#3fbcff', '#ff7f3f']

nstyle_inner = NodeStyle()
nstyle_inner['size'] = 0
nstyle_inner['fgcolor'] = 'black'

nstyle = NodeStyle()
#    nstyle[species_list[i]]["hz_line_color"] =  color_list[i]
#    nstyle[species_list[i]]["fgcolor"] = color_list[i]
nstyle["shape"] = 'circle'

#nstyle["bgcolor"]="Moccasin"


for node in t.traverse():
    if node.is_leaf():
        sp=re.match(r'(.*?)_', node.name)
        if(sp):
            node.set_style(nstyle)
    else: 
        node.set_style(nstyle_inner)
#            N = AttrFace("name", fsize=5, fgcolor=nstyle[sp[1]]['fgcolor'])
#            node.add_face(N, 0, position = "aligned")      

#t.show(tree_style=ts)
t.render(sys.argv[1] + ".pdf", tree_style = ts)

