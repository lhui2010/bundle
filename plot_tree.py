from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle
import re   
t=Tree('/lustre/home/liuhui/project/buzzo/interproscan/orthofinder/total_tps.aln.format.nwk')
ts = TreeStyle()
ts.show_leaf_name=True
#ts.mode = "c"

species_list = ['B73', 'W22', 'A188','Mo17', 'PH207', 'sorghum',]
color_list = ['#ff3f7c', '#ffa7fc', '#ffe2f4', '#7fff3f', '#3fbcff', '#ff7f3f']

nstyle_inner = NodeStyle()
nstyle_inner['size'] = 0
nstyle_inner['fgcolor'] = 'black'

nstyle = {}
for  i in range(0,6):
    nstyle[species_list[i]] = NodeStyle()
    nstyle[species_list[i]]["hz_line_color"] =  color_list[i]
    nstyle[species_list[i]]["fgcolor"] = color_list[i]
    nstyle[species_list[i]]["shape"] = 'circle'

nstyle['B73']["bgcolor"]="Moccasin"


for node in t.traverse():
    if node.is_leaf():
        sp=re.match(r'(.*?)_', node.name)
        if(sp):
            node.set_style(nstyle[sp[1]])
    else: 
        node.set_style(nstyle_inner)
#            N = AttrFace("name", fsize=5, fgcolor=nstyle[sp[1]]['fgcolor'])
#            node.add_face(N, 0, position = "aligned")      

#t.show(tree_style=ts)
t.render("mytree.pdf", tree_style = ts)

