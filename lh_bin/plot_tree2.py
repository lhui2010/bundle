from ete3 import Tree, faces, AttrFace, TreeStyle, NodeStyle, TextFace, add_face_to_node
import sys
import re   
t=Tree(sys.argv[1])
ts = TreeStyle()
#ts.show_leaf_name=True
ts.mode = "c"
ts.show_leaf_name = False
ts.legend.add_face(TextFace("0.5 support"), column=1)

#color_list = ['#ff3f7c', '#7fff3f', '#3fbcff']
color_list = ['#fbb4ae', '#b3cde3', '#ccebc5', '#decbe4', '#fed9a6', '#ffffcc']

def my_layout(node):
    color_list = ['#fbb4ae', '#b3cde3', '#ccebc5', '#decbe4', '#fed9a6', '#ffffcc']
    key_list = [ r'LOC_Os10g34500|LOC_Os05g41530|LOC_Os12g05920|AK069174|LOC_Os01g14840|LOC_Os09g27650|LOC_Os08g36390|LOC_Os03g13400|LOC_Os10g28330|LOC_Os01g09850|LOC_Os01g14010|LOC_Os08g44050|LOC_Os09g38340|LOC_Os02g45054|LOC_Os04g47860|LOC_Os01g39110|LOC_Os01g70870|LOC_Os02g31890|LOC_Os07g39310|LOC_Os03g10140|LOC_Os12g18150|LOC_Os01g68104|LOC_Os06g47840|LOC_Os10g18364|LOC_Os03g40710|LOC_Os04g39050', r'LOC_Os07g40080|LOC_Os01g62130|LOC_Os03g50850|LOC_Os05g02390|LOC_Os03g32220|LOC_Os09g38790|LOC_Os05g19970|LOC_Os09g39270|LOC_Os02g44130|LOC_Os04g46680|LOC_Os03g05690|LOC_Os06g44200|LOC_Os09g03500', r'AK110194|LOC_Os03g05490|LOC_Os06g48530|LOC_Os02g05610|LOC_Os11g47620|LOC_Os08g40560|LOC_Os07g39970|LOC_Os04g39520|AK110488', r'LOC_Os01g28230|LOC_Os11g25610|LOC_Os03g05480|LOC_Os08g19890|LOC_Os04g08060|LOC_Os02g57790|AK109981|LOC_Os08g20580|AK110095|LOC_Os02g02424|LOC_Os04g08600|LOC_Os09g39660|LOC_Os07g38240|LOC_Os09g25420|LOC_Os07g40300|LOC_Os12g39220|LOC_Os02g19180|LOC_Os09g25430|LOC_Os09g29750|LOC_Os03g11680|LOC_Os08g36110|LOC_Os09g27320|LOC_Os02g53530|LOC_Os03g13600|LOC_Os07g05900|LOC_Os03g60570|LOC_Os05g38620|LOC_Os02g10240|LOC_Os01g04120|LOC_Os01g32920|LOC_Os09g25440|LOC_Os02g01090|LOC_Os07g01180|LOC_Os05g20930|LOC_Os11g48000|LOC_Os04g50070|LOC_Os06g07020|LOC_Os01g62190|LOC_Os05g38600|LOC_Os02g44120|LOC_Os04g46670|LOC_Os01g66570|LOC_Os07g40950|LOC_Os02g08510|LOC_Os08g17640|LOC_Os06g10470|LOC_Os04g36650|LOC_Os03g60540|LOC_Os03g60560|LOC_Os11g47630|LOC_Os03g17150|LOC_Os03g55540|LOC_Os07g10870|LOC_Os12g39400|LOC_Os03g41390|LOC_Os05g14130|LOC_Os03g41110|LOC_Os08g44190|LOC_Os09g26050|LOC_Os09g26100|LOC_Os09g26200|LOC_Os06g49490|LOC_Os09g38610|LOC_Os03g57240|LOC_Os07g23450|LOC_Os05g05970|LOC_Os12g42250|AK119817|LOC_Os01g03840|LOC_Os06g51140|LOC_Os04g48520|LOC_Os12g38940|LOC_Os12g31840|LOC_Os01g68160|LOC_Os02g39580|LOC_Os02g34680|LOC_Os05g01550|LOC_Os09g31150|LOC_Os10g07470|LOC_Os01g01770|LOC_Os02g47920|LOC_Os03g32230|LOC_Os07g39960|LOC_Os10g40660|LOC_Os08g37904|LOC_Os08g37920|LOC_Os11g39580|LOC_Os11g39630|LOC_Os12g38960|LOC_Os04g08034|LOC_Os09g13630|AK119806|LOC_Os01g35040|LOC_Os11g38330', r'LOC_Os07g40780|LOC_Os07g44640|LOC_Os03g31240|LOC_Os02g57550|AK110195|LOC_Os02g35460|LOC_Os12g13130|LOC_Os04g08290|LOC_Os08g44830|LOC_Os12g07280|LOC_Os01g65080|LOC_Os03g62230|LOC_Os02g36360|LOC_Os09g13680|LOC_Os08g39390|LOC_Os06g40960|LOC_Os09g10980|LOC_Os01g63980|LOC_Os05g37190|LOC_Os10g35190|LOC_Os02g30180|LOC_Os09g21430|LOC_Os08g34310|LOC_Os11g34700|LOC_Os01g59450|LOC_Os11g30484|LOC_Os11g36930|LOC_Os01g40110|LOC_Os05g03020|LOC_Os03g61640|LOC_Os06g20020|AK110102|LOC_Os05g51830|LOC_Os01g62460|AK110182|LOC_Os06g22420|LOC_Os10g17740|LOC_Os01g57650|LOC_Os06g46910|LOC_Os11g06840|LOC_Os04g59380|LOC_Os04g02510|LOC_Os03g15790|LOC_Os08g16940|LOC_Os01g67970|LOC_Os03g49132|LOC_Os09g26210' ]
    if node.is_leaf():
         # If terminal node, draws its name
        face_color_index = -1
        name_face = AttrFace("name")
        for k in range(0,5):
            if re.match(key_list[k], node.name):
                face_color_index = k
                name_face.inner_background.color=color_list[face_color_index]
                #name_face.fgcolor=color_list[face_color_index]
                break
        faces.add_face_to_node(name_face, node, column=0, aligned = True, position="branch-right")
        #print( node.name + str(face_color_index))
    # Adds the name face to the image at the preferred position

ts.layout_fn = my_layout


nstyle_inner = NodeStyle()
nstyle_inner['size'] = 0
nstyle_inner['fgcolor'] = 'black'
#nstyle_inner["hz_line_width"] = 3
#nstyle_inner["vt_line_width"] = 3

nstyle = {}
for  i in range(0,6):
    nstyle[i] = NodeStyle()
#    nstyle[i]["hz_line_color"] =  color_list[i]
#    nstyle[i]["bgcolor"] = color_list[i]
#    nstyle[i]["node_bgcolor"] = color_list[i]
    nstyle[i]["shape"] = 'circle'
    nstyle[i]["size"] = 0
    nstyle[i]["hz_line_color"] = color_list[i]
    nstyle[i]["hz_line_width"] = 3
    nstyle[i]["vt_line_color"] = color_list[i]
    nstyle[i]["vt_line_width"] = 3

#nstyle['B73']["bgcolor"]="Moccasin"

color_count = 0

key=r'LOC_Os04g39050|LOC_Os09g03500|AK110488|LOC_Os11g38330'

#for node in t.traverse():
for node in t.traverse(strategy='preorder'):
#    print(node.name)
#    if node.is_leaf():
    if (not node.is_root() and not node.get_ancestors()[0].is_root()
        and (not node.get_ancestors()[0].get_ancestors()[0].is_root() or color_count == 2) ):
        node.set_style(nstyle[color_count])
        sp=re.match(key, node.name)
#        textface = TextFace(node.name)
#        if(re.match('LOC_Os07g05900', node.name)):
#            textface = TextFace(node.name, bold = True)
#        node.add_face(textface, 0, aligned = True)
        if(sp):
            print(sp)
            color_count +=1
    else: 
        node.set_style(nstyle_inner)
        if(not node.is_root() and not node.get_ancestors()[0].is_root()):
            flag = False
            for n in node.get_children()[0].get_children():
                if('AK110194' == n):
                    flag = True
            if(flag):
                node.set_style(nstyle[color_count])
#            N = AttrFace("name", fsize=5, fgcolor=nstyle[sp[1]]['fgcolor'])
#            node.add_face(N, 0, position = "aligned")      

#t.show(tree_style=ts)
t.render("mytree.pdf", tree_style = ts)

