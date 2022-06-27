import os
from ete3 import Tree
nwk_file=[i for i in os.listdir('.') if i.endswith('treefile')]


outgroup='Cechi'
spA='Caaus'
spB='Phlun'

com_print = ""
ind_print = ""
other_print = ""
for i in nwk_file:
    #path='tree/' + i
    path=i
    t=Tree(path)
    t_leave_name=[ii for ii in t.iter_leaf_names()]
    out_group_name=[ii for ii in t_leave_name if ii.endswith(outgroup)][0]####找到外类群######
    t.set_outgroup(out_group_name)####设置外类群，改变系统发育树的拓扑结构#####
    list_ind=[]
    list_com=[]
    for k in t.traverse():
        if k.is_leaf():
            continue
        else:
            child=k.get_children()[0].name + k.get_children()[1].name####获取节点名字#####
            if spA in child and spB in child:
                list_ind.append(child)
            elif spA in k.get_children()[0].name and spA in k.get_children()[1].name:
                list_com.append(child)
            elif spB in k.get_children()[0].name and spB in k.get_children()[1].name:
                list_com.append(child)
    if len(list_ind) ==2 :
        ind_print += t.write() + "\n"
    elif len(list_com) ==2 :
        com_print += t.write() + "\n"
    else:
        other_print += t.write() + "\n"

with open("{}.{}.{}.indWGD".format(outgroup,spA,spB), "w") as fh_ind, \
  open("{}.{}.{}.comWGD".format(outgroup,spA,spB), "w") as fh_com, \
  open("{}.{}.{}.otherWGD".format(outgroup,spA,spB), "w") as fh_other:
    fh_ind.write(ind_print)
    fh_com.write(com_print)
    fh_other.write(other_print)
