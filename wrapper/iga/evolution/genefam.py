"""
Function of this scripts
"""

import logging
import coloredlogs

from iga.apps.base import emain

logger = logging.getLogger(__name__)
coloredlogs.install(level='DEBUG', logger=logger)

"""

# Step1. Install
git clone https://github.com/hahnlab/CAFE5.git
mamba create -n cafe -c bioconda cafe

#Step 2. Transforming orthofinder result to cafe input
sed "s/^/(null)\t/;s/^.*Orthogroup/Desc\tFamily ID/" Orthogroups.GeneCount.tsv > cafe.genefam.txt

#Step 3. Running cafe5

cafe5 -i mammal_gene_families.txt -t mammal_tree.txt -k 2

# Input format

1. mammal_gene_families.txt need to be type of  unix, use :set ff=unix to change it

# Example 
Example: mammal_gene_families.txt (tab seprated for different column)
Desc	Family ID	human	chimp	orang	baboon	gibbon	macaque	marmoset rat	mouse	cat	horse	cow
ATPase	ORTHOMCL1	 52	 55	 54	 57	 54	  56	  56	 53	 52	57	55	 54
(null)	ORTHOMCL2	 76	 51	 41	 39	 45	  36	  37	 67	 79	37	41	 49
Example: tree (The newick tree, nothing else
((((cat:68.710507,horse:68.710507):4.566782,cow:73.277289):20.722711,(((((chimp:4.444172,human:4.444172):6.682678,orang:11.126850):2.285855,gibbon:13.412706):7.211527,(macaque:4.567240,baboon:4.567240):16.056992):16.060702,marmoset:36.684935):57.315065):38.738021,(rat:36.302445,mouse:36.302445):96.435575);




########################################## 
# https://zhuanlan.zhihu.com/p/413101212
# 查看文件内容
#第一列是使用的时间分歧树内容
Tree:(Vvinifera:1.0176,(Fvesca:1.01265,(Csinensis:0.888532,((Dzibethinus:0.26756,Tcacao:0.26756):0.517603,Cpapaya:0.785163):0.103369):0.124118):0.004947)ivaana
#推断的出生-死亡参数λ预测值
Lambda: 0.384269
Lambda tree:    (1,(1,(1,((1,1)1,1)1)1)1)
#每个节点编号
# IDs of nodes:(Vvinifera<0>,(Fvesca<2>,(Csinensis<4>,((Dzibethinus<6>,Tcacao<8>)<7>,Cpapaya<10>)<9>)<5>)<3>)ivaana<1>
# Output format for: ' Average Expansion', 'Expansions', 'No Change', 'Contractions', and 'Branch-specific P-values' = (node ID, node ID): (0,3) (2,5) (4,9) (6,8) (7,10)
# Output format for 'Branch cutting P-values' and 'Likelihood Ratio Test': (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#平均每个基因家族中扩张的基因数目，负数表示基因家族收缩
Average Expansion:      (-0.0456955,-0.0292223) (1.53621,0)     (0.201884,0)    (0.446966,0.0637459)    (0,0.000175246)
#发生了扩张的基因家族数目
Expansion :     (3777,0)        (11643,0)       (3901,0)        (8934,3795)     (0,3380)
#没有发生改变的基因家族数目
nRemain :       (10365,22158)   (6051,22825)    (11427,22825)   (6762,12401)    (22825,11162)
#发生了收缩的基因家族数目
nDecrease :     (8683,667)      (5131,0)        (7497,0)        (7129,6629)     (0,8283)
#每个基因家族具体扩张和收缩的情况，ID、树的信息（物种下划线后是该基因家族的数量）、基因家族扩张收缩的pvalue、每个分枝的基因家族扩张的显著性
'ID'    'Newick'        'Family-wide P-value'   'Viterbi P-values'      'cut P-value'   'Likelihood Ratio'
OG0000000       (Vvinifera_0:1.0176,(Fvesca_3:1.01265,(Csinensis_6:0.888532,((Dzibethinus_0:0.26756,Tcacao_1:0.26756)_0:0.517603,Cpapaya_374:0.785163)_0:0.103369)_0:0.124118)_0:0.004947)ivaana_1      0.5005  ((-,-),(-,-),(-,-),(-,-),(-,-))

# 提取结果
#使用python2运行，如果环境没python2可以安装。
python2 $CafePATH/python_scripts/cafetutorial_report_analysis.py -i out.cafe  -o out.report
#形成文件
#out.report_anc.txt #和GeneCounts.tsv内容相似，格式不一样
#out.report_fams.txt #每个节点发生扩张收缩的基因家族，有*是具有显著性。
#out.report_pub.txt #每个物种收缩、扩张、获得、丢失的情况
#out.report_node.txt #每个节点基因家族收缩、扩张和快速进化的家族数

# 画树
# 
# 
# 
# #也是使用安装包内的python脚本
# #但如果使用的不是具有界面化的服务器，建议按照以下操作修改脚本
# vi $CafePATH/python_scripts/cafetutorial_draw_tree.py
# ## 在import 模块后
# ## 添加一行不输出图片在桌面
# plt.switch_backend('agg')
# #运行，也是用python2
# python2 $CafePATH/python_scripts/cafetutorial_draw_tree.py \
# -i out.report_node.txt  -y Expansions  \
# -t '(Vvinifera:1.0176,(Fvesca:1.01265,(Csinensis:0.888532,((Dzibethinus:0.26756,Tcacao:0.26756):0.517603,Cpapaya:0.785163):0.103369):0.124118):0.004947)'  \
# -d '(Vvinifera<0>,(Fvesca<2>,(Csinensis<4>,((Dzibethinus<6>,Tcacao<8>)<7>,Cpapaya<10>)<9>)<5>)<3>)<1>' \
# -o Expansionstree.png
-i 输入的信息文件
-y是对应输入文件的标题选择展示 可选：Expansions/Contractions/Rapid
-t 是输入树文件 （在cafe输出文件能找到）
-d 是树的结构文件（在cafe输出文件能找到）
-o 命名，脚本默认是png；
不过由于是由python来写，还是能修改脚本的命令调整输出为svg



#将这句的格式调整svg输出或者pdf都可
fig.savefig(output_file, format='svg', bbox_inches='tight', dpi=300)

https://zhuanlan.zhihu.com/p/413101212

"""

def abc(foo=None, bar=None):
    """
    Remember:
     1. Positional args are recognized by software with None value
     2. All Non-positional args need to be given a default value, like ''
     2. no % is allowed in function docs.
    Returns:
    """

if __name__ == "__main__":
    emain()
