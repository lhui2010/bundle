
from iga.apps.base import sh, emain
from collections import defaultdict
import os
import re
import pandas as pd
import logging

#Medicago_truncatula__v__Averrhoa_carambola.tsv.m-to-m.nodgenes
def get_gene_age(dir=None, topology_file=None):
    """
    """
    gene_dict = {}
    singlecopy_dict = defaultdict(dict)
    duplicate_dict  = defaultdict(dict)
    cmd = r"""
    cd {}
    for i in *tsv
    do
        for j in 1-to-1 1-to-m m-to-1 m-to-m
        do
            perl ortholog2ortho.$j.pl $i > $i.$j    
            #perl selectItem.pl -k nod_genes.txt.addinfo Medicago_truncatula__v__$i.tsv.$j > Medicago_truncatula__v__$i.tsv.$j.nodgenes
        done
    done
    """.format(dir)
    single_extensions = ['tsv.1-to-1', 'tsv.1-to-m']
    duplicate_extensions = ['tsv.m-to-1']
#, 'm-to-m']
    #sh(cmd)
# Step1. Prepare input
    for f in os.listdir(dir):
        for extension in single_extensions:
            qryname = re.sub("\..*", "", f)
            qryname = re.sub(".*__v__", "", qryname)
            if f.endswith(extension):
                logging.info(f)
                f_content = pd.read_table(f, sep='\t', index_col = 1)
                singlecopy_dict[qryname] = f_content[qryname].to_dict()
                for ok in singlecopy_dict[qryname]:
                    #record gene names
                    gene_dict[ok] = 1
                continue
                #print(singlecopy_dict)
                #exit(1)
        for extension in duplicate_extensions:
            qryname = re.sub("\..*", "", f)
            qryname = re.sub(".*__v__", "", qryname)
            if f.endswith(extension):
                logging.info(f)
                f_content = pd.read_table(f, sep='\t', index_col = 1)
                raw_dict = f_content[qryname].to_dict()
# A,B C
                for old_k in raw_dict:
                    for ok in old_k.split(','):
                        duplicate_dict[qryname][ok] = raw_dict[old_k]
                        gene_dict[ok] = 1
# Step2. Get pairs from topology
#Oryza_sativa,Aquilegia_coerulea,|Vitis_vinifera,Arabidopsis_thaliana,Averrhoa_carambola,Populus_trichocarpa,|Polygala_tatarinowii,|Cercis_chinensis,Bauhinia_variegata,Lysidice_rhodostegia,Goniorrhachis_marginata,Sindora_glabra,Eperua_falcata|,Duparquetia_orchidacea,|Dialium_schlechtneri,Zenia_insignis,|Umtiza_listeriana,Mimosa_pudica,Faidherbia_albida,Senna_septemtrionalis,Chamaecrista_pumila,|Dipteryx_alata,Castanospermum_australe,Angylocalyx_braunii,|Styphnolobium_japonicum,Cladrastis_platycarpa,|Zollernia_splendens|,Ammopiptanthus_nanus,Lupinus_albus,Nissolia_schottii,Dalbergia_odorifera,Aeschynomene_evenia,Arachis_duranensis,|Andira_inermis|,Abrus_precatorius,Spatholobus_suberectus,Cajanus_cajan,Amphicarpaea_edgeworthii,Glycine_soja,Lablab_purpureus,Phaseolus_lunatus,Vigna_unguiculata,|Lotus_japonicus,Glycyrrhiza_uralensis,Astragalus_sinicus,Cicer_arietinum,Pisum_sativum,Trifolium_subterraneum,Melilotus_albus,Medicago_truncatula
    with open(topology_file) as fh:
        topo = fh.readline().rstrip()
    topo_list = topo.split('|')
    prefix_group = topo_list[0]
    for i in range(1, len(topo_list)):
        prefix_group_list = prefix_group.split(',')
        basal_group = topo_list[i]
        basal_group_list = basal_group.split(',')

        line_tag = "{}|{}".format(prefix_group, basal_group)
        for gene in gene_dict:
            flag = 0
#Find duplicate in prefix but become 1-to-1 ortholog in basal groups. i.e. Find it's duplication date
            for pg in prefix_group_list:
                if gene in duplicate_dict[pg]:
                    flag += 1
                    break
            for bg in basal_group_list:
                if gene in singlecopy_dict[bg]:
                    flag += 10
                    break
            if flag == 11:
                print("{}\t{}".format(line_tag, gene))
        prefix_group += basal_group


if __name__ == "__main__":
    emain()
