#!/usr/bin/env python
import csv
from collections import defaultdict
from pandas import DataFrame
import pandas as pd
import sys


#https://stackoverflow.com/questions/3269769/convert-three-column-text-file-to-matrix

input_eg = """
Abrus_precatorius       Aeschynomene_evenia     0.0426034775458634
Abrus_precatorius       Ammopiptanthus_nanus    0.0559874089113069
Abrus_precatorius       Amphicarpaea_edgeworthii        0.138162578272697
Abrus_precatorius       Arachis_duranensis      0.0480154924397725
Abrus_precatorius       Arachis_ipaensis        0.0636948241826619
Abrus_precatorius       Astragalus_sinicus      0.0796985704560489
Abrus_precatorius       Cajanus_cajan   0.141483699729247
Abrus_precatorius       Cercis_canadensis       0.565177073751075
Abrus_precatorius       Cercis_chinensis        0.547207578503591
Abrus_precatorius       Chamaecrista_fasciculata        -0.036514147682558
Abrus_precatorius       Chamaecrista_pumila     -0.0655096690833022
Abrus_precatorius       Cicer_arietinum 0.0753694144466679
Abrus_precatorius       Cicer_reticulatum       0.0889864190848149
Abrus_precatorius       Cladrastis_platycarpa   0.0137808184619186
"""
output_eg= """
qry     Abrus_precatorius       Aeschynomene_evenia     Ammopiptanthus_nanus    Amphicarpae
Abrus_precatorius               0.0426034775458634      0.0559874089113069      0.138162578
Aeschynomene_evenia     0.0422907461434441              0.0332686379916175      0.060979212
Ammopiptanthus_nanus    0.0564306943766211      0.0327144924337111              0.050354751
Amphicarpaea_edgeworthii        0.13859765623362        0.0609376460252394      0.050979614
Arachis_duranensis      0.0476858263494873      0.15016825082003        0.0315538069144682
Arachis_ipaensis        0.0637113088624425      0.186930350166827       0.0748319085746745
Astragalus_sinicus      0.0789823462270908      0.0644071458728237      0.0164877817217287
Cajanus_cajan   0.141655429068037       0.0603454377324713      0.0570601106824453      0.1
Cercis_canadensis       0.567323264425995       0.555122795463815       0.581443858106267
Cercis_chinensis        0.546652940654533       0.550721482249997       0.576780871830883
Chamaecrista_fasciculata        -0.0353384175850417     -0.0156072551307903     -0.04300891
Chamaecrista_pumila     -0.0659285208472855     -0.0776206466502875     -0.069030731118341
Cicer_arietinum 0.0761314231124298      0.0691930091965604      0.0394360806976979      0.0
Cicer_reticulatum       0.0888704966113206      0.0759917804993717      0.0528643303055777
"""

if len(sys.argv) <= 1:
    print("stat2matrix.py input > output\n####Input####\n{}\n####Output####\n{}".format(input_eg, output_eg))
    exit(1)

rdr = csv.reader(open(sys.argv[1]), delimiter='\t', skipinitialspace=True)
datacols = defaultdict(list)

# skip header
# rdr.next()
for qry, ref, result in rdr:
    datacols['qry'].append(qry)
    datacols['ref'].append(ref)
    datacols['result'].append(float(result))

df = DataFrame(datacols)
#https://stackoverflow.com/questions/13838405/custom-sorting-in-pandas-dataframe/13839029#13839029
df['qry'] = pd.Categorical(df['qry'], ["Abrus_precatorius", "Aeschynomene_evenia", "Ammopiptanthus_nanus", "Amphicarpaea_edgeworthii", "Arachis_duranensis", "Arachis_ipaensis", "Astragalus_sinicus", "Cajanus_cajan", "Cicer_arietinum", "Cicer_reticulatum", "Cladrastis_platycarpa", "Dalbergia_odorifera", "Glycine_max", "Glycine_soja", "Glycyrrhiza_uralensis", "Lablab_purpureus", "Lotus_japonicus", "Lupinus_albus", "Lupinus_angustifolius", "Medicago_polymorpha", "Medicago_ruthenica", "Medicago_truncatula", "Melilotus_albus", "Nissolia_schottii", "Phaseolus_acutifolius", "Phaseolus_lunatus", "Phaseolus_vulgaris", "Pisum_sativum", "Trifolium_subterraneum", "Vigna_angularis", "Vigna_radiata", "Vigna_subterranea", "Vigna_trilobata", "Vigna_unguiculata", "Chamaecrista_fasciculata", "Chamaecrista_pumila", "Faidherbia_albida", "Mimosa_pudica", "Senna_septemtrionalis", "Senna_tora", "Zenia_insignis", "Duparquetia_orchidacea", "Lysidice_rhodostegia", "Cercis_canadensis", "Cercis_chinensis"])

df['ref'] = pd.Categorical(df['ref'], ["Abrus_precatorius", "Aeschynomene_evenia", "Ammopiptanthus_nanus", "Amphicarpaea_edgeworthii", "Arachis_duranensis", "Arachis_ipaensis", "Astragalus_sinicus", "Cajanus_cajan", "Cicer_arietinum", "Cicer_reticulatum", "Cladrastis_platycarpa", "Dalbergia_odorifera", "Glycine_max", "Glycine_soja", "Glycyrrhiza_uralensis", "Lablab_purpureus", "Lotus_japonicus", "Lupinus_albus", "Lupinus_angustifolius", "Medicago_polymorpha", "Medicago_ruthenica", "Medicago_truncatula", "Melilotus_albus", "Nissolia_schottii", "Phaseolus_acutifolius", "Phaseolus_lunatus", "Phaseolus_vulgaris", "Pisum_sativum", "Trifolium_subterraneum", "Vigna_angularis", "Vigna_radiata", "Vigna_subterranea", "Vigna_trilobata", "Vigna_unguiculata", "Chamaecrista_fasciculata", "Chamaecrista_pumila", "Faidherbia_albida", "Mimosa_pudica", "Senna_septemtrionalis", "Senna_tora", "Zenia_insignis", "Duparquetia_orchidacea", "Lysidice_rhodostegia", "Cercis_canadensis", "Cercis_chinensis"])

df2 = df.pivot(index='qry', columns='ref', values='result')
df2 = df2.sort_values(axis = 0, by='qry')
df2 = df2.sort_values(axis = 1, by='ref')
df2.to_csv(sys.argv[1] + '.xls', sep='\t')

