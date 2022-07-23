require('ape')

z <- "((Aquilegia_coerulea:0.363336,(((Averrhoa_carambola:0.229881,Arabidopsis_thaliana:0.104743)N4:0.0597203,(Populus_trichocarpa:0.254181,((((Cercis_chinensis:0.0392068,Bauhinia_variegata:0.193284)N10:0.0834982,((Goniorrhachis_marginata_Det:0.0692086,Lysidice_rhodostegia:0.059051)N13:0.00446505,(Eperua_falcata:0.0647717,Sindora_glabra:0.140978)N14:0.0256284)N11:0.111491)N8:0.00774706,(Duparquetia_orchidacea:0.306583,(((((Castanospermum_australe:0.0597794,Angylocalyx_braunii_Pap:0.105806)N22:0.0170228,Dipteryx_alata:0.120242)N19:0.00644757,((Cladrastis_platycarpa:0.050087,Styphnolobium_japonicum:0.0661515)N23:0.015931,(Zollernia_splendens_Pap:0.065509,(((((Dalbergia_odorifera:0.0926784,Aeschynomene_evenia:0.123786)N36:0.0273695,Arachis_duranensis:0.176067)N33:0.00668278,Nissolia_schottii:0.117445)N30:0.052343,(Ammopiptanthus_nanus:0.0735248,Lupinus_albus:0.195145)N31:0.0404402)N28:0.00328628,(Andira_inermis_Pap:0.0642007,(((((Cicer_arietinum:0.0910488,(((Melilotus_albus:0.0493188,Medicago_truncatula:0.0653878)N48:0.0302107,Trifolium_subterraneum:0.0936658)N46:0.00879787,Pisum_sativum:0.239168)N43:0.0173324)N41:0.0659768,Astragalus_sinicus:0.170014)N39:0.0163204,Glycyrrhiza_uralensis:0.0875438)N37:0.0192061,Lotus_japonicus:0.171052)N34:0.0188284,(Abrus_precatorius:0.191197,(Spatholobus_suberectus:0.0545558,(((Lablab_purpureus:0.0598149,(Phaseolus_lunatus:0.0576938,Vigna_unguiculata:0.0573969)N47:0.0113987)N44:0.0584329,(Amphicarpaea_edgeworthii:0.117141,Glycine_soja:0.0687862)N45:0.0203646)N42:0.0109456,Cajanus_cajan:0.100762)N40:0.00628519)N38:0.0393983)N35:0.0335326)N32:0.0245228)N29:0.00467139)N27:0.011817)N24:0.0166959)N20:0.0125437)N17:0.0471633,(Umtiza_listeriana_Cae:0.0680843,((Chamaecrista_pumila:0.123298,Senna_septemtrionalis:0.082298)N25:0.017063,(Mimosa_pudica:0.184517,Faidherbia_albida:0.0806019)N26:0.129488)N21:0.025528)N18:0.0404782)N15:0.00469256,(Zenia_insignis:0.0895134,Dialium_schlechtneri_Dia:0.0997121)N16:0.146227)N12:0.00772078)N9:0.00297183)N7:0.0829144,Polygala_tatarinowii:0.559894)N6:0.121425)N5:0.057966)N3:0.0521563,Vitis_vinifera:0.583657)N2:0.142019)N1:0.134598,Oryza_sativa:0.134598)N0;"
tz <- read.tree(text = z)
#http://blog.phytools.org/2012/02/printing-tree-to-file-without-branch.html
#pdf("tree_with_node_labels.pdf", width=28, height=14)
#plot(tz, use.edge.length = FALSE, main = "Without branch lengths", show.node.label=T)
#nodelabels()
#edgelabels()
#dev.off()


#a<- read.table('gene_pav.nodulation.txt', header = T)
a<- read.table('gene_pav.txt', header = T)
b<- a[a$OrthoType=="single",]
c<- na.omit(b)

#https://stackoverflow.com/questions/37575081/r-for-loop-for-all-groups-of-rows-with-the-same-value-in-column-do
for(i in unique(c$Gene)){
  # chek if there is a one for 
  # each group
  f <- c[c$Gene==i,]
  mrca_node <- getMRCA(tz, f$Species)
  print(paste(i, tz$node.label[mrca_node-Ntip(tz)]), sep="\t")
}


