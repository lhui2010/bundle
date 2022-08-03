#!/usr/bin/env perl

my @f = @ARGV;
# Nissolia_schottii.Zenia_insignis.indWGD
# Nissolia_schottii.Zenia_insignis.comWGD
#

my %hash;


my @sp_list = ("Medicago_truncatula",
"Melilotus_albus",
"Trifolium_subterraneum",
"Pisum_sativum",
"Cicer_arietinum",
"Astragalus_sinicus",
"Glycyrrhiza_uralensis",
"Lotus_japonicus",
"Vigna_unguiculata",
"Phaseolus_lunatus",
"Lablab_purpureus",
"Glycine_soja",
"Amphicarpaea_edgeworthii",
"Cajanus_cajan",
"Spatholobus_suberectus",
"Abrus_precatorius",
"Arachis_duranensis",
"Aeschynomene_evenia",
"Dalbergia_odorifera",
"Nissolia_schottii",
"Lupinus_albus",
"Ammopiptanthus_nanus",
"Cladrastis_platycarpa",
"Styphnolobium_japonicum",
"Castanospermum_australe",
"Dipteryx_alata",
"Chamaecrista_pumila",
"Senna_septemtrionalis",
"Faidherbia_albida",
"Mimosa_pudica",
"Zenia_insignis",
"Duparquetia_orchidacea",
"Eperua_falcata",
"Sindora_glabra",
"Lysidice_rhodostegia",
"Bauhinia_variegata",
"Cercis_chinensis");

my %subfam_hash = ("Medicago_truncatula" => "Papilionoideae", 
"Melilotus_albus" => "Papilionoideae", 
"Trifolium_subterraneum" => "Papilionoideae", 
"Pisum_sativum" => "Papilionoideae", 
"Cicer_arietinum" => "Papilionoideae", 
"Astragalus_sinicus" => "Papilionoideae", 
"Glycyrrhiza_uralensis" => "Papilionoideae", 
"Lotus_japonicus" => "Papilionoideae", 
"Vigna_unguiculata" => "Papilionoideae", 
"Phaseolus_lunatus" => "Papilionoideae", 
"Lablab_purpureus" => "Papilionoideae", 
"Glycine_soja" => "Papilionoideae", 
"Amphicarpaea_edgeworthii" => "Papilionoideae", 
"Cajanus_cajan" => "Papilionoideae", 
"Spatholobus_suberectus" => "Papilionoideae", 
"Abrus_precatorius" => "Papilionoideae", 
"Arachis_duranensis" => "Papilionoideae", 
"Aeschynomene_evenia" => "Papilionoideae", 
"Dalbergia_odorifera" => "Papilionoideae", 
"Nissolia_schottii" => "Papilionoideae", 
"Lupinus_albus" => "Papilionoideae", 
"Ammopiptanthus_nanus" => "Papilionoideae", 
"Cladrastis_platycarpa" => "Papilionoideae", 
"Styphnolobium_japonicum" => "Papilionoideae", 
"Castanospermum_australe" => "Papilionoideae", 
"Dipteryx_alata" => "Papilionoideae", 
"Chamaecrista_pumila" => "Caesalpinioideae", 
"Senna_septemtrionalis" => "Caesalpinioideae", 
"Faidherbia_albida" => "Caesalpinioideae", 
"Mimosa_pudica" => "Caesalpinioideae", 
"Zenia_insignis" => "Dialioideae", 
"Duparquetia_orchidacea" => "Duparquetioideae", 
"Eperua_falcata" => "Detarioideae", 
"Sindora_glabra" => "Detarioideae", 
"Lysidice_rhodostegia" => "Detarioideae", 
"Bauhinia_variegata" => "Cercidoideae", 
"Cercis_chinensis" => "Cercidoideae"); 

# %subfam_hash = ("Medicago_truncatula" => "grey0", 
# "Melilotus_albus" => "grey0", 
# "Trifolium_subterraneum" => "grey0", 
# "Pisum_sativum" => "grey0", 
# "Cicer_arietinum" => "grey0", 
# "Astragalus_sinicus" => "grey0", 
# "Glycyrrhiza_uralensis" => "grey0", 
# "Lotus_japonicus" => "grey0", 
# "Vigna_unguiculata" => "grey0", 
# "Phaseolus_lunatus" => "grey0", 
# "Lablab_purpureus" => "grey0", 
# "Glycine_soja" => "grey0", 
# "Amphicarpaea_edgeworthii" => "grey0", 
# "Cajanus_cajan" => "grey0", 
# "Spatholobus_suberectus" => "grey0", 
# "Abrus_precatorius" => "grey0", 
# "Arachis_duranensis" => "grey0", 
# "Aeschynomene_evenia" => "grey0", 
# "Dalbergia_odorifera" => "grey0", 
# "Nissolia_schottii" => "grey0", 
# "Lupinus_albus" => "grey0", 
# "Ammopiptanthus_nanus" => "grey0", 
# "Cladrastis_platycarpa" => "grey0", 
# "Styphnolobium_japonicum" => "grey0", 
# "Castanospermum_australe" => "grey0", 
# "Dipteryx_alata" => "grey0", 
# "Chamaecrista_pumila" => "cyan2", 
# "Senna_septemtrionalis" => "cyan2", 
# "Faidherbia_albida" => "cyan2", 
# "Mimosa_pudica" => "cyan2", 
# "Zenia_insignis" => "dodgerblue1", 
# "Duparquetia_orchidacea" => "green3", 
# "Eperua_falcata" => "magenta2", 
# "Sindora_glabra" => "magenta2", 
# "Lysidice_rhodostegia" => "magenta2", 
# "Bauhinia_variegata" => "gold1", 
# "Cercis_chinensis" => "gold1"); 
for my $f(@f)
{
    my @e=split/\./, $f;
    #print $e[-1];
    my $count = `wc -l $f`;
    $count =~s/\s.*//;
    chomp $count;
    if($e[-1] eq "indWGD")
    {
        $hash{$e[0]}{$e[1]}{"ind"} = $count;
    }
    elsif($e[-1] eq "comWGD")
    {
        $hash{$e[0]}{$e[1]}{"com"} = $count;
        #    print "COMMON\n\n";
    }
    else
    {
        print $e[-1];
    }
}

print("QRY\tREF\tCount\tTreeType\tsubfam\n");
for my $k(sort keys %hash)
{
    #for my $l (sort keys %{$hash{$k}})
    #print $sp_list[0];
    for my $l (@sp_list)
    {
        #print $l;
        if(exists $hash{$k}{$l})
        {
            print($k, "\t", $l, "\t", $hash{$k}{$l}{"ind"}, "\tind\t", $subfam_hash{$l}, "\n");
            print($k, "\t", $l, "\t", $hash{$k}{$l}{"com"}, "\tcom\t", $subfam_hash{$l}, "\n");
        }
    }
}




