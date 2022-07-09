#!/usr/bin/env perl

# ==> Cercis_chinensis.Bauhinia_variegata.1x2 <==
# 8,9,13,14,25,27,37,47,48,49,50,84,87,90,91,92,93,94,96,98,99,106,108,110,111,129,130,138,160,161,162,163,164,169,197,198,199,202,203,210,216,219,220,225,226,228,229,237,239,241,242,243,268,269,271,275,276,277,278,279,282,284,287,289,291,294,295,296,298,321,331,332,338,342,344,345,346,347,404,405,406,409,410,411,415,439,441,442,446,450,451,469,490,491,500,501,509,531,536,537,557,558,575,576,577,583,584
# 
# ==> Cercis_chinensis.Bauhinia_variegata.anchors <==
# ### Alignment 0: score=1929.0 e_value=9.3e-153 N=43 Chr10_Bavar&chr01_Cechi plus
# CECHI00012316-t1_Cechi	BAUVAR_g17813.t2_Bavar	1e-122
# CECHI00012318-t1_Cechi	BAUVAR_g17824.t1_Bavar	3e-33
# CECHI00012319-t1_Cechi	BAUVAR_g17825.t1_Bavar	7e-43
# CECHI00012321-t1_Cechi	BAUVAR_g17831.t1_Bavar	3e-41
# CECHI00012324-t2_Cechi	BAUVAR_g17839.t1_Bavar	2e-71
# CECHI00012327-t1_Cechi	BAUVAR_g17842.t1_Bavar	1e-30
# CECHI00012332-t1_Cechi	BAUVAR_g17851.t1_Bavar	5e-36
# CECHI00012333-t1_Cechi	BAUVAR_g17855.t1_Bavar	0
# CECHI00012336-t2_Cechi	BAUVAR_g17859.t1_Bavar	1e-152

my $ID = shift;

open IN, $ID  or die;
my %key_blocks;
while(<IN>)
{
    chomp;
    my @e=split/,/, $_;
    for my $e(@e)
    {
        $key_blocks{$e} = 1;
    }
}

close IN;

my $block_id = -1;
while(<>)
{
    if(/###/)
    {
        $block_id+=1;
        $block_hash{$block_id} .=$_;
    }
    else
    {
        $block_hash{$block_id} .=$_;
    }
}

for my $k (sort {$a<=>$b} keys %key_blocks)
{
    print $block_hash{$k};
}

