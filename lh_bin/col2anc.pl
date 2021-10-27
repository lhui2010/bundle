#!/usr/bin/perl -w

# !Important!
# filename fasta, gff, and gene name must be suffixed with species abbreviation.
# like chr01_arath. arath.peach.collinearity
# Perl script for converting MCScanX's collinearity file to MCScan(jcvi)'s anchor file

# ############### Parameters ###############
# # MATCH_SCORE: 50
# # MATCH_SIZE: 5
# # GAP_PENALTY: -1
# # OVERLAP_WINDOW: 5
# # E_VALUE: 1e-05
# # MAX GAPS: 25
# ############### Statistics ###############
# # Number of collinear genes: 1155, Percentage: 5.31
# # Number of all genes: 21751
# ##########################################
# ## Alignment 0: score=277.0 e_value=6.4e-10 N=6 GWHBEBS00000009_Arfim&GWHBEBS00000009_Arfim minus
#  0-  0:        AfUnG005400.1_Arfim     AfUnG009800.1_Arfim       4e-12
#  0-  1:        AfUnG005600.1_Arfim     AfUnG009500.1_Arfim       1e-36
#  0-  2:        AfUnG005700.1_Arfim     AfUnG009200.1_Arfim       6e-22
#  0-  3:        AfUnG006400.1_Arfim     AfUnG009100.1_Arfim       4e-33
#  0-  4:        AfUnG007400.1_Arfim     AfUnG008600.1_Arfim       6e-70
#  0-  5:        AfUnG007900.1_Arfim     AfUnG008100.1_Arfim       8e-31


sub extract_ortho_from_line;

# last five charactor of left species. defined once on first line 
my $filename = $ARGV[0];
my $left_suffix=$filename;
$left_suffix =~s/\..*//;
my ($genera, $species) = split/_/, $left_suffix;
$left_suffix = substr($genera, 0, 2).substr($species, 0,3);


my @buff_list;
# where this is a reverse alignment
my $minus_flag = 0;
# where left and right species and reversed
my $reversed_sp_found=0;
while(<>)
{
    my @e=split;
    if(/## Alignment/)
    {
        my $buff = join("", @buff_list);
        print $buff;
        print '#'.$_;

        @buff_list = ();
        if(/minus/)
        {
            $minus_flag = 1;
        }
        else
        {
            $minus_flag = 0;
        }
        my @header = split;
        my $chr_ortho = $header[-2];
        $chr_ortho=~s/&.*//;
        my $this_left_suffix = substr($chr_ortho, -5, 5);
#        if($left_suffix eq "")
#        {
#            $left_suffix = $this_left_suffix;
#            $reversed_sp_found = 0;
#        }
#        else
#        {
        if($left_suffix eq $this_left_suffix)
        {
            $reversed_sp_found = 0;
        }
        else
        {
            $reversed_sp_found = 1;
            #    print STDERR "$left_suffix\t".$this_left_suffix;
        }
    } 
    elsif(/#/)
    {
        next;
    }
    else
    { 
        my $ortho_line = extract_ortho_from_line($_, $reversed_sp_found);
        if($reversed_sp_found && $minus_flag)
        {
            unshift @buff_list, $ortho_line;
        }
        else
        {
            push @buff_list, $ortho_line;
        }
    }
}

my $buff = join("", @buff_list);
print $buff;

sub extract_ortho_from_line()
{
    my ($line, $reverse) = @_;
    $line=~s/.*?:\s+//;
    my @list = split/\s+/, $line;
    my $score = pop @list;
    my ($orthoA, $orthoB) = @list;
    if($reverse)
    {
        ($orthoA, $orthoB) = ($orthoB, $orthoA);
    }
    my $result = $orthoA."\t".$orthoB."\t".$score."\n";
    return $result;
}

