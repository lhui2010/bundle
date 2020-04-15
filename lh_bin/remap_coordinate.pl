#!/lustre/user/liuhui/bin/Assembly/software/perl/bin/perl 

use 5.12.0;
use Bio::FeatureIO;#THE MODULE GFF.pm was modified to bypass iregular Terms
use Bio::SeqFeature::Annotated;
use Bio::SeqIO;

use warnings;
use Data::Dumper;

use Bio::RangeI;
use Bio::Coordinate::GeneMapper;

use Time::HiRes;


#compare the derived allele frequency based on concensus genotype of the four populations
#Tue Jan  8 12:48:22 CST 2013


# vars sub-routine needed
#I hate to add those vars here, but if evolving into SNP_stat is too hard


our @DescriptionForCodonStat = (qw/Synonymous NonSynonymous/, ('Radical')x15);


#use Bio::Tools::CodonTable;
#use Bio::Align::DNAStatistics;


#my $myCodonTable   = Bio::Tools::CodonTable->new();


if(@ARGV <2)
{
	die "usage:  $0.pl  gff_of_superscaffolds gff_to_be_transformed \n";
}


my $gff=Bio::FeatureIO->new(-file => "$ARGV[0]", -format => 'GFF', -version => 3);
#my $fa = Bio::SeqIO->new(-file => "$ARGV[1]", -format => 'Fasta');
open my $gff_file, $ARGV[1] or die;
#open my $out_file, ">$ARGV[2].stat.via_ancestor" or die;
#open $snp, "chromosome01.genotype.snp" or die;

#warn "Reading fasta...";
#open $gff_file, "tmp.snp" or die;
#warn "Done\n";



#######################read fasta##########################
#my %seq;
#while(my $tmp = $fa->next_seq())
#{
#	        $seq{$tmp->display_id()} = $tmp;
#}


warn "Reading gff...\n";



######################read gff into memory##################


##THESE TWO VARs IS ALSO USED IN SUBROUTINE AS GLOBAL VAR
#like a index
##THESE TWO VARs IS ALSO USED IN SUBROUTINE AS GLOBAL VAR
#like a index
our %array_of_gene_start;#used for determin the loci's nearest gene, sorted by ascending;
our %start2end;#return end of the gene given its start and chr
our %start2gene;#return the gene of this start

my $chr;
my %gene; #containing total genes from this gff file
while ( my $feature = $gff->next_feature() ) {


#       no strict 'refs';
#       print "Instance METHOD IS  " . Dumper( \%{ref ($seq)."::" }) ;exit;
        if($feature->type ->name eq "chromosome")
        {
                $start2gene{$feature->seq_id}{$feature -> start} = $feature -> get_Annotations ('Name')->value;
                $start2end{$feature->seq_id}{$feature -> start}  = $feature -> end;


        }

        next unless($feature-> type ->name eq "scaffolds");

        push  @{$gene{($feature -> get_Annotations ('Name'))->value}}, $feature;

	$chr = $feature->seq_id;

}
for my $sid( %start2gene)
{
        $array_of_gene_start{$sid} = [sort {$a<=>$b} (keys %{$start2gene{$sid}})];
}

warn "Done\nBuilding gene mapper...";


###################Building Gene Mapper#####################

use Storable 'dclone';

my (%chr_mapper, %cds_mapper);

for my $this_gene (sort keys  %gene)
{
	$chr_mapper{$this_gene} = Bio::Coordinate::GeneMapper->new(
		-in  => 'chr',
		-out => 'cds',
		-exons => [
			map {
				Bio::Location::Simple->new(
					-start  => $_ ->start,#location->{_start}, 
					-end =>    $_ ->end, #location->{_end},
					-strand => $_ ->strand, #location->{_strand},
					-seq_id => $_ ->seq_id,
					)
			} @{$gene{$this_gene}}
			],
		);
	$cds_mapper{$this_gene} =  dclone($chr_mapper{$this_gene});
	$cds_mapper{$this_gene} -> swap;

#	$chr = $this_gene;
#	$hash_return{"strand"} = $$chr_mapper{$this_gene} ->map($loci_chr) ->strand;
}



warn "Done\nReading SNP file...\n\n\n";
##############read SNP file#############
#print $_, "\t" for (qw/chromosome strand location_chr ref_base snp_base location_cds phase ref_codon snp_codon mutation_type/);
#print "\n";


#eg of snp file
#indica japonica nivara rufipogon
#chromosome01    1       C       C8 G1      C9      C4      C3
my $i;
my %frequency;

my %strand_hash = (
	"+" => 1,
	"-" => -1,
	);
my %hash_strand = (
	1 =>"+",
	-1 => "-",
	);

while(<$gff_file>)
{
#scaffold78      GLEAN   mRNA    282169  288267  0.833564        +       .       ID=Olong01m10024664.1;
	chomp;
	my @e=split;
	my $scaf = $e[0];
	my $start = $e[3];
	my $end = $e[4];
	my $gene_strand = $strand_hash{$e[6]};
	
	my $loci =  Bio::Location::Simple->new(
	        -start => $start,
	        -end => $end,
		-strand => $gene_strand
				        );
	my $trans_form_loci = $cds_mapper{$scaf}->map($loci);
#	my $chr_strand = $trans_form_loci->{"strand"};
eval{
	$e[0] = $chr;
	$e[3] = $trans_form_loci->start;
	$e[4] = $trans_form_loci->end;
	$e[6] = $hash_strand{$trans_form_loci->strand};
};
	$_ = join"\t", @e;
	$_.="\n";
	print;

#	print Data::Dumper->Dump($loci);
#	print Data::Dumper->Dump($trans_form_loci);
#	exit;

}#end reading snp

warn "Done\nProgram finished..Outputing result into $ARGV[2].stat\n";

warn "Clearing memory...\n";

#############sub routine###############


#accept: chromosome_id  coordinate_on_chr Ref-base SNP-base \%chr_mapper, \%cds_mapper
#return an a hash containing the following elements
#chromosome strand location_chr ref_base snp_base location_cds phase ref_codon snp_codon mutation_type





