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
use Bio::Align::RadicalChanges;#DIY PACKAGE, cloned DNAStatistics except get_syn_changes


our %CodonStat = Bio::Align::RadicalChanges ->new()->get_syn_changes;#$CodonStat{'ATG'}{'ATG'} returns 0, after Description becomes Synony..
our @DescriptionForCodonStat = (qw/Synonymous NonSynonymous/, ('Radical')x15);


#use Bio::Tools::CodonTable;
#use Bio::Align::DNAStatistics;


#my $myCodonTable   = Bio::Tools::CodonTable->new();


if(@ARGV <3)
{
	die "usage:  $0.pl  gff fasta snp \n";
}


my $gff=Bio::FeatureIO->new(-file => "$ARGV[0]", -format => 'GFF', -version => 3);
my $fa = Bio::SeqIO->new(-file => "$ARGV[1]", -format => 'Fasta');
open my $snp_file, $ARGV[2] or die;
#open my $out_file, ">$ARGV[2].stat.via_ancestor" or die;
#open $snp, "chromosome01.genotype.snp" or die;

warn "Reading fasta...";
#open $snp_file, "tmp.snp" or die;
warn "Done\n";



#######################read fasta##########################
my %seq;
while(my $tmp = $fa->next_seq())
{
	        $seq{$tmp->display_id()} = $tmp;
}


warn "Reading gff...\n";



######################read gff into memory##################


##THESE TWO VARs IS ALSO USED IN SUBROUTINE AS GLOBAL VAR
#like a index
##THESE TWO VARs IS ALSO USED IN SUBROUTINE AS GLOBAL VAR
#like a index
our %array_of_gene_start;#used for determin the loci's nearest gene, sorted by ascending;
our %start2end;#return end of the gene given its start and chr
our %start2gene;#return the gene of this start

my %gene; #containing total genes from this gff file
while ( my $feature = $gff->next_feature() ) {


#       no strict 'refs';
#       print "Instance METHOD IS  " . Dumper( \%{ref ($seq)."::" }) ;exit;
        if($feature->type ->name eq "mRNA")
        {
                $start2gene{$feature->seq_id}{$feature -> start} = $feature -> get_Annotations ('Name')->value;
                $start2end{$feature->seq_id}{$feature -> start}  = $feature -> end;
        }

        next unless($feature-> type ->name eq "CDS");

        push  @{$gene{($feature -> get_Annotations ('Parent'))->value}}, $feature;

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

while(<$snp_file>)
{
#preparint SNP
#my $start = Time::HiRes::time();
	chomp;
	next if($_ eq "");
	my ($chr, $location, %population, $nipponbare_base);
	($chr, $location, $nipponbare_base, $population{indica}{input}, $population{japonica}{input}, $population{nivara}{input}, $population{rufipogon}{input}) = split /\t/, $_;
#warn (Time::HiRes::time()-$start); $start = Time::HiRes::time();


	for my $species (keys %population)
	{
		my @SNP = split /\s/, $population{$species}{input};
		for my $tmp(@SNP)
		{
			$population{$species}{base_count}{substr( $tmp, 0, 1)}=substr( $tmp, 1);
			$population{$species}{total_count}+=substr( $tmp, 1);
		}
	}

#call concensus genotype
        my (%combined_count, $con_genotype, $count);
        $count = 0;
#combine the two populations' data
	for my $species (keys %population)
	{
                        for my $base(qw/A T C G/) #sort keys %{$population{$sp_ref}{base_count}})
                        {
                                if(exists $population{$species}{base_count}{$base})
                                {
                                        $combined_count{$base} += $population{$species}{base_count}{$base};
                                }
#find the concensus genotype
			}
	}
	
	for my $base(qw/A T C G/)
	{
                                if(exists $combined_count{$base})
                                 {
                                         if($combined_count{$base} >$count)#found one max genotype
                                        {
                                                $con_genotype = $base;
                                                $count = $combined_count{$base};
                                        }
                                 }
	}

#warn (Time::HiRes::time()-$start);
#$start = Time::HiRes::time();
	for my $species (keys %population)
	{

##########################
# count SNPs by ancestor_snp's genotype###
###########################
		for my $wild_base (sort keys %{$population{$species}{base_count}})
		{
			last unless ( defined $ancestor_snp[$location]);
#warn (Time::HiRes::time()-$start); $start = Time::HiRes::time();


			my @tmp_array = &SNP_stat($chr, $location, 
				$con_genotype, 
				$wild_base, 
				\%chr_mapper, \%cds_mapper, \%seq);

			for my $tmp(@tmp_array)
			{
			next if ($$tmp{strand} eq "NA");
#warn (Time::HiRes::time()-$start); $start = Time::HiRes::time();
			$frequency{$species}{$$tmp{mutation_type}} {
				int(
				1000*#0.001 scale
				$population{$species}{base_count}{$wild_base} /
				$population{$species}{total_count}
				)/1000
			}++;

			}
#warn (Time::HiRes::time()-$start); $start = Time::HiRes::time();exit;
		}
	}
#warn (Time::HiRes::time()-$start); $start = Time::HiRes::time();exit;
#warn (Time::HiRes::time()-$start); $start = Time::HiRes::time();
#exit;
#	print $tmp{$_}, "\t" for (qw/chromosome strand location_chr ref_base snp_base location_cds phase ref_codon snp_codon mutation_type/);
#	print "\n";

}#end reading snp

warn "Done\nProgram finished..Outputing result into $ARGV[2].stat\n";

for my $species (keys %frequency)
{

	open my $out_file, ">$ARGV[2].stat.via_concensus.$species" or die;

		for my $key (sort keys %{$frequency{$species}})
		{
			print $out_file $key, "\n";
			print $out_file $_,"\t",$frequency{$species}{$key}{$_},"\n" for (sort keys %{$frequency{$species}{$key}});
		}
}

warn "Clearing memory...\n";

#############sub routine###############


#accept: chromosome_id  coordinate_on_chr Ref-base SNP-base \%chr_mapper, \%cds_mapper
#return an a hash containing the following elements
#chromosome strand location_chr ref_base snp_base location_cds phase ref_codon snp_codon mutation_type

sub SNP_stat
{
#my $start = Time::HiRes::time();

#print @DescriptionForCodonStat;
#print %CodonStat;
#exit;
        #all results will be stored in here
        #chromosome strand location_chr ref_base snp_base location_cds phase ref_codon snp_codon mutation_type
        my @results;
        my %hash_return;
        $hash_return{$_} = "NA" for(qw/chromosome strand location_chr ref_base snp_base location_cds phase ref_codon snp_codon mutation_type ref_pep snp_pep gene_name/);


        my ($chr_mapper, $cds_mapper, $seq);
        ($hash_return{"chromosome"}, $hash_return{"location_chr"}, $hash_return{"ref_base"}, $hash_return{"snp_base"}, $chr_mapper, $cds_mapper, $seq) = @_;

        #like: chromosome01 100 A G \%chr_mapper \%_cds_mapper

#warn (Time::HiRes::time()-$start);
#$start = Time::HiRes::time();
        return @results if($hash_return{"ref_base"} eq $hash_return{"snp_base"} or !exists $hash_return{"ref_base"} or !exists $hash_return{"snp_base"});


##Sear for nearest gene##
        if($hash_return{"location_chr"} <$array_of_gene_start{$hash_return{"chromosome"}}[0])
        {
                return @results;# this SNP gets up too early, no gene is evolved yet!
        }
#from the first to the last, no one escapes!
        my (@nearest_gene);#, $nearest_gene_start);
        for my $search_for_nearest_gene( @{$array_of_gene_start{$hash_return{"chromosome"}}})
        {
                if($hash_return{"location_chr"} >= $search_for_nearest_gene and
                        $hash_return{"location_chr"}<= $start2end{$hash_return{"chromosome"}}{$search_for_nearest_gene})#found it's stage, next is future , now is most important
                {
#                       last;
                        push @nearest_gene, $start2gene{$hash_return{"chromosome"}}{$search_for_nearest_gene};
                }
        #       else
        #       {
        #               $nearest_gene_start = $search_for_nearest_gene;
        #       }
        }
#TODO 改为正链找一次，然后负链找一次


        my $loci_chr = Bio::Location::Simple->new(
                -start => $hash_return{"location_chr"},
                -end => $hash_return{"location_chr"},
                -seq_id => $hash_return{"chromosome"}
                );
        for my $this_gene(@nearest_gene)
        {
#clear cache
                $hash_return{$_} = "NA" for(qw/strand location_cds phase ref_codon snp_codon mutation_type ref_pep snp_pep gene_name/);
                $hash_return{gene_name} = $this_gene;
                next unless (exists $$chr_mapper{$this_gene});

                my $loci_cds = $$chr_mapper{$this_gene}->map($loci_chr);

#found the corresponding gene
                if(defined $loci_cds ->start)
                {
                        $hash_return{"location_cds"} = $loci_cds ->start;

                        $hash_return{"strand"} = $$chr_mapper{$this_gene} ->map($loci_chr) ->strand;

                        my @codon_phase_delay;
                        my (@ref_codon, @snp_codon);
                        my (@codon_loci, @codon_base);

                        $hash_return{"phase"}=$hash_return{"location_cds"}%3 -1;#get phase and transform from 1 2 0 to 0 1 2 ;
                        $hash_return{"phase"} +=3 if($hash_return{"phase"} <0);#0 1 2

#                       if($loci_start->strand >0) seems free of strand problem
#                               {
                                 if($hash_return{"phase"} == 0)
                                {
                                        @codon_phase_delay=(0, 1, 2);
                                }
                                elsif($hash_return{"phase"} == 1)
                                {
                                        @codon_phase_delay = (-1, 0, 1);
                                }
                                elsif($hash_return{"phase"} == 2)
                                {
                                        @codon_phase_delay = (-2, -1, 0);
                                }

 #       eval
#        {
                                @codon_loci = map {
                                        $$cds_mapper{$this_gene}->map(
                                        Bio::Location::Simple->new(
                                        -start=>$hash_return{"location_cds"}+$_,
                                        -end=>$hash_return{"location_cds"}+$_,
                                        )
                                        ) ->start
                                        }@codon_phase_delay;

                                        my $tmp_ref = $hash_return{"ref_base"};
                                        my $tmp_snp = $hash_return{"snp_base"};

                                if($hash_return{"strand"} > 0)
                                {
                                        @codon_base = map{uc($$seq{$hash_return{"chromosome"}}->subseq($_,$_))}@codon_loci;
                                }
                                else
                                {
                                        @codon_base = map{ uc($$seq{$hash_return{"chromosome"}}->subseq($_,$_))}@codon_loci;
                                        $_ =~ tr/ATCG/TAGC/ for(@codon_base);

                                        $tmp_ref =~ tr/ATCG/TAGC/;
                                        $tmp_snp =~ tr/ATCG/TAGC/;
                                }

                                @ref_codon = @codon_base;#TODO use Seq object instead
                                $ref_codon[$hash_return{"phase"}] = $tmp_ref;
                                @snp_codon = @codon_base;
                                $snp_codon[$hash_return{"phase"}] = $tmp_snp;
  #      };
#                       }

                        $hash_return{"ref_codon"} = join("", @ref_codon);
                        $hash_return{"snp_codon"} = join("", @snp_codon);

                        $hash_return{"mutation_type"} = $DescriptionForCodonStat[$CodonStat{$hash_return{"ref_codon"}}{$hash_return{"snp_codon"}}];

                        use Bio::Tools::CodonTable;
                        my $myCodonTable   = Bio::Tools::CodonTable->new();

                        $hash_return{"ref_pep"} = $myCodonTable->translate($hash_return{"ref_codon"});
                        $hash_return{"snp_pep"} = $myCodonTable->translate($hash_return{"snp_codon"});

#               return %hash_return; #ignoring overlapping genes
                        push @results, {%hash_return};
                }#endif
        }#end for
#warn (Time::HiRes::time()-$start); $start = Time::HiRes::time();exit;
#       return %hash_return;#no corresponding geen found
        return @results;
}




