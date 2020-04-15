#!/usr/bin/perl -w
# Parsing BLAST reports with BioPerl's Bio::SearchIO module
# WI Biocomputing course - Bioinformatics for Biologists - October 2003

# See help at http://www.bioperl.org/HOWTOs/html/SearchIO.html for all data that can be extracted

use Bio::SearchIO;

# Prompt the user for the file name if it's not an argument
# NOTE: BLAST file must be in text (not html) format


$extend = 0;

#$faname = $ARGV[0];

$inFile = $ARGV[0];

$report = new Bio::SearchIO(
         -file=>"$inFile",
              -format => "blast"); 

#print "QueryAcc\tHitDesc\tHitSignif\tHSP_rank\t\%ID\teValue\tHSP_length\n";

my (%homolog_gene, %seq);
# Go through BLAST reports one by one              
while($result = $report->next_result) 
{
   # Go through each each matching sequence
   while($hit = $result->next_hit) 
   { 
      # Go through each each HSP for this sequence
        while ($hsp = $hit->next_hsp)
         { 
            # Print some tab-delimited data about this HSP
		$start = $hsp->start('hit');
		$end = $hsp->end('hit');
		($start, $end) = ($end, $start) if ($start >$end);
	  	$homolog_gene{$hit->name}{$result->query_name}{'start'} = $start if(!exists $homolog_gene{$hit->name}{$result->query_name}{'start'} or $homolog_gene{$hit->name}{$result->query_name}{'start'} >$start);
	    	$homolog_gene{$hit->name}{$result->query_name}{'end'} = $end if (!exists $homolog_gene{$hit->name}{$result->query_name}{'end'} or $homolog_gene{$hit->name}{$result->query_name}{'end'} < $end);
	    
            
            print $result->query_name, "\t";
	    print $hsp->start('query'), "\t";
	    print $hsp->end('query'), "\t";
#            print $hit->description, "\t";
	    print $hit->name, "\t";
	    print $hsp->start('hit'), "\t";
	    print $hsp->end('hit'), "\t";
 #           print $hit->significance, "\t";
 #          print $hsp->rank, "\t";
            print $hsp->percent_identity, "\t";
	    print $hsp->evalue, "\t";
	    print $hsp->end('query') - $hsp->start('query'), "\t";
	    print $hsp->end('hit') - $hsp->start('hit'), "\n";
  #          print $hsp->hsp_length, "\n";
      } 
   } 
}
exit;
open FA, $faname or die;
$/ = '>';
$_=<FA>;
while(<FA>)
{
        chomp;
       ($name, $seq) = split /\n/, $_, 2;
        $seq =~ s/\n//g;
        $seq{$name} = $seq;

        next unless (defined $seq and $seq ne "");
}


for $scaffold(sort keys %homolog_gene)
{
	'for $rgene( keys %{$homolog_gene{$scaffold}})
	{
	#	print $rgene, "\t";
	#	print "$scaffold\t";
	#	print $homolog_gene{$scaffold}{$rgene}{start}, "\t";
	#	print $homolog_gene{$scaffold}{$rgene}{end}, "\t";
	#	print "\n";
		for $another_rgene( keys %{$homolog_gene{$scaffold}})
		{
			next if ($rgene eq $another_rgene);
			next unless (exists $homolog_gene{$scaffold}{$rgene} and exists $homolog_gene{$scaffold}{$another_rgene});
#sense strand
			if ($homolog_gene{$scaffold}{$rgene}{start} - $homolog_gene{$scaffold}{$another_rgene}{end} <5000 and $homolog_gene{$scaffold}{$rgene}{start} - $homolog_gene{$scaffold}{$another_rgene}{end} >0)
			{
				$combined_name = $another_rgene."+".$rgene;
				$combined_start = $homolog_gene{$scaffold}{$another_rgene}{start};
				$combined_end = $homolog_gene{$scaffold}{$rgene}{end};
				delete $homolog_gene{$scaffold}{$rgene};
				delete $homolog_gene{$scaffold}{$another_rgene};
				$homolog_gene{$scaffold}{$combined_name}{start} = $combined_start;
				$homolog_gene{$scaffold}{$combined_name}{end} = $combined_end;
			}
		}
	}
        for $rgene( keys %{$homolog_gene{$scaffold}})
        {
        #       print $rgene, "\t";
        #       print "$scaffold\t";
        #       print $homolog_gene{$scaffold}{$rgene}{start}, "\t";
        #       print $homolog_gene{$scaffold}{$rgene}{end}, "\t";
        #       print "\n";
                for $another_rgene( keys %{$homolog_gene{$scaffold}})
                {
                        next if ($rgene eq $another_rgene);
                        next unless (exists $homolog_gene{$scaffold}{$rgene} and exists $homolog_gene{$scaffold}{$another_rgene});
#sense strand
                        if ($homolog_gene{$scaffold}{$rgene}{start} - $homolog_gene{$scaffold}{$another_rgene}{end} <5000)
                        {
                                $combined_name = $another_rgene."+".$rgene;
                                $combined_start = $homolog_gene{$scaffold}{$another_rgene}{start};
                                $combined_end = $homolog_gene{$scaffold}{$rgene}{end};
                                delete $homolog_gene{$scaffold}{$rgene};
                                delete $homolog_gene{$scaffold}{$another_rgene};
                                $homolog_gene{$scaffold}{$combined_name}{start} = $combined_start;
                                $homolog_gene{$scaffold}{$combined_name}{end} = $combined_end;
                        }
                }
        }

        for $rgene( keys %{$homolog_gene{$scaffold}})
        {
        #       print $rgene, "\t";
        #       print "$scaffold\t";
        #       print $homolog_gene{$scaffold}{$rgene}{start}, "\t";
        #       print $homolog_gene{$scaffold}{$rgene}{end}, "\t";
        #       print "\n";
                for $another_rgene( keys %{$homolog_gene{$scaffold}})
                {
                        next if ($rgene eq $another_rgene);
                        next unless (exists $homolog_gene{$scaffold}{$rgene} and exists $homolog_gene{$scaffold}{$another_rgene});
#sense strand
                        if ($homolog_gene{$scaffold}{$rgene}{start} - $homolog_gene{$scaffold}{$another_rgene}{end} <5000)
                        {
                                $combined_name = $another_rgene."+".$rgene;
                                $combined_start = $homolog_gene{$scaffold}{$another_rgene}{start};
                                $combined_end = $homolog_gene{$scaffold}{$rgene}{end};
                                delete $homolog_gene{$scaffold}{$rgene};
                                delete $homolog_gene{$scaffold}{$another_rgene};
                                $homolog_gene{$scaffold}{$combined_name}{start} = $combined_start;
                                $homolog_gene{$scaffold}{$combined_name}{end} = $combined_end;
                        }
                }
        }exit;
';
	for $rgene( sort keys %{$homolog_gene{$scaffold}})
	{
	#	$rgene, "\t";
	#	print "$scaffold\t";
		$start =  $homolog_gene{$scaffold}{$rgene}{start} -$extend;
		$end =  $homolog_gene{$scaffold}{$rgene}{end}+$extend;
		$start = 0 if ($start <0);
		$len = $end - $start +1;
		print ">$scaffold|$start-$end|Homolog_of_$rgene\n";
		&formated_print (substr($seq{$scaffold}, $start, $len));
		#print "\n";
	}
}

sub formated_print
{
	        my $length = 70;
		my $seq_p =shift;
                for ( my $pos = 0 ; $pos < length($seq_p) ; $pos += $length ) {
                     print substr($seq_p, $pos, $length), "\n";
                }
}
