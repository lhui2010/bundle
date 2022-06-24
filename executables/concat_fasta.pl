#!/usr/bin/env perl 
#https://github.com/santiagosnchez/concat_fasta/blob/master/concat_fasta.pl
# ©Santiago Sanchez-Ramirez, University of Toronto

if (($ARGV[0] =~ m/^-+he{0,1}l{0,1}p{0,1}$/) or (length($ARGV[0]) == 0)){
    die "
Try:
perl concat_fasta.pl [ --suffix | --prefix ] your_pattern --outfile output [ --filter filter_pattern ] [ -n | -p ]
     --suffix, --prefix         A string of characters present only in the FASTA files you wish to concatenate.
                                An example would be \".fasta\" or \".fas\" (without quotation marks) in the case of the suffix,
                                or \"xxx\" or \"sample\" at the begining of the file in the case of the prefix.
     --outfile                  The name of file where you are saving the concatenated data.
                                By default the format is FASTA, unless the following flags are stated.
     --filter    [Optional]     Will include only the terminals that include the pattern.
                                This can also be a text file with a list of filter patterns.
     --no-concat                Use only if you just want to filter records while keeping individual gene alignments.
     -n          [Optional]     The output format is NEXUS.
                                The partitions will be printed at the end of the file in NEXUS style.
     -p          [Optional]     The output format is PHYLIP.
                                Partitions are printed to the screen or can be printed to a file by using \">\"
Example:
perl concat_fasta.pl --suffix .fasta --outfile concat.phylip --filter sp3 -p > part.txt\n\n";
}

my $suffix;
my $prefix;
my $outfile;
my $filter;
my $nexfl=0;
my $phyfl=0;
my @list=();
my $nc=0;

for my $i (0..$#ARGV){
    if ($ARGV[$i] =~ m/^--suf/){
        $suffix = $ARGV[$i+1];
    }
    if ($ARGV[$i] =~ m/^--pref/){
        $prefix = $ARGV[$i+1];
    }
    if ($ARGV[$i] =~ m/^--out/){
        $outfile = $ARGV[$i+1];
    }
    if ($ARGV[$i] =~ m/^--no-conc/){
        $nc = 1;
    }
    if ($ARGV[$i] =~ m/^--fil/){
        $filter = $ARGV[$i+1];
        if (open(FIL, "<", $filter)){
            while(<FIL>) { chomp($_); push @list, $_ }
        } else {
            print "Filter $filter is not a file.\n";
        }
    }
    if ($ARGV[$i] =~ m/^-n/){
        $nexfl = 1;
    }
    if ($ARGV[$i] =~ m/^-p/){
        $phyfl = 1;
    }
}

my @blocks = ();
my @master = ();
my @lengths = ();
my @allLab = ();
my $nchar;

my @files = ();

if ((length($prefix) == 0) and (length($suffix) > 0)){
    @files = `ls *$suffix`;
}
elsif ((length($suffix) == 0) and (length($prefix) > 0)){
    @files = `ls $prefix*`;
}
elsif ((length($suffix) > 0) and (length($prefix) > 0)){
    @files = `ls $prefix*$suffix`;
} else {
    die "Suffix/Prefix not specified. Try: concat_fasta.pl -h\n";
}

if (scalar(@files) == 0){ die "Suffix/Prefix not found in fasta file names\n"; }

foreach(@files){
    my $data;
    open FASTA, "<$_";
    while(<FASTA>){
        next if (/^$/);
        $data .= $_;
    }
    close FASTA;
    push @blocks, $data;
}

for (my $i=0; $i<scalar(@blocks); ++$i){
    my @tempArr = split "\n", $blocks[$i];
    my @tempLab=();
    my $allSeqs;
    my %tempHash=();
    for (my $j=0; $j<scalar(@tempArr); ++$j){
        if ($tempArr[$j] =~ m/^>/){
            push @tempLab, $tempArr[$j];
            push @allLab, $tempArr[$j];
        } else {
            $allSeqs .= $tempArr[$j];
        }
    }
    my $genlen = length($allSeqs)/scalar(@tempLab);
    my @tempSeq = $allSeqs =~ m/(.{$genlen})/g;
    if (scalar(@tempSeq) == scalar(@tempLab)){
        push @lengths, $genlen;
        for my $j (0..$#tempSeq){
            $tempHash{$tempLab[$j]} = $tempSeq[$j];
        }
    } else {
        die "Format error...\n".scalar(@tempSeq)."\n".scalar(@tempLab);
    }
    push @master, {%tempHash};
}

$nchar = eval join('+', @lengths);

if ((length($prefix) == 0) and (length($suffix) > 0)){
    @files = `ls *$suffix`;
}
elsif ((length($suffix) == 0) and (length($prefix) > 0)){
    @files = `ls $prefix*`;
}
elsif ((length($suffix) > 0) and (length($prefix) > 0)){
    @files = `ls $prefix*$suffix`;
} else {
    die "Suffix/Prefix not specified. Try: concat_fasta.pl -h\n";
}

if (($nc != 0) and (length($filter) != 0)){
    if (scalar(@list) == 0){
        `mkdir ../filtered_$filter`;
        print "Filtering taxa without concatenation.\nFiles will be stored in ../filtered_$filter\n";
        for my $i (0 .. $#master){
            chomp($files[$i]);
            open(OUT, ">", "../filtered_$filter/$files[$i]");
            %tempH = %{ $master[$i] };
            @k = keys %tempH;
            my @kf=();
            @kf = grep { /$filter/ } @k;
            foreach(@kf){
                print OUT $_ . "\n" . $tempH{$_} . "\n";
            }
            close OUT;
        }
    } else {
        for $list (@list){
            `mkdir ../filtered_$list`;
            print "Filtering taxa without concatenation.\nFiles will be stored in ../filtered_$list\n";
            for my $i (0 .. $#master){
                chomp($files[$i]);
                open(OUT, ">", "../filtered_$list/$files[$i]");
                %tempH = %{ $master[$i] };
                @k = keys %tempH;
                my @kf=();
                @kf = grep { /$list/ } @k;
                foreach(@kf){
                    print OUT $_ . "\n" . $tempH{$_} . "\n";
                }
                close OUT;
            }
        }
    }
} else {

my @singLab = sort {$a cmp $b} uniq(@allLab);
my $ntax = scalar(@singLab);

my @fil=();

if (length($filter) > 0){
    if (scalar(@list) > 0){
        for $pat (@list){
            @temp = grep { /$pat/ } @singLab;
            push @fil, @temp;
        }
    } else {
        @fil = grep { /$filter/ } @singLab;
    }
    $ntax = scalar(@fil);
}

if ((length($prefix) == 0) and (length($suffix) > 0)){
    @files = `ls *$suffix`;
}
elsif ((length($suffix) == 0) and (length($prefix) > 0)){
    @files = `ls $prefix*`;
}
elsif ((length($suffix) > 0) and (length($prefix) > 0)){
    @files = `ls $prefix*$suffix`;
} else {
    die "Suffix/Prefix not specified. Try: concat_fasta.pl -h\n";
}

#print STDERR $outfile;
open OUTFILE, ">$outfile" or die;

my $labellength = max(@singLab);

if ($nexfl==1) {
    print OUTFILE "#NEXUS\n\nBegin DATA;\n\tDimensions ntax=$ntax nchar=$nchar;\n\tFormat datatype=DNA gap=- missing=?;\n\tMatrix\n";
    if (scalar(@fil) > 0){
        for my $i (0..$#fil){
            for my $j (0..(scalar(@master)-1)){
                if ($j == 0){
                    print OUTFILE "\t" . substr($fil[$i],1) . " " x ($labellength-length($fil[$i])+1);
                }
                if (${@master[$j]}{$fil[$i]}){
                    print OUTFILE ${@master[$j]}{$fil[$i]};
                } else {
                    print OUTFILE '?' x $lengths[$j];
                }
            }
            print OUTFILE "\n";
        }
        print OUTFILE ";\nEnd;\n\n";
    } else {
        for my $i (0..$#singLab){
            for my $j (0..(scalar(@master)-1)){
                if ($j == 0){
                    print OUTFILE "\t" . substr($singLab[$i],1) . " " x ($labellength-length($singLab[$i])+1);
                }
                if (${@master[$j]}{$singLab[$i]}){
                    print OUTFILE ${@master[$j]}{$singLab[$i]};
                } else {
                    print OUTFILE '?' x $lengths[$j];
                }
            }
            print OUTFILE "\n";
        }
        print OUTFILE ";\nEnd;\n\n";
    }
} elsif ($phyfl==1) {
    print OUTFILE "\t$ntax\t$nchar\n";
    if (scalar(@fil) > 0){
        for my $i (0..$#fil){
            for my $j (0..(scalar(@master)-1)){
                if ($j == 0){
                    print OUTFILE "\t" . substr($fil[$i],1) . " " x ($labellength-length($fil[$i])+1);
                }
                if (${@master[$j]}{$fil[$i]}){
                    print OUTFILE ${@master[$j]}{$fil[$i]};
                } else {
                    print OUTFILE '?' x $lengths[$j];
                }
            }
            print OUTFILE "\n";
        }
        close OUTFILE;
    } else {
        for my $i (0..$#singLab){
            for my $j (0..(scalar(@master)-1)){
                if ($j == 0){
                    print OUTFILE "\t" . substr($singLab[$i],1) . " " x ($labellength-length($singLab[$i])+1);
                }
                if (${@master[$j]}{$singLab[$i]}){
                    print OUTFILE ${@master[$j]}{$singLab[$i]};
                } else {
                    print OUTFILE '?' x $lengths[$j];
                }
            }
            print OUTFILE "\n";
        }
        close OUTFILE;
    }
} else {
    if (scalar(@fil) > 0){
        for my $i (0..$#fil){
            for my $j (0..(scalar(@master)-1)){
                if ($j == 0){
                    print OUTFILE $fil[$i] . "\n";
                }
                if (${@master[$j]}{$fil[$i]}){
                    print OUTFILE ${@master[$j]}{$fil[$i]};
                } else {
                    print OUTFILE '?' x $lengths[$j];
                }
            }
            print OUTFILE "\n";
        }
        close OUTFILE;
    } else {
        for my $i (0..$#singLab){
            for my $j (0..(scalar(@master)-1)){
                if ($j == 0){
                    print OUTFILE $singLab[$i] . "\n";
                }
                if (${@master[$j]}{$singLab[$i]}){
                    print OUTFILE ${@master[$j]}{$singLab[$i]};
                } else {
                    print OUTFILE '?' x $lengths[$j];
                }
            }
            print OUTFILE "\n";
        }
        close OUTFILE;
    }
}

print "Partitions:\n";

for (my $i=0; $i<scalar(@lengths); ++$i){
    $files[$i] =~ s/\.[fF][a-zA-Z]*$//;
    $files[$i] =~ s/[\.-]/_/g;
    $x += $lengths[$i];
    chomp($files[$i]);
    if ($phyfl==1){
        print "DNA, ";
    }
    print $files[$i] . " = ";
    print $x-($lengths[$i]-1) . "-" . $x . "\n";
}

my $x;
if ($nexfl==1){
    print OUTFILE "Begin MrBayes;\n";
    for (my $i=0; $i<scalar(@lengths); ++$i){
        $x += $lengths[$i];
        chomp($files[$i]);
        print OUTFILE "\tcharset " . $files[$i] . " = ";
        print OUTFILE $x-($lengths[$i]-1) . "-" . $x . ";" . "\n";
    }
    $ngenes = scalar(@files);
    print OUTFILE "\n\tpartition all = $ngenes:";
    for (my $i=0; $i<scalar(@files); ++$i){
        chomp($files[$i]);
        $files[$i] =~ s/\.(fas|fasta)//;
        if ($i != $#files){
            print OUTFILE $files[$i] . ",";
        } else {
            print OUTFILE "$files[$i];\n";
        }    
    }
    print OUTFILE "\n\t[Change models settings!]\n\n";
    print OUTFILE "\tlset applyto=(all) rates=invgamma nst=6;\n\tprset applyto=(all) ratepr=variable;\n\tunlink shape=(all) pinvar=(all) revmat=(all);\n\n";
    print OUTFILE "\tmcmc ngen=5000000 printfreq=1000 samplefreq=1000 diagnfreq=10000;\nEnd;\n\n";
    close OUTFILE;        
}

}
# Subroutines

sub max {
    my @res=();
    foreach(@_){
        push @res, length($_);
    }
    my @res2 = sort {$b <=> $a} @res; 
    return(@res2[0]);
}

sub uniq {
    my %seen;
    return grep !$seen{$_}++, @_;
}

sub filtering {
    my @fil=();
    if (scalar(@list) > 0){
        for $pat (@list){
            push @fil, grep { /$pat/ } @_;
        }
    } else {
        push @fil, grep { /$filter/ } @_;
    }
    return(@fil);
}

