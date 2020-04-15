sub Find_singlecopy {
    my $category_file = shift;
    my $orthomcl_out = shift;
#   print "$orthomcl_out\n";
#   my $cds_id = shift;

    open IN,"$category_file" || die "$!";
    my @specs;
    while (<IN>){
        chomp;
        my ($spec, undef) = split /\s+/,$_;
        push @specs,$spec;
    }
    close IN;
    my $muscle_shell = "$Outdir/muscle_shell.sh";
    my $pep2cds = "$Outdir/pep2cds.sh";

    open MUS,">$muscle_shell" || die "$!";
#   open PEP,">$pep2cds" || die "$!";

    open MCL,"$orthomcl_out" || die "$!";
    $/ = "\n";
    my @gene_seq;
    while (<MCL>){
        my ($id,$cluster) = split /\:\s+/,$_;
        my @clus_id = split /\s+/,$cluster;
        my @gene_num = ();
        my @gene_id = ();
        for (my $i=0;$i<@specs;$i++){
            my $num = 0;
            my $id = "";
            for (my $j=0;$j<@clus_id;$j++){
                if ($clus_id[$j] =~ /$specs[$i]/){
                    $num ++;
                    $clus_id[$j] =~ s/\($specs[$i]\)$//;
                    $id .= $clus_id[$j];
                }
            }
            push @gene_num,$num;
            push @gene_id,$id;
        }
        my $flag = "y";
        for (my $i=0;$i<@specs;$i++){
            if ($gene_num[$i] != 1){
                $flag = "n";
                next;
            }
        }
        if ($flag eq "y") {
            my $cds_dir = &Make_dir_cds(@gene_id); ## put the single_copy into family
            my $pep_file = &Cds_aa($cds_dir);
            my $pep_muscle_file = $pep_file;
            $pep_muscle_file .= '.muscle';
            print MUS "$muscle -in $pep_file -out $pep_muscle_file; perl $Bin/pepMfa_to_cdsMfa.pl $pep_muscle_file $cds_dir > $cds_dir.muscle;\n";
#           `muscle -in $pep_file -out $pep_muscle_file -quiet`;
#           my $cds_muscle_file = &Aa_cds($pep_muscle_file,$cds_dir);
#           my $kaks_file = &Calculate_kaks($philip_file);
        }
    }
    close MCL;
    close MUS;
#   close PEP;
    ($muscle_shell);
}

