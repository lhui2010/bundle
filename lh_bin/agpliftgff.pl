#!/usr/bin/perl

=head1 NAME  
  
  agpliftgff.pl

=head1 DESCRIPTION

  convert GFF locations from old to new assembly that differs
  in minor ways identified by oldAgp, newAgp scaffold assembly tables.
  after Jim Kent's agpLift
  
  *** NOTE: This version NOW ONLY handles case of changing N spacer fragments.

  my $ok=&GetOptions( 
  "oldagp=s" => \$oldagp,
  "newagp=s" => \$newagp,
  "dirnew=s" => \$dirnew,
  "gff=s" => \@gff,
  "test!" => \$test,
  "verbose!" => \$verbose,
    );
  
=cut  

use strict;
use Getopt::Long;
use Digest::MD5;

my $test = 0;
my $verbose = 1;
my $refmap= 1;
my $domd5 = 1;
my $oldagp= undef;
my $newagp= undef;
my $dirnew= undef;
my @gff= ();

warn "This ONLY handles case of changing N spacer fragments\n";

my $ok=&GetOptions( 
  "oldagp=s" => \$oldagp,
  "newagp=s" => \$newagp,
  "dirnew=s" => \$dirnew,
  "gff=s" => \@gff,
  "test!" => \$test,
  ## "refmap!" => \$refmap,
  "verbose!" => \$verbose,
  );


push(@gff, @ARGV);
die "See perldoc $0 for help \n"
  unless($ok && $oldagp && $newagp && ($test || @gff));

$dirnew =~ s,/$,, if($dirnew);
unless(-d $dirnew) { warn "mkdir $dirnew\n"; mkdir($dirnew); }

my $agpold= readAgp($oldagp);
my $agpnew= readAgp($newagp);

foreach my $gff (@gff) {
  (my $gffnew= $gff) =~ s/\.gz//; $gffnew.=".lft";
  $gffnew =~ s,^.*/,$dirnew/, if($dirnew);
  liftByAgp($agpold,$agpnew,$gff,$gffnew);
}
exit;

#--------------------------------------------------------------------------

sub readAgp {
  my($agpfile)= @_;
  my %agp; 
  my($ngap,$ncontig);
  open(F,$agpfile) or die "cant read $agpfile";
  while(<F>){
    next unless(/^\w/);
    chomp;
    my @v= split"\t";
    # WN_TYPE == W,N
    my ($scafid, $scaf_b, $scaf_e, $agpi, $WN_TYPE, 
        $ctgid, $ctgstart, $ctglen, $gaps, $FRAGMENT, $YES, $PLUS);
    if($v[4] eq 'N' || $v[6] eq 'fragment') {
      ($scafid, $scaf_b, $scaf_e, $agpi, $WN_TYPE, $gaps, $FRAGMENT, $YES)= @v;
      $agpi--; # 0-origin
      my $agpr= [ $scaf_b, $scaf_e, $WN_TYPE, $gaps, $agpi, 0, $scafid, ];
      $agp{$scafid}[$agpi]= $agpr;
      $ngap++;
      # push( @{$scafs{$scafid}}, $agpr);
    } elsif($v[4] eq 'W') {
      ($scafid, $scaf_b, $scaf_e, $agpi, $WN_TYPE, $ctgid, $ctgstart, $ctglen, $PLUS)= @v;
      $agpi--; # 0-origin
      my $agpr= [  $scaf_b, $scaf_e, $WN_TYPE, $ctglen, $agpi, $ctgid, $scafid, ];
      $agp{$scafid}[$agpi]= $agpr;
      $agp{$ctgid}= $agpr; # for crossref by ctgid
      $ncontig++;
    } else {
      warn "Error unknown line: $_\n"; next;
    }
  }
  warn "# read $agpfile: n_scaf=",scalar(keys(%agp))," n_contig=",$ncontig," ngap=",$ngap,"\n" if ($verbose);
  return \%agp;
}

##  FIXME: case where new has dropped (added?) N/fragment lines
##  need to match ctgid's
# microbe% grep -c fragment *.agp
# dana_caf051209.agp:6800
# dana_caf060210.assembly.agp:6783
# dere_caf051209.agp:2497
# dere_caf060210.assembly.agp:2486
# dgri_caf051209.agp:6724
# dgri_caf060210.assembly.agp:6717
# dmoj_caf051209.agp:5042
# dmoj_caf060210.assembly.agp:5033
# dvir_caf051209.agp:4871
# dvir_caf060210.assembly.agp:4852


sub anewIndex {
  my($ar,$aoldr,$anew)= @_;
  ## new has diff index :((((
  my $ari= $$ar[4];
  # my $arnew= $$anewr[$ari]; ## BAD index often
  my $ctgid= $$ar[5];
  unless($ctgid) { 
    my $ar1= $$aoldr[$ari+1]; # next in list always is contig ?
    $ctgid= $$ar1[5];
    }
  return unless($ctgid);
  return $anew->{$ctgid};
}

sub findOldNewAgp {
   my($aold,$anew,$ref,$start,$stop)= @_;
   my($astart,$astop);
   

   ## dgg, jun07: fixme to lift from scaffold.agp to chrom.agp : diff refs **
   ## use contigids to map from scaf to chr id?

   my $aoldr = $aold->{$ref};
   # my $anewr = $anew->{$ref};
   return (undef,undef,$start,$stop) unless(defined $aoldr); #?? && defined $anewr);
   
   ## also check for -dna ends far past old scaff end .. correct here?
   # my $stopOffBy1 = $stop;
   my $aend= $$aoldr[-1];
   $stop  = $$aend[1] if($$aend[1] < $stop);
   $start = $$aend[1] - 1 if($$aend[1] <= $start); # got this bad one also !
   
   foreach my $ar ( @$aoldr ) {
      # $ar == [  $scaf_b, $scaf_e, $WN_TYPE, $ctglen, $agpi, $ctgid, ];
      my($b,$e)= ($$ar[0],$$ar[1]);
      if($start < $e) {
      
        if($start >= $b && $stop <= $e) { ## contained in
          my $avo= [ $b, $b, $e, $$ar[2], $$ar[3], $$ar[6] ];
          my $arnew= anewIndex($ar, $aoldr, $anew);
          my $avn= [ $$arnew[0], $$arnew[0], $$arnew[1], $$arnew[2], $$arnew[3], $$arnew[6] ];
          return($avo,$avn, $start, $stop);
          }
        elsif (!$astart) { $astart= $ar; }
        
        if ($stop <= $e) {
          my $alen= $e - $$astart[0] + 1;
          my $avo= [ $$astart[0], $b, $e, "C", $alen, $$ar[6] ];

          my $arnew0= anewIndex($astart, $aoldr, $anew);
          my $arnew = anewIndex($ar, $aoldr, $anew);
          my $avn= [ $$arnew0[0], $$arnew[0], $$arnew[1], $$arnew[2], $$arnew[3], $$arnew[6] ];
          
          return($avo,$avn, $start, $stop);
        } else {
          # in between; save what?
        }
      }
    }
    
  return (undef,undef,$start,$stop);  # error ??
}


sub offsetBy {
  my($aold,$anew,$start,$stop)= @_;
  # return($start,$stop) unless(defined $aold && defined $anew);
  
  # $start= $$aold[0] if($start < $$aold[0]); # fix bad data ??
  # $stop = $$aold[2] if($stop > $$aold[2]); # fix bad data ??
  
  $start = $start - $$aold[0] + $$anew[0];
    ## tricky: needs last match not first
  $stop = $stop - $$aold[1] + $$anew[1];
  my $refnew= $$anew[5];
  return ($start,$stop,$refnew);
}

sub liftByAgp {
  my($agpold,$agpnew,$gff,$gffnew)= @_;
 
  warn "# liftByAgp $gff\n" if ($verbose);
  my ($nchange,$nread,$ok);
  if($gff =~ /\.gz/ && -f $gff) { $ok= open(F,"gunzip -c $gff|"); } 
  else { $ok= open(F,$gff); }
  die "cant read $gff" unless ($ok);
  open(OUT,">$gffnew") or die "cant write $gffnew";
  my $outh= *OUT;
  while(<F>){
    unless(/^\w/){ print $outh $_; next; }
    my @v= split"\t";
    my($ref,$gffsource,$type,$start,$stop,$eval,$strand,$offs,$attr)= @v;
    
    ## need to correct some bad endpoints past old scaffolds ?* -dmel-dna aligns
    my ($aold,$anew,$start,$stop)= findOldNewAgp($agpold,$agpnew,$ref,$start,$stop);
    my($nref,$nstart,$nstop)= ($ref,$start,$stop);
    
    if(defined $aold && defined $anew) {
      ($nstart,$nstop,$nref)= offsetBy($aold,$anew,$start,$stop);
    } else {
      warn "# cant find $ref:$start,$stop\n" ;
    }
#     my $aold= findAgp($agpold,$ref,$start,$stop);
#     my $anew= findAgp($agpnew,$ref,$start,$stop);
#     warn "# cant find $ref,$start\n" unless(defined $aold && defined $aold);
#     my($nstart,$nstop)= offsetBy($aold,$anew,$start,$stop);

    $nchange++ if($nstart != $start || $nstop != $stop);
    $nread++;
    print $outh join("\t", $nref,$gffsource,$type,$nstart,$nstop,$eval,$strand,$offs,$attr);
    if($verbose && $nread % 1000 == 0){ print STDERR "."; print STDERR "\n" if($nread % 50000 == 0); }
  }
  close(F); close(OUT);
  warn "\n# wrote $gffnew: nchange=",$nchange," nread=",$nread,"\n" if ($verbose);
}


## FIXME; need cumulative offset old -> new
## make it findOldNewAgp; use ref,agpi index into new
# sub findAgp {
#    my($agp,$ref,$start,$stop)= @_;
#    my($astart,$astop);
#    
#    foreach my $ar (@{$agp->{$ref}} ) {
#       # [  $scaf_b, $scaf_e, $WN_TYPE, $ctglen, $ctgid, $ctgstart, $agpi,  ];
#       my($b,$e)= ($$ar[0],$$ar[1]);
#       if($start < $e) {
#         if($start >= $b && $stop <= $e) { ## contained in
#           # return $ar; #?
#           return [ $b, $b, $e, $$ar[2], $$ar[3] ];
#           }
#         elsif (!$astart) { $astart= $ar; }
#         if ($stop <= $e) {
#           # $astop = $ar;
#           my $alen= $e - $$astart[0] + 1;
#           # my $av= [ $$astart[0], $e, "C", $alen ];
#           my $av= [ $$astart[0], $b, $e, "C", $alen ];
#           return $av;
#         } else {
#           # in between; save what?
#         }
#       }
#     }
#     
#   return undef;  # error ??
# }

__END__

# AGP OLD
scaffold_1696	1	2098	1	W	contig_1918	1	2098	+
scaffold_1696	2099	2130	2	N	32	fragment	yes
scaffold_1696	2131	7778	3	W	contig_1919	1	5648	+
scaffold_1696	7779	7803	4	N	25	fragment	yes
scaffold_1696	7804	12220	5	W	contig_1920	1	4417	+

# AGP NEW
scaffold_1696	1	2098	1	W	contig_1918	1	2098	+
scaffold_1696	2099	2123	2	N	25	fragment	yes	    ** lost 7
scaffold_1696	2124	7771	3	W	contig_1919	1	5648	+
scaffold_1696	7772	7839	4	N	68	fragment	yes	    ** gained cum. 36
scaffold_1696	7840	12256	5	W	contig_1920	1	4417	+

*** 535,546 **** OLD
  scaffold_1683 DGIL_SNP        exon    1985    3136    15.865  -       .       Parent=GG_DGIL_SNP_28100197
  scaffold_1696 DGIL_SNP        gene    854     1633    .       -       .       ID=GG_DGIL_SNP_28100198
  scaffold_1696 DGIL_SNP        exon    854     1633    23.557  -       .       Parent=GG_DGIL_SNP_28100198
! scaffold_1696 DGIL_SNP        gene    4543    5259    .       +       .       ID=GG_DGIL_SNP_28100199
! scaffold_1696 DGIL_SNP        exon    4543    5259    -29.981 +       .       Parent=GG_DGIL_SNP_28100199
! scaffold_1696 DGIL_SNP        gene    9287    10111   .       +       .       ID=GG_DGIL_SNP_28100200
! scaffold_1696 DGIL_SNP        exon    9287    9369    -2.787  +       .       Parent=GG_DGIL_SNP_28100200
! scaffold_1696 DGIL_SNP        exon    9440    9725    3.004   +       .       Parent=GG_DGIL_SNP_28100200
! scaffold_1696 DGIL_SNP        exon    10079   10111   10.237  +       .       Parent=GG_DGIL_SNP_28100200
  scaffold_1704 DGIL_SNP        gene    911     1114    .       -       .       ID=GG_DGIL_SNP_28100201
  scaffold_1704 DGIL_SNP        exon    911     1114    2.787   -       .       Parent=GG_DGIL_SNP_28100201
  scaffold_1705 DGIL_SNP        gene    358     807     .       -       .       ID=GG_DGIL_SNP_28100202
--- 535,546 ---- NEW
  scaffold_1683 DGIL_SNP        exon    1985    3136    15.865  -       .       Parent=GG_DGIL_SNP_28100197
  scaffold_1696 DGIL_SNP        gene    854     1633    .       -       .       ID=GG_DGIL_SNP_28100198
  scaffold_1696 DGIL_SNP        exon    854     1633    23.557  -       .       Parent=GG_DGIL_SNP_28100198
! scaffold_1696 DGIL_SNP        gene    4536    5252    .       +       .       ID=GG_DGIL_SNP_28100199
! scaffold_1696 DGIL_SNP        exon    4536    5252    -29.981 +       .       Parent=GG_DGIL_SNP_28100199
   .. lost 7 ^^
! scaffold_1696 DGIL_SNP        gene    9323    10147   .       +       .       ID=GG_DGIL_SNP_28100200
! scaffold_1696 DGIL_SNP        exon    9323    9405    -2.787  +       .       Parent=GG_DGIL_SNP_28100200
   .. gained 36
! scaffold_1696 DGIL_SNP        exon    9476    9761    3.004   +       .       Parent=GG_DGIL_SNP_28100200
! scaffold_1696 DGIL_SNP        exon    10115   10147   10.237  +       .       Parent=GG_DGIL_SNP_28100200
  scaffold_1704 DGIL_SNP        gene    911     1114    .       -       .       ID=GG_DGIL_SNP_28100201
  scaffold_1704 DGIL_SNP        exon    911     1114    2.787   -       .       Parent=GG_DGIL_SNP_28100201
  scaffold_1705 DGIL_SNP        gene    358     807     .       -       .       ID=GG_DGIL_SNP_28100202
