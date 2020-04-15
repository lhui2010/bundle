#!/usr/bin/perl

#modified by Liu Hui to adjust Multiple N sequence
use strict;
use Getopt::Long;

my %opts;
my $version="1.0 Alvin Chen 2011-06-10";
my $cap3para="-a 11 -i 21 -o 16 -s 251 -p 66";

GetOptions(\%opts,"i=s","m=s","t=s","d=i","o=s","help");
if((!defined $opts{i}) || (!defined $opts{o}) || (!defined $opts{t}) ){
    &Usage();
}

my $parInputfile=$opts{i};
my $parMethod=(defined $opts{m}) ? $opts{m} : "phrap";
my $parOutfile=$opts{o};
my $parDirCorrection=(defined $opts{d}) ? $opts{d} : 0;
my $parTempdir=$opts{t};
my $parLogfile="$parOutfile.log";

if($parMethod ne "phrap" && $parMethod ne "cap3"){
	print "Error: the method for assembling sequences must be \"phrap\" or \"cap3\"\n";
	exit(1);
}

my $iserror;
&Filecheck($parInputfile,0);

open (my $fh_infile,"$parInputfile");
open (my $fh_outfile,">$parOutfile");
open (my $fh_logfile,">$parLogfile");

my $seqN;
my $seqs;
my $seqsDir;
my $seqNs;

print $fh_logfile "Seq_Assembly $version\nJob start at:";
print $fh_logfile time();
print $fh_logfile "\n parameters: -i $parInputfile -m $parMethod -d $parDirCorrection -t $parTempdir\n\n";

while(<$fh_infile>){
	if(/^\#/){
		if($seqs ne ""){
			print $fh_logfile "\#$seqN\n$seqNs\t$seqsDir\n";
			&Assembly($seqN,$seqNs,$seqsDir,$seqs);
		}
		my $seqName=$_;
		$seqName=~s/\#//g;
		my @seqNinfo=split(/\t/,$seqName);
		if($parDirCorrection == 1 && $seqNinfo[1] eq ""){
			print "Warning:No direction info for $seqNinfo[0]\n";
		}
		$seqN=$seqNinfo[0];
		$seqN=~s/\n//g;
		$seqsDir=$seqNinfo[1];
		$seqsDir=~s/\n//g;
		$seqs="";
		$seqNs="";
	}
	else{
		if(/^>(\S+)/){
			if($seqNs eq ""){
				$seqNs=$1;
			}
			else{
				$seqNs.=",$1";
			}
		}
		$seqs.=$_;
	}
}
&Assembly($seqN,$seqNs,$seqsDir,$seqs);
print "All jobs finished at:".time()."\n";

sub Assembly(){
	my ($seqname,$seqNs,$seqsdir,$sequences)=@_;
	my $seqnum=$sequences=~s/\>/\>/g;
	if($seqnum<2){
		my @seqs=split(/\n/,$sequences);
		my $sequence=dealseqid($seqname,@seqs);
		$sequence=formatseq($sequence);
		print $fh_outfile "$sequence";
		print $fh_logfile "singlets: $seqNs\n\n";
	}
	else{
		my $tempfile="$parTempdir\/temp.fas";
		open(my $fh_tempfile,">$tempfile");
		print $fh_tempfile $sequences;
		close $fh_tempfile;
		my @contigs;
		my @contigstrand;
		if ($parMethod eq "cap3"){ ##cap3 assembly 
			my $result=`cap3 $parTempdir\/temp.fas $cap3para`;
			my $contigfile="$parTempdir\/temp.fas.cap.contigs";
			open(my $confile,$contigfile) || die "Can't open file $contigfile\n";
			@contigs=<$confile>;
			close $confile;
			@contigstrand=&detectcap3strand($result);
		}
		else{
			my $result=`phrap -trim_qual 10 $parTempdir\/temp.fas`;
			my $contigfile;
			if($parMethod eq "phrap"){
				$contigfile="$parTempdir\/temp.fas.contigs";
			}
			else{
				$contigfile="$parTempdir\/temp.fas.cap.contigs";
			}
			open(my $confile,$contigfile) || die "Can't open file $contigfile\n";
			@contigs=<$confile>;
			close $confile;
			my $logfile="$parTempdir\/temp.fas.log";
			open(my $logf,$logfile) || die "Can't open file $logfile\n";
			my @logs=<$logf>;
			close $logf;
			@contigstrand=&detectphrapstrand(@logs);
		}
		if($parDirCorrection==0){
			my $contignum=scalar(@contigstrand);
			if($contignum==0){
				my @seqs=split(/\n/,$sequences);
				foreach my $i(@seqs){
					$i=$i."\n";
				}
				my $sequence=dealseqid($seqname,@seqs);
				$sequence=formatseq($sequence);
				print $fh_outfile "$sequence";
				print $fh_logfile "singlets: $seqNs\n\n";
			}
			else{
				my $sequence=dealseqid($seqname,@contigs);
				$sequence=formatseq($sequence);
				print $fh_outfile "$sequence";
				my $contigcount=1;
				foreach my $i(@contigstrand){
					next if($i eq "");
					my @contiginfo=split(/\t/,$i);
					print $fh_logfile "contig_$contigcount\t$contiginfo[0]\t$contiginfo[1]\n";
					$contigcount++;
				}
				my $singlets=getsinglets();
				if($singlets ne ""){
					print $fh_logfile "singlets: $singlets\n\n";
				}
				else{
					print $fh_logfile "\n";
				}
			}
		}
		else{
			my $contignum=scalar(@contigstrand);
			if($contignum==0){
				my @seqs=split(/\n/,$sequences);
				foreach my $i(@seqs){
					$i=$i."\n";
				}
				my $sequence=dealseqid($seqname,@seqs);
				$sequence=formatseq($sequence);
				my $seqswithdir=dealdirection($sequence,$seqsdir);
				print $fh_outfile "$seqswithdir";
				print $fh_logfile "singlets: $seqNs\n\n";
			}
			else{
				my $contigdir=dealcontigdir($seqsdir,$seqNs,\@contigstrand);
				my $sequence=dealseqid($seqname,@contigs);
				my $seqswithdir=dealdirection($sequence,$contigdir);
				print $fh_outfile "$seqswithdir";
				my $contigcount=1;
				foreach my $i(@contigstrand){
					next if($i eq "");
					my @contiginfo=split(/\t/,$i);
					print $fh_logfile "contig_$contigcount\t$contiginfo[0]\t$contiginfo[1]\n";
					$contigcount++;
				}
				my $singlets=getsinglets();
				if($singlets ne ""){
					print $fh_logfile "singlets: $singlets\n\n";
				}
				else{
					print $fh_logfile "\n";
				}
			}
		}
	}
}

sub dealcontigdir{
	my ($seqDs,$seqNs,$refcontigstrand)=@_;
	#print "%%%%$seqDs\t$seqNs\n";
	my @contigstrands=@{$refcontigstrand};
	my @seqDs=split(/,/,$seqDs);
	my @seqNs=split(/,/,$seqNs);
	my %seqsdir;
	for(my $k=0;$k<@seqNs;$k++){
		$seqsdir{$seqNs[$k]}=$seqDs[$k];
	}
	my $contigdir;
	foreach my $i(@contigstrands){
		next if($i eq "");
		#print "!!!$i\n";
		my @strandinfo=split(/\t/,$i);
		my @ids=split(/,/,$strandinfo[0]);
		my @dirs=split(/,/,$strandinfo[1]);
		my $transscore=0;
		my $cisscore=0;
		for(my $j=0;$j<@ids;$j++){
			#print "$ids[$j]\t$dirs[$j]:$seqsdir{$ids[$j]}\n";
			if($dirs[$j] ne $seqsdir{$ids[$j]}){
				#print "AA\n";
				$transscore++;
			}
			else{
				#print "BB\n";
				$cisscore++;
			}
		}
		if($cisscore>=$transscore){
			if($contigdir eq ""){
				$contigdir="+";
			}
			else{
				$contigdir.=",+";
			}
		}
		else{
			if($contigdir eq ""){
				$contigdir="-";
			}
			else{
				$contigdir.=",-";
			}
		}
	}
	#print "contigdir:$contigdir\n";
	return $contigdir;
}

sub formatseq{
	my $Seq=shift;
	#print "SEQ:\n$Seq\nEND\n";
	my @Seqs=split(/\n/,$Seq);
	my $formatedseq;
	my $seqNa;
	my $seq;
	for(my $i=0;$i<@Seqs;$i++){
		if($Seqs[$i]=~/^>(\S+)/){
			if($seq ne ""){
				my $modulus = length($seq) % 80;
				$seq=~s/(.{80})/$1\n/ig;
				if ($modulus!=0) {
					$seq=$seq."\n";
				}
				$formatedseq.=">$seqNa\n$seq";
				#print "FORMATED:\n$formatedseq\END\n";
			}
			$seqNa=$1;
			$seq="";
		}
		else{
			$seq.=$Seqs[$i];
			#print "SEQ:$seq\n";
		}
	}
	if($seq ne ""){
		my $modulus = length($seq) % 80;
		$seq=~s/(.{80})/$1\n/ig;
		if ($modulus!=0) {
			$seq=$seq."\n";
		}
		$formatedseq.=">$seqNa\n$seq";
	}
	return $formatedseq;
}


sub dealdirection{
	my ($Seq,$Dirs)=@_;
	#print "$Dirs\n$Seq\n";
	my @seqs=split(/\n/,$Seq);
	my @Dirs=split(/,/,$Dirs);
	my $seqcount=-1;
	my $seq;
	my $contigs;
	my $seqname;
	for(my $i=0;$i<@seqs;$i++){
		$seqs[$i]=~s/\n//g;
		if($seqs[$i]=~/^>(\S+)/){
			if($seq ne ""){
				if($Dirs[$seqcount] eq "-"){
					$seq=dna_reverser($seq);
				}
				my $modulus = length($seq) % 80;
				$seq=~s/(.{80})/$1\n/ig;
				if ($modulus!=0) {
					$seq=$seq."\n";
				}
				$contigs.=">$seqname\n$seq";
			}
			$seqname=$1;
			$seq="";
			$seqcount++;
		}
		else{
			$seq.=$seqs[$i];
		}
	}
	if($seq ne ""){
		if($Dirs[$seqcount] eq "-"){
			$seq=dna_reverser($seq);
		}
		my $modulus = length($seq) % 80;
		$seq=~s/(.{80})/$1\n/ig;
		if ($modulus!=0) {
			$seq=$seq."\n";
		}
		$contigs.=">$seqname\n$seq";
	}

	return $contigs;
}

sub dna_reverser {
    my($Seq) = @_;
    my $Rev_Seq = reverse($Seq);
    $Rev_Seq =~ tr/[atgc]/[tacg]/;
    $Rev_Seq =~ tr/[ATGC]/[TACG]/;
    return($Rev_Seq);
}

sub dealseqid{
    my @conseq=@_;
    my $locusid=shift(@conseq);
    my $count=1;
    my $fasseqs;
    for(my $i=0;$i<@conseq;$i++){
		if($conseq[$i]=~/^>\S+/){
		    $fasseqs.=">$locusid\_contig$count\n";
		    $count++;
		}
		else{
		    $fasseqs.="$conseq[$i]";
		}
    }
    return $fasseqs;
}

sub detectphrapstrand{
	my @loginfo=@_;
	my $contigNum=0;
	my $seqids;
	my $dirs;
	my $ismatch=0;
	my @dirinfo;
	my $isstart=0;
	for(my $i=0;$i<@loginfo;$i++){
		if($loginfo[$i]=~/^Contig\s(\d+):\s+\S+/){
			$isstart=1;
			#print "$loginfo[$i]\n";
			if($seqids eq ""){
				$ismatch=1;
			}
			else{
				my $dirinfo="$seqids\t$dirs";
				#print "contigNum:$contigNum\t$dirinfo\n";
				$dirinfo[$contigNum]=$dirinfo;
			}
			$contigNum=$1;
			$seqids="";
			$dirs="";
		}
		else{
			if($isstart==1){
				if($loginfo[$i]=~/^\s+\d+\s+C\s(\S+)\s+\d+/){
					#print "AA$loginfo[$i]\n";
					if($seqids eq ""){
						$seqids=$1;
						$dirs="-";
					}
					else{
						$seqids.=",$1";
						$dirs.=",-";
					}
				}
				elsif($loginfo[$i]=~/^\s+\d+\s+(\S+)\s+\d+/){
					#print "BB$loginfo[$i]\n";
					if($seqids eq ""){
						$seqids=$1;
						$dirs="+";
					}
					else{
						$seqids.=",$1";
						$dirs.=",+";
					}
				}
			}
		}
	}
	if($seqids ne ""){
		my $dirinfo="$seqids\t$dirs";
		#print "contigNum:$contigNum\t$dirinfo\n";
		$dirinfo[$contigNum]=$dirinfo;
	}
	return @dirinfo;
}

sub detectcap3strand{
    my $assinfo=shift;
    my @assinfos=split(/\n/,$assinfo);
    my @dirinfo;
    my $seqids;
    my $dirs;
    my $ismatch=0;
    #print "$assinfo\n";
    for(my $i=0;$i<@assinfos;$i++){
		if($assinfos[$i]=~/\*+\s(\S+\s\d+)\s\*+/){
		    if($seqids eq ""){
				$ismatch=1;
				next;
		    }
		    else{
				my $dirinfo="$seqids\t$dirs";
				push(@dirinfo,$dirinfo);
		    }
		    $seqids="";
		    $dirs="";
		}
		else{
		    if($ismatch==1){
				if($assinfos[$i]=~/(\S+)([\+\-])/){
				    if($seqids eq ""){
						$seqids=$1;
						$dirs=$2;
				    }
				    else{
						$seqids.=",$1";
						$dirs.=",$2";
				    }
				}
				elsif($assinfos[$i]=~/\s+(\S+)([\+\-])/){
				    if($seqids eq ""){
						$seqids=$1;
						$dirs=$2;
				    }
				    else{
						$seqids.=",$1";
						$dirs.=",$2";
				    }
				}
				elsif($assinfos[$i]=~/DETAILED DISPLAY OF CONTIGS/){
				    last;
				}
		    }
		}
    }
    if($seqids ne ""){
		my $dirinfo="$seqids\t$dirs";
		push(@dirinfo,$dirinfo);
    }
    return @dirinfo;
}

sub Filecheck{
    my ($Filename,$mode)=@_;
    if($mode==0){
	if(!(-e $Filename)){
	    print stderr "Can't find file: $Filename\n";
	    $iserror=1;
	}
    }
    else{
	if(-e $Filename){
	    print stderr "File $Filename already exists\n";
	    $iserror=1;
	}
    }
}

sub getsinglets{
	my $singlefile;
	if($parMethod eq ""){
		$singlefile="$parTempdir\/temp.fas.singlets";
	}
	else{
		$singlefile="$parTempdir\/temp.fas.cap.singlets";
	}
	open (my $fh_singlefile,$singlefile) || return "";
	my @singlets=<$fh_singlefile>;
	my $singletsN;
	for(my $i=0;$i<@singlets;$i++){
		if($singlets[$i]=~/^>(\S+)/){
		    if($singletsN eq ""){
		    	$singletsN=$1;
		    }
		    else{
		    	$singletsN.=",$1";
		    }
		}
    }
	return $singletsN;
}

sub Usage(){
  print << "    Usage";

	Usage:  $0 (version $version)

	<options>
		-i     Input file for assembly
		       file format:
		       				#seqidA +,-
		       				>seqnameA
		       				XXXXX
		       				>seqnameB
		       				XXXXX
		       				#seqidB...
		-m     Method for assembling sequences (phrap or cap3), default: phrap
		-d     Direction correction for sequences (0:off, 1:on, default:0)
		-t     Temporary director
		-o     Output file

    Usage
	exit(1);
}
