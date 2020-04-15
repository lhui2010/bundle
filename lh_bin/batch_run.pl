#!env perl
use warnings;
use Cwd;

#DB table 
my %GENOME_HASH=(
    "maize" => "/lustre/local/database/GENOME/zea_mays_agp_v4/Zea_mays.AGPv4.dna_sm.toplevel.fa",
    );
my %GTF_HASH=(
    "maize" => "/lustre/local/database/GENOME/zea_mays_agp_v4/Zea_mays.AGPv4.32.gff3",
    );



#Running Config
my $threads=8;



my $dir = getcwd;
my @list = `ls -d */`;

for my $user (@list)
#for my $user (qw/test/)
{
    my (%left, %right);
    my $sp_name=`ls -d $user/*/ |grep input_`;
    $sp_name=~s/.*input_//;
    $sp_name=~s/\/\n//;
    #print $sp_name;exit;
    my $GM=$GENOME_HASH{$sp_name};
    my $GTF=$GTF_HASH{$sp_name};

    my $PE_PATH="$dir/$user/input_$sp_name";
    #print $PE_PATH;exit;
    chdir ($PE_PATH) or die "cannot change: $!\n";

    my $cdir = getcwd;
    #print $cdir; exit;

    open IN, "../name.txt" or die;
    while(<IN>)
    {
        my ($sampleID, $newID) = split;
        my $tmp=$sampleID."_";
        $left{$newID} = `ls $tmp*_1.clean.fq.gz`;
        $right{$newID} = `ls  $tmp*_2.clean.fq.gz`;
        $left{$newID} =~s/\n//;
        $right{$newID}=~s/\n//;
        print ("qsub -cwd -j n -b y -V -N $newID.tuxedo -pe smp $threads  \"hisat2 --mp 3,1 -p $threads -x $GM -1 $PE_PATH/$left{$newID} -2 $PE_PATH/$right{$newID} -S $newID.sam   && samtools sort -@ $threads -o $newID.bam $newID.sam && stringtie -e -b $newID-st-abun-out -p $threads -G $GTF -o $newID.with_novel.gtf -A $newID.exp_table $newID.bam \"") ;
        #system ("qsub -cwd -j n -b y -V -N $newID.tuxedo -pe smp $threads  \"hisat2 --mp 3,1 -p $threads -x $GM -1 $PE_PATH/$left{$newID} -2 $PE_PATH/$right{$newID} -S $newID.sam   && samtools sort -@ $threads -o $newID.bam $newID.sam && stringtie -e -b $newID-st-abun-out -p $threads -G $GTF -o $newID.with_novel.gtf -A $newID.exp_table $newID.bam \"") ;
    }
}

