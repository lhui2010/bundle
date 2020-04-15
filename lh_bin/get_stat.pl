 use File::Basename;

my $count = 1;
while(<>)
{
chomp;
my $gz_name = $_;
my $stat_name = $gz_name.".stat";
#my $dir_name = dirname($gz_name);
open SHELL, ">ol.qual.$count.sh" or die;
print SHELL "#!/bin/sh
#PBS -N qual_stat.$count
#PBS -p 1023
#PBS -q longlong
#PBS -l nodes=1:node-32:ppn=1

echo Working directory is $PBS_O_WORKDIR
echo Start date is `date`
gzip -dc $gz_name |perl /lustre/user/dongy/project/HJB/00.data/Sample_0-HJB/shell/check_qual.pl >$stat_name;
echo End date is `date`
";
close SHELL;

$count++;
}
