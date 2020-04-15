#!/usr/bin/perl

use strict;
use warnings;

my $file=shift;

my %rec;
open IN,"< $file" or die $!;
while(<IN>){
        chomp;
        my $ori_name=$_;
        my ($path,$filename)=($1,$2) if($ori_name=~/^(.*)\/([^\/]+)$/);
#       chop $path;
        my ($tag,$temp,$lane,$mate)=($1,$2,$3,$4) if($filename=~/^([^\_]+)\_(.*)\_(L\d+)\_R(\d)/);
        my $new_name="$path/130902_SN350_$temp\_$lane\_$tag\_$mate.fq.gz";
        #my $new_name="$path/130226_SN350_$temp\_$lane\_$tag\_$mate.fq.gz";
        `ln -s $ori_name  $new_name`;
        push @{$rec{$tag}}, $new_name;
}
close IN;
foreach my $key(keys %rec){
        my @array=@{$rec{$key}};
        #my $path=$1 if($array[0]=~/^(.*)\/130226_SN350_/);
        my $path=$1 if($array[0]=~/^(.*)\/130902_SN350_/);
       chdir $path;
        my $insertsize=$1 if($key=~m/(\d+)/);
        #$insertsize-=120;
        open OUT,">  $key.lib.lst" or die $!;
        print OUT "$key $insertsize";
        close OUT;
        open OUT,"> $key.lane.lst" or die $!;
        foreach(my $i=0;$i<@array;$i+=2){
                print OUT "$array[$i] 0 0 40\n$array[$i+1] 0 0 10\n";
        }
        close OUT;
        `perl /lustre/user/dongy/Mybin/Assembly/Filter_data/run_filter.pl   $key.lane.lst  $key.lib.lst  `;
        `chmod 755 $key.lane.lst.filter.pbs`;
#       `perl  /lustre/user/dongy/Mybin/common-bin/submit.pl  $key.lane.lst.filter.pbs  long`;
}
