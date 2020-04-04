#!/usr/bin/perl -w
print <<"usage" if (@ARGV < 2);

                =======================
                ||                   ||
                ||  select_fasta.pl  ||
                ||                   ||
                ||      by: liuhui   ||
                =======================

FUNCTION:               select target section from fasta file

HOWTO:                  .pl limition_list fasta_file >result

Thu Mar 10 15:41:40 CST 2011

usage

exit if (@ARGV <2);


$list_name=$ARGV[0];
#$fasta_name=$ARGV[1];
my %list;
$list{$list_name}=1;
#open FASTA, $fasta_name or die;

open IN, $ARGV[1] or die;
while(<IN>)
{
        if(/>/)
        {
                chomp;
                if(defined $fa_info and &iswithinlist_hash )#exists $list{$fa_info} )
                {
                        print $fa_content;
                }

                $fa_info = $_;
                $fa_content = $fa_info."\n";
        }
        else
        {
                $fa_content.=$_ ;#if(exists $list{$fa_info});
        }
}

print $fa_content if (&iswithinlist);#exists $list{$fa_info});

sub iswithinlist_hash
{
        return 0 if ($fa_info eq "");

        my $tmp = (split /\s+/, $fa_info)[0];

        $tmp =~ s/>//;

        if(exists $list{$tmp})
        {
                return 1;
        }
        else
        {
                return 0;
        }
}

sub iswithinlist
{
        my @e=split /\s+/, $fa_info;
        my $tmp=substr ($e[0], 1);
        if(exists $list{$tmp})
        {
                return 1;
        }
        return 0;
}
