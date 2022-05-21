
$_=<>;
print;
while(<>)
{
    s/\n//;
    s/\r//;
    s/"//g;
    my @e=split/\t/, $_;
    print $e[0];
    shift @e;
    for my $e(@e)
    {
        my $count = 0;
        $e=~s/\$//;
        if($e ne "")
        {
            $count ++;
            my $c=()= $e =~ /,/g;
            $count +=$c;
        }
        print "\t$count";
    }
    print "\n";
}
