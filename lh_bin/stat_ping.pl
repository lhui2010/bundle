#--- logozy.info ping statistics ---
#2 packets transmitted, 2 packets received, 0.0% packet loss
#round-trip min/avg/max/stddev = 227.903/229.030/230.157/1.127 ms
#

for my $f (@ARGV)
{
    open IN, $f or die;
    while(<IN>)
    {
        if(/ (\d+.\d%) packet loss/)
        {
            $pack = $1;
        }
        if(/= .*\/(.*)\/.*\/.* ms/)
        {
            $ms = $1;
        }
    }

    print "$f\t$ms\t$pack\n";
}
