for i in `cat list `; do ping -c 40 $i >$i.ping & done
perl stat_ping.pl *ping
