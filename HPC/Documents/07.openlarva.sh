cd /opt/openlava-2.2/
ll
cd etc/
ll
vi lsf.cluster.openlava
sudo vi lsf.cluster.openlava
lsid
sudo vi /etc/hosts
##
#192.168.119.10  yi01
# 192.168.119.9   yi02
# 192.168.119.11  yi03
# 192.168.119.12  yi04
##
ll
vi lsf.cluster.openlava
ping yimiddle
more /etc/sysconfig/iptables
more /etc/sysconfig/iptables-config
sudo iptables -I INPUT -s 192.168.119.0/28 -p udp --match multiport --dports 1024:65535 -j ACCEPT
sudo iptables -I INPUT -s 192.168.119.0/28 -p tcp --match multiport --dports 1024:65535 -j ACCEPT
sudo iptables-save|sudo tee /etc/sysconfig/iptables
ls
/etc/init.d/openlava restart
sudo /etc/init.d/openlava restart
hostnamectl status
/etc/init.d/openlava restart

