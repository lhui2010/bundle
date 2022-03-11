#sudo systemctl enable firewalld
sudo firewall-cmd --new-zone=lsf --permanent
sudo firewall-cmd --zone=lsf --add-source=192.168.119.0/28  --permanent 
sudo firewall-cmd --zone=lsf --add-port=1024-65535/tcp --permanent 
sudo firewall-cmd --zone=lsf --add-port=1024-65535/udp --permanent 
