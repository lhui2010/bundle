sudo yum install autofs
sudo vim /etc/auto.master
sudo vim /etc/fstab
sudo vim /etc/auto.master
sudo vim /etc/auto.nfsdb
sudo systemctl start autofs.service
cd /nfs/
ls
ll
ll liuhui/
systemctl enable autofs.service
sudo systemctl enable autofs.service

