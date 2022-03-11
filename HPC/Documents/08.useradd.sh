#!/bin/bash

set -euo pipefail

USAGE="Usage:\n\t$0 userid disk_soft_limit disk_hard_limit \nEg:\n\t./add_user.sh xiahouzuoying yitingshuang 2T 3T "

Group='yitingshuang'
NfsQuota='2T'
NfsHardQuota='3T'
# IQUOTA='10^8'
# IHQUOTA='10^8'

if [ "$#" -eq 0 ]; then
    echo -e $USAGE
    exit
elif [ $# -eq 1 ]; then
    echo "Using default group (yitingshuang) and default quota (2T 3T 10^8 10^8)"
    User=$1
elif [ $# -eq 2 ]; then
    User=$1
    Group=$2
    echo "Using group $Group and default quota (2T 3T 10^8 10^8)"
elif [ $# -eq 4 ]; then
    User=$1
    Group=$2
    NfsQuota=$3
    NfsHardQuota=$4
else
    echo -e $USAGE
    exit
fi


if [ `getent group |grep ${Group}` ]
then
        echo "Group existing, continueing.."
else
        echo "Creating Group $Group"
        sudo groupadd $Group
        sleep 1s
fi

echo "Adding User $User to Group $Group"
sleep 1s
sudo useradd $User -G $Group
echo "Setting Password"
Password=`openssl rand -base64 16`
echo "${Password}" | sudo passwd "$User" --stdin 


echo "Synchronizing"
#不再同步/etc/shadow，仅通过sge提交任务
for host in yi03 yi04
do
    sudo scp /etc/passwd ${host}:/etc/passwd
    sudo scp /etc/group ${host}:/etc/group
done

for host in yi02
do
    sudo scp -P 201 /etc/passwd ${host}:/etc/passwd
    sudo scp -P 201 /etc/group ${host}:/etc/group
done


#改权限
# chown ${User}:${User} /lustre/home/${User} # mkdir /lustre/home/${User}/.ssh
# chown ${User}:${User} /lustre/home/${User}/.ssh
# touch /lustre/home/${User}/.ssh/authorized_keys
# chmod 600 /lustre/home/${User}/.ssh/authorized_keys


#500M, 1024 block size = 524288
#1G, 1048576
#setquota -u -F vfsv0 bob 1048576 1048576 1000000 1000000 /
#首先设置本地quota
echo "Setting Quota"
sudo ssh yi04 "xfs_quota -x -c \"limit bsoft=$NfsQuota bhard=$NfsHardQuota $User\" /nfs"

echo "Add user ${User} successfully"
echo "Password is ${Password}"

echo "${User}，
你好，你的超算账号已经开通，登陆方式如下：
ssh地址: 192.168.119.10/210.72.89.32
ssh端口: 201
用户名： ${User}
密码： ${Password}

注意事项：
超算账号仅限本人使用，严禁外借、泄漏登陆方式，作业使用LSF调度系统提交，禁止在登陆节点直接运行
使用时请遵守超算管理规定，共同爱护超算的环境 

祝好!
"


