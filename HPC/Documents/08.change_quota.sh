#!/bin/bash

set -euo pipefail

USAGE="Usage:\n\t$0 userid disk_soft_limit disk_hard_limit \nEg:\n\t$0 xiahouzuoying  2T 3T "

Group='yitingshuang'
NfsQuota='2T'
NfsHardQuota='3T'
# IQUOTA='10^8'
# IHQUOTA='10^8'

if [ "$#" -eq 0 ]; then
    echo -e $USAGE
    exit
elif [ $# -eq 3 ]; then
    User=$1
    NfsQuota=$2
    NfsHardQuota=$3
else
    echo -e $USAGE
    exit
fi
echo "Setting Quota"
sudo ssh yi04 "xfs_quota -x -c \"limit bsoft=$NfsQuota bhard=$NfsHardQuota $User\" /nfs"

echo "Quota for $User set to soft limit of $NfsQuota and hard limit of $NfsHardQuota"
