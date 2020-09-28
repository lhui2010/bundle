
for i in `seq -w 02 05`
    do
        #bsub -o output.%J -e error.%J -q Q64C1T_X4 -m node$i "rm -fr /tmp/lh_config &&  cp -fr  /ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/config /tmp/lh_config"
        #bsub -o output.%J -e error.%J -q Q64C1T_X4 -m node$i "rm -fr /tmp/lh_config && cp -fr  /ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/config /tmp/lh_config"
        bsub  -q Q64C1T_X4 -m node$i "touch /tmp/lh_config && rm -fr /tmp/lh_config && cp -fr  /ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/config /tmp/lh_config"
        #bsub  -q Q64C1T_X4 -m node$i "rm -fr /tmp/maker_*"
    done

for i in `seq 10 13`
    do
        #bsub -o output.%J -e error.%J -q Q104C512G_X4 -m node$i "rm -fr /tmp/lh_config && cp -fr  /ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/config /tmp/lh_config"
        bsub  -q Q104C512G_X4 -m node$i "touch /tmp/lh_config && rm -fr /tmp/lh_config && cp -fr  /ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/config /tmp/lh_config"
        #bsub  -q Q104C512G_X4 -m node$i "rm -fr /tmp/maker_*"
        #bsub -o output.%J -e error.%J -q Q104C512G_X4 -m node$i "rm -fr /tmp/lh_config &&  cp -fr  /ds3200_1/users_root/yitingshuang/lh/bin/maker3/exe/augustus-3.3.3/augustus-3.3.3/config /tmp/lh_config"
    done
#export AUGUSTUS_CONFIG_PATH=/tmp/lh_config

