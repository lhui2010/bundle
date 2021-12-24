
while [ `pgrep scp` ]
do
    echo scp online
    sleep 1s
done

# if [ `pgrep scp2` ]
# then
#     echo scp online
# else
#     echo scp offline
# fi
