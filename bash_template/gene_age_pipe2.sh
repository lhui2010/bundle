#!/bin/bash
cat $1 | while read line
do
    DIR_NAME=`echo $line |awk '{print $1}'`
    GENES=`echo $line | awk '{print $2}'`

    # touch ${DIR_NAME}
    # rm -rf ${DIR_NAME}
    mkdir -p ${DIR_NAME}
    pushd ${DIR_NAME}
    bsub -q Q104C512G_X4 -o output.%J -e error.%J -J ${DIR_NAME} -n 20 "cd $PWD && python -m iga.project.nitfix correct_gene_age --threads 20 ${GENES} "
    sleep 1s
    # echo dir:$DIR_NAME
    # echo gnees:$GENES
    popd
    #deb ug
    #break
done
