#!/bin/zsh

watch_file=$1
cmd_to_run=$2

while true    
    do
        hash_1=`md5sum ${watch_file}`
        if [[ "$hash_1" != "$hash_2" ]]
            then    
                eval $cmd_to_run
                hash_2=$hash_1
        fi
    sleep 5
done
