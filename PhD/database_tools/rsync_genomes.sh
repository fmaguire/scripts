#!/bin/bash

rsync -n -zPav --delete midgard:~/Desktop/genomes/cider/ /genomes/cider

read -p "Does the above look right? " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
    rsync -zPav --delete midgard:~/Desktop/genomes/cider/ /genomes/cider
    ls -1 /genomes/cider/*.fas.psq > species_list.txt
    sed -i 's/.fas.psq//g' species_list.txt
    sed -i 's/\/genomes\/cider\///g' species_list.txt
    sort species_list.txt > species_list_sorted.txt
    mv species_list_sorted.txt species_list.txt
    ./get_taxonomy_for_species_list.py species_list.txt taxonomy_list.txt
fi

