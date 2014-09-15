#!/bin/bash

rsync -n -zPav --size-only --delete midgard:~/Desktop/genomes/cider/ /genomes/cider

read -p "Does the above look right? " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
    echo "Retrieving files"
    rsync -zPav --size-only --delete midgard:~/Desktop/genomes/cider/ /genomes/cider
    echo "Generating species list"
    ls -1 /genomes/cider/*.fas.psq | grep -oP '(?<=\/genomes\/cider\/).*?(?=\.fas\.psq|$)' \
        > species_list.txt
    echo "Sorting species list"
    sort --parallel=10 species_list.txt > species_list_sorted.txt
    mv species_list_sorted.txt species_list.txt
    echo "Generating taxonomy list"
    split -n l/10 species_list.txt split_
    for i in split_*; do
        ./get_taxonomy_for_species_list.py ${i} ${i}_taxonomy_list.txt &
    done
    wait
    cat split_*_taxonomy_list.txt > taxonomy_list.txt
    rm split_*
    sort --parallel=10 taxonomy_list.txt > taxonomy_list_sorted.txt
    mv taxonomy_list_sorted.txt taxonomy_list.txt
fi

