#!/bin/bash
# usage: ./clone_user.sh $USERNAME

usr=$1
url="https://api.github.com/users/${usr}/repos"
num=$(curl -sI "$url?page=1&per_page=100" | sed -nr 's/^Link:.*page=([0-9]+)&per_page=100>; rel="last".*/\1/p')

for i in $(seq 1 $num)
do 
    for j in $(curl -s "$url?page=${i}&per_page=100" | grep "clone_url"  | sed 's/.*clone_url\"://' | tr -d "\"" | tr -d "," | sed 's/https:\/\//git@/' | sed 's/github.com\//github.com:/')
    do 
        git clone $j
    done
done  



#| sed -nr 's/.*clone_url": "(.*)",/\1/p')
