#!/bin/bash

base=$(pwd)
echo "" > out.html
for i in Sample*;
    do j=${i}/trimmomatic_output; 
    cd $j;
        for zipped_fastqc in *.zip;
            do unzip -o -qq $zipped_fastqc;
        done;
        for folder in *fastqc;
            do echo "<h2><a href=\"$PWD/${folder}.html\">$folder</h2>" >> $base/out.html 
                ls -1 $PWD/$folder/Images/*.png | perl -ne 'chomp; printf "<img src=\"%s\" width=\"400px\" height=\"400px\"\>\n" , $_' >> $base/out.html;
        done; 
    cd $base;
done 
