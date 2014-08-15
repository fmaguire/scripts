#!/usr/bin/perl -w
use strict;

print "\n\nThis is a small script for extracting sequences with a specified numbers of predicted transmembrane helices from a tmhmm-2.0c output file\n\n";

print "\n\n tmhmm-2.0c output filename?\n\n";
my $filename=<>;

chomp $filename;

open FILE, "<$filename" or die $!;
open OUT, ">$filename.transmembrane" or die $!;

print "\n\n What is your minimum number of predicted helices? (Input a positive integer greater than 0)\n\n";
my $helix_number=<>;
chomp $helix_number;

if ($helix_number =~ /^\d+$/ ) { 
        print "\n\n Searching for sequences with ", $helix_number, " predicted helices\n\n";
    } else { 
            print "\n\n Sorry input must be a positive integer greater than 0\n\n";
    die; 
    }

my $helices;
my @file=<FILE>;
my $file;
foreach $file(@file) {
        if ($file =~ m/.*?(PredHel)(=)(\d+)/is){
                $helices=$3;
                  if ($helices>=$helix_number){
                    print OUT $file;
                  } 
        }
}
chomp $filename;
print "\n\n Results outputted into ", $filename,".transmembrane\n\n";

close FILE;
close OUT;
