#!/usr/bin/python
# Author:Finlay Maguire
# Date:09/11/2012

# Usage: python logdet_replicate_p4.py -i input file

# import all of p4 module
from p4 import *
import sys
from optparse import OptionParser

# commandline options
parser = OptionParser()

parser.add_option("-i", "--input",
                  action="store",
                  type="string",
                  dest="inputfilename",
                  help="Input file containing alignment")

parser.add_option("-v", "--verbose",
                  action="store_true",
                  dest="verbose",
                  default=False,
                  help="Highly verbose stdout for debugging, default=False")

parser.add_option("-c", "--correction",
                  action="store",
                  type="string",
                  dest="correction",
                  default="TK02",
                  help="Matrix correction - TK02 or L94, default=TK02")

parser.add_option("-m", "--missing_char",
                  action="store",
                  type="string",
                  dest="missing_char_strat",
                  default="fudge",
                  help="Missing characters - refuse or fudge, default=fudge")

parser.add_option("-d", "--nonpostive_determiant",
                  action="store",
                  type="string",
                  dest="nPosDetStrat",
                  default="invert",
                  help="Non-positive determinant - refuse or invert, default=invert")

parser.add_option("-g", "--generations",
                  action="store",
                  type="int",
                  dest="generations",
                  default=1000,
                  help="Generations to run, default=1000")

(options, args) = parser.parse_args()

# validate the input
validate=0

if options.generations<=0:
    print "\n\n-----ERROR-----\n\n"
    print "Generations must be a positive integer"
    print "\n\n-----ERROR-----\n\n"
    validate=+1

if options.correction not in ['TK02','L94']:
    print "\n\n-----ERROR-----\n\n"
    print "Matrix correction must equal TK02 or L94"
    print "\n\n-----ERROR-----\n\n"
    validate=+1

if options.missing_char_strat not in ['refuse','fudge']:
    print "\n\n-----ERROR-----\n\n"
    print "Missing character strategy must be set to: refuse or fudge"
    print "\n\n-----ERROR-----\n\n"
    validate=+1

if options.nPosDetStrat not in ['refuse','invert']:
    print "\n\n-----ERROR-----\n\n"
    print "Non-positive determinant strategy must be set to: refuse or invert"
    print "\n\n-----ERROR-----\n\n"
    validate=+1

if validate >= 1:
    parser.print_help()
    sys.exit("Execution failed")

# load alignment from argument when invoking
read(options.inputfilename)
in_alignment = var.alignments[0]

# output filename
output_all_trees_nexus=options.inputfilename+"_all_trees_file.nxs"
consensus_tree_output=options.inputfilename+"_consensus_tree.nxs"

# init lists
bootstrapped_alignments=[]
log_det_distance_matrices=[]
nj_trees=[]

# init taxnames
species_list=in_alignment.taxNames

# pInvar of Constant
alignment_without_gaps_or_ambig=in_alignment.noGapsOrAmbiguitiesCopy()
alignment_pinvar_constants=alignment_without_gaps_or_ambig.constantSitesProportion()

# loop for specified generations
for index in range(options.generations):

    print "\n\nReplicate",index,":"

    # bootstrap replicate alignment
    bootstrapped_alignments.append(in_alignment.bootstrap())

    if options.verbose==True:
        print "\n\n-----ALIGNMENT-----\n\n"
        print bootstrapped_alignments[index]
        bootstrapped_alignments[index].writeNexus(fName=None)

    print "Bootstrap generated"

    # make logdet distance matrix for each replicate
    ld_mat=bootstrapped_alignments[index].logDet(
    correction=options.correction,
    missingCharacterStrategy=options.missing_char_strat,
    nonPositiveDetStrategy=options.nPosDetStrat,
    doPInvarOfConstants=True,
    pInvarOfConstants=alignment_pinvar_constants,
    )

    log_det_distance_matrices.append(ld_mat)

    if options.verbose==True:
        print "\n\n-----MATRIX-----\n\n"
        print log_det_distance_matrices[index]
        log_det_distance_matrices[index].writeNexus(fName=None)

    print "LogDet distance matrix generated"

    # build nj tree for each distance matrix
    tree=log_det_distance_matrices[index].bionj()
    nj_trees.append(tree)

    if options.verbose==True:
        print "\n\n-----PHYLOGENY----\n\n"
        print nj_trees[index]
        nj_trees[index].writeNexus(fName=None)

    print "NJ tree generated"

    # write tree to nexus for safe keeping
    nj_trees[index].writeNexus(fName=output_all_trees_nexus, append=1)


# make trees object containing all trees
all_trees_as_trees_obj=Trees(nj_trees,taxNames=species_list)

# make tree partition object from all trees
all_tree_partitions_as_tp_object=TreePartitions(all_trees_as_trees_obj)

# build consenus tree from tree partition object
consensus_tree = all_tree_partitions_as_tp_object.consensus()

# transfer support values to tree from support object
for n in consensus_tree.iterInternalsNoRoot():
    n.name = "%.0f" % (100. * n.br.support)

# output consensus tree to stdout
print "\n\n\n\n\n\n Consensus Tree: \n\n\n\n"
consensus_tree.draw()

# print consensus tree to file
consensus_tree.writeNexus(consensus_tree_output)

# print run settings to log file

run_command=str(sys.argv)

log_file_name=options.inputfilename+".log"
log_file=open(log_file_name, "w")
log_file.write("-----LogDeterminant Analysis Settings-----")
log_file.write("\n\nInput file: "+options.inputfilename)
log_file.write("\n\nCommandline input and arguments: "+run_command)
log_file.write("\n\nCorrections: "+options.correction)
log_file.write("\n\nmissingCharacterStrategy: "+options.missing_char_strat)
log_file.write("\n\nnonPositiveDetStrategy: "+options.nPosDetStrat)
log_file.write("\n\npInvarOfConstants: "+ str(alignment_pinvar_constants))
log_file.write("\n\nOutput file of all trees: "+output_all_trees_nexus)
log_file.write("\n\nOutput file of consensus tree: "+consensus_tree_output)
log_file.close()

log_file=open(log_file_name,"r")
log=log_file.read()
print "\n\n"
print log
log_file.close()


