#!/opt/anaconda/bin

import glob
import argparse
import subprocess
import os
from Bio import SeqIO
from Bio.Blast import NCBIXML


def system_call(command_string):
    '''wrapper for the system call that catches
    non-zero returncodes logs and raises an error'''
    with open('/dev/null', 'w') as devnull:
        # devnull to suppress all incidental process output
        return_code = subprocess.call(command_string.split(),
                                      stderr=devnull,
                                      stdout=devnull)
        if return_code != 0:
            assert False


def parse_blast(unparsed_blast_hit_file, proteome):
    hit_id_set = set()
    hit_records = []
    blast_hits = NCBIXML.parse(open(unparsed_blast_hit_file, 'r'))
    # returns an iterator of blast record objects over the xml blastoutput file
    for blast_record in blast_hits:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                hit_id_set.add(alignment.hit_def)
    for hit_id in hit_id_set:
        # search proteome for sequence
        proteome_parse = SeqIO.parse(open(proteome+'.fas', 'r'), 'fasta')
        proteome_sequences = {seqrec.description: seqrec for seqrec
                              in proteome_parse}
        sequence_record = proteome_sequences[hit_id]
        hit_records.append(sequence_record)
    return hit_records


def blast(input_seed_seq, mode, proteome_paths, evalue):
    '''Blast seed sequence against db using BLASTP'''
    parsed_output_fh = open(seed_file+'_output', 'a')
    for proteome in proteome_paths:
        # for each of the proteomes we want to blast against
        unparsed_blast_hit_file = "{0}_v_{1}.bo".format(str(input_seed_seq),
                                                        str(proteome).replace(
                                                            '/',
                                                            '_')
                                                        )
        blast_cmd = "{4} -query {0} -db {1} -out {2} -evalue {3} -outfmt 5".format(
                                     input_seed_seq,
                                     proteome,
                                     unparsed_blast_hit_file,
                                     evalue,
                                     mode)
        # blast a seed seq against that proteome

        system_call(blast_cmd)

        hits = parse_blast(unparsed_blast_hit_file, proteome)
        blast_summary_file = "{0}_{1}_hit_summary".format(input_seed_seq, mode)
        with open(blast_summary_file, 'a') as blast_summary_fh:
            num_hits = str(len(hits))
            blast_summary_fh.write("{0}_hits: {1}\n".format(proteome, num_hits))
        SeqIO.write(hits, parsed_output_fh, 'fasta')
        os.unlink(unparsed_blast_hit_file)

    parsed_output_fh.close()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Pipeline to rapidly generate \
                                     phylogenies from a list of seed \
                                     sequences")

    parser.add_argument("--protein_seeds", "-p",
                        dest='prot_seeds',
                        type=str,
                        default=None,
                        help="File containing all the protein sequences \
                             desired in fasta format")

    parser.add_argument("--nucleotide_seeds", "-n",
                        dest='nucl_seeds',
                        type=str,
                        default=None,
                        help="File containing all the nucleotide seed \
                             sequences desired in fasta format")

    parser.add_argument('--proteome_folder', '-g',
                        dest='proteomes',
                        type=str,
                        required=True,
                        help='Folder containing the proteomes and blastdbs \
                              you want to search')

    parser.add_argument('--output_folder', '-o',
                        dest='out_folder',
                        type=str,
                        default=os.getcwd(),
                        help='Folder to output results')

    parser.add_argument('--evalue', '-e',
                        dest='evalue',
                        type=float,
                        default=1e-4,
                        help='E-value to use in BLAST searches')

    args = parser.parse_args()

    if not bool(args.prot_seeds) ^ bool(args.nucl_seeds):
        assert False, "can't run both modes at the same time"
    if bool(args.nucl_seeds):
        mode = 'blastx'
        seeds = args.nucl_seeds
    if bool(args.prot_seeds):
        mode = 'blastp'
        seeds = args.prot_seeds

    proteome_glob = glob.glob(args.proteomes+'/*.fas')
    proteome_paths = [path.replace('.fas', '') for path in proteome_glob]

    for seq in SeqIO.parse(open(seeds, 'r'), 'fasta'):
        seq_name = seq.description.replace(' ', '_')
        seed_file = '{0}/{1}.fas'.format(args.out_folder, seq_name)
        SeqIO.write(seq, open(seed_file, 'w'), 'fasta')
        blast(seed_file, mode, proteome_paths, str(args.evalue))
        os.unlink(seed_file)
