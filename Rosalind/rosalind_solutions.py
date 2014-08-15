#!/usr/bin/python

import sys

def read_in_fasta_sequences(input_file,input_file_format):
    """Function that reads in a file of sequences
    Using biopython SeqIO and returns a list of SeqIO seq record objects"""

    import Bio.SeqIO

    sequence_record_list=[]

    for seq_record in Bio.SeqIO.parse(input_file,input_file_format):
        sequence_record_list.append(seq_record)
    return sequence_record_list

def consensus_sequence(seq_record_list):
    """Function that takes in sequence list of equal length (from read_in_fasta_sequences_function)
    and builds a character matrix out of them to determine the consensus sequence (based on most common
    nt at that posistion)"""

    running_max_index=0
    acount_list, ccount_list, gcount_list, tcount_list, raw_seq_list, character_Matrix, index_consensus, consensus_string_list=[],[],[],[],[],[],[],[]

    for index in range(len(seq_record_list)):
        raw_seq_list.append(seq_record_list[index].seq)
    for sequence in raw_seq_list:
        character_Matrix.append([character for character in sequence])

    for j in range(len(character_Matrix[0])):
        acount, tcount, gcount, ccount = 0,0,0,0
        for i in range(len(character_Matrix)):
            if character_Matrix[i][j]=="T":
                tcount+=1
            elif character_Matrix[i][j]=="A":
                acount+=1
            elif character_Matrix[i][j]=="C":
                ccount+=1
            elif character_Matrix[i][j]=="G":
                gcount+=1
            else:
                print "Error building matrix"
        acount_list.append(acount)
        ccount_list.append(ccount)
        gcount_list.append(gcount)
        tcount_list.append(tcount)

    profile=[acount_list,ccount_list,gcount_list,tcount_list]

    for j in range(len(profile[0])):
        running_max=profile[0][j]
        for i in range(len(profile)):
            if profile[i][j]>=running_max:
                running_max=profile[i][j]
                running_max_index=i
            if i==(len(profile)-1):
                index_consensus.append(running_max_index)

    for i in range(len(index_consensus)):
        if index_consensus[i]==0:
            consensus_string_list.append("A")
        elif index_consensus[i]==1:
            consensus_string_list.append("C")
        elif index_consensus[i]==2:
            consensus_string_list.append("G")
        elif index_consensus[i]==3:
            consensus_string_list.append("T")

    consensus_string=(''.join(consensus_string_list))
    #print "A:",(' '.join(map(str,profile[0])))
    #print "C:",(' '.join(map(str,profile[1])))
    #print "G:",(' '.join(map(str,profile[2])))
    #print "T:",(' '.join(map(str,profile[3])))

    return(consensus_string)


def directed_adjacency_graph(sequences):
    """Build a directed adjacency graph from a list of seqrecord objects"""

    sequence_dictionary={}
    adjacency_graph_list=[]

    for index in range(len(sequences)):
        sequence_dictionary.update({sequences[index].id:str(sequences[index].seq)})

    for keys in sequence_dictionary:
        for keys2 in sequence_dictionary:
            overlap=min(len(sequence_dictionary[keys]),len(sequence_dictionary[keys]))

            while overlap>0:
                if (sequence_dictionary[keys])[-overlap:] == (sequence_dictionary[keys2])[:overlap]:
                    break
                overlap -= 1

            if overlap > 1 and keys != keys2:
                adjacency_tuple=(keys, keys2)
                adjacency_graph_list.append(adjacency_tuple)

    return adjacency_graph_list

def hamming_distance(sequences):
    """Calculates the hamming distance between sequences"""

    import itertools

    hamming_distance_list=[]

    pairwise_combinations=list(itertools.combinations(sequences,2))
    for index in range(len(pairwise_combinations)):
        sequence_record_a=pairwise_combinations[index][0]
        sequence_record_b=pairwise_combinations[index][1]
        hamm=sum(ch1 != ch2 for ch1, ch2 in zip(sequence_record_a.seq, sequence_record_b.seq))
        hamming_tuple=(sequence_record_a.id, sequence_record_b.id, hamm)
        hamming_distance_list.append(hamming_tuple)

    return hamming_distance_list

if __name__=="__main__":
    sequences=read_in_fasta_sequences(sys.argv[1],"fasta")
    print(consensus_sequence(sequences))
