#!/opt/anaconda/bin/python

import time
import urllib2
from HTMLParser import HTMLParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from BeautifulSoup import BeautifulSoup

def possible_url_generator_builder(leader_id_string, database_id, main_url):
    '''Produces a generator for 00000-99999 with the appropriate
    sequence_id_string'''

    possible_urls = ('{0}{1}{2:05d}{3}'.format(main_url,
                                               leader_id_string,
                                               x,
                                               database_id) for x in range(1, 35000)) #40000
    return possible_urls

def query_and_retrieve_sequences(possible_urls):

    fail_index = 0
    sequences = []
    for url in possible_urls:

        try:
            page = str(urllib2.urlopen(url).read())
        except:
            print "Failed to open URL: {0}".format(url)
            time.sleep(60)

        # open page and convert to str
        fail_index, possible_seq = parse_page(page, fail_index)

        if fail_index == 5:
            print "Failed to retrieve 5 consecutive sequences, assuming max number for database has been reached"
            return sequences

        if possible_seq is not None:
            split_seq = possible_seq.split()
            seq_id = split_seq[0].replace('>', '')
            print seq_id
            seq_desc = split_seq[1]
            seq = ''.join(split_seq[2:])
            seq_rec = SeqRecord(Seq(seq, IUPAC.protein), id=seq_id, description=seq_desc)
            sequences.append(seq_rec)
        else:
            continue

    return sequences

def parse_page(page, fail_index):
    soup = BeautifulSoup(page)
    try:
        possible_seq = soup.find('pre').getText('pre')
    except:
        fail_index += 1
        print soup
        print "Number of failures: {0}".format(fail_index)
        return fail_index, None

    return fail_index, possible_seq

def write_sequences(sequences, output_fh):
    SeqIO.write(sequences, output_fh, 'fasta')


if __name__=='__main__':

    main_url = 'http://paramecium.cgm.cnrs-gif.fr/cgi/tool/get_sequence?dbtype=protein&id='

    database_information = {'PPRIMP' : '&db=primaurelia_Ir4-2_annotation_v1.protein',
                            'PSEXPNG' : '&db=sexaurelia_AZ8-4_annotation_v1.protein',
                            'PCAUDP' : '&db=caudatum_43c3d_annotation_v1.protein',
                            'PMMNP' : '&db=multimicronucleatum_MO3c4_annotation_v1.protein',
                            'PBIGNP' : '&db=biaurelia_V1-4_annotation_v1.protein'}

    index = 0
    for id_leader, database_id in database_information.items():

        possible_urls = possible_url_generator_builder(id_leader, database_id, main_url)

        sequences = query_and_retrieve_sequences(possible_urls)

        write_sequences(sequences, id_leader+'.fas')

        index += 1
        print 'database searched {0}, {1} of {2}'.format(id_leader, index, len(database_information))



