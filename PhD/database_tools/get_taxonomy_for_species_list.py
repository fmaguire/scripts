#!/opt/anaconda/bin/python

import sys
from Bio import Entrez
import Queue
import threading

def parse_species_list(species_list_file):
    with open(species_list_file, 'r') as fh:
        species_list = [x.replace('_', ' ').rstrip() for x in fh.readlines()]
    return species_list

"""j
def fetch_parallel(species_names):
    results_queue = Queue.Queue()

    threads = [threading.Thread(target = taxonomy_lookup,
                                args = (species_name, results_queue)) for species_name in species_list]

    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
"""

def taxonomy_lookup(taxa_name):
    '''
    use NCBI taxonomy API to look up taxonomy of species_list
    '''

    Entrez.email = "finlaymaguire@gmail.com"

    esearch_handle = Entrez.esearch(db='Taxonomy', term=taxa_name)
    esearch_record = Entrez.read(esearch_handle)
    if esearch_record is None:
        print "{0} failed to get entrez ID".format(taxa_name)
        return results_queue

    entrez_species_id = esearch_record['IdList'][0]

    efetch_handle = Entrez.efetch(db='Taxonomy', id=entrez_species_id, retmode='xml')
    efetch_record = Entrez.read(efetch_handle)
    if efetch_record is None:
        print "{0} failed to get a Taxonomy".format(entrez_species_id)
        return results_queue

    entrez_species_lineage = efetch_record[0]['Lineage']

    return entrez_species_lineage

if __name__=='__main__':

    species_list_file = sys.argv[1]
    taxonomy_list_file = sys.argv[2]

    species_list = parse_species_list(sys.argv[1])

    taxa_list = []

    taxonomy_list_fh = open(taxonomy_list_file, 'w')

    for species in species_list:
        taxonomy = taxonomy_lookup(species)
        out_str = "{0}:-    {1}\n".format(species, taxonomy)
        taxonomy_list_fh.write(out_str)

    taxonomy_list_fh.close()
