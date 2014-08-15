#!/usr/bin/python
import re

def read_in_files_to_lists(input_file):
	"""read files in"""
	input_file_list=[]
	input_fh=open(input_file)
	for input_file_item in input_fh:
		input_file_list.append(input_file_item)
	return input_file_list

def search_combine_annotations_and_full_details(full_accession_list,annotated_accession_list,output_filename):
	"""Searches through and combines data from full and annotated accessions"""
	output_fh=open(output_filename,"w") 
	"""annotations accession files
	>m.XXXX|ANNOTATION
	full annotations accession files
	(>m.227585)( g.227585  ORF g.227585 m.227585 )(type:)(complete)( len:)(189)( )(comp10000_c0_seq1)(:)(355-921(-))
	     1                    2                       3       4        5     6  7           8         9     10   11
	1=orf ID, 4=type of ORF 6=length 8=contig orf is from 10=location in contig 11=read direction"""
	full_acession_regex=re.compile('(>)(m\.[0-9]+)(.+)(type:)(\w+)(\s+len:)(\d+)(.+)(comp\w+)(:)(\d+\-\d+)(\(.\))')
	annotation_regex=re.compile('(>)(m\.[0-9]+)(.*)')
	for full_accession in full_accession_list:
		full_match=full_acession_regex.match(full_accession)
		for annotated_accession in annotated_accession_list:
      			annotated_match=annotation_regex.match(annotated_accession)
			if annotated_match.group(2)==full_match.group(2):
            			print (full_match.group(9)+","+full_match.group(2)+","+full_match.group(5)+","+full_match.group(7)+","+full_match.group(11)+","+full_match.group(12)+","+(annotated_match.group(3))[1:]+",\n")
            			output_fh.write(full_match.group(9)+","+full_match.group(2)+","+full_match.group(5)+","+full_match.group(7)+","+full_match.group(11)+","+full_match.group(12)+","+(annotated_match.group(3))[1:]+",\n")
	
	out_fh.close()
	
search_combine_annotations_and_full_details(read_in_files_to_lists("para_trin_uni_cds_full_accessions"),read_in_files_to_lists("para_trin_uni_best_cds_annotations_accessions"),"uni_combined_annotations_and_accessions")
search_combine_annotations_and_full_details(read_in_files_to_lists("para_trin_tet_cds_full_accessions"),read_in_files_to_lists("para_trin_tet_best_cds_annotations_accessions"),"tet_combined_annotations_and_accessions")

