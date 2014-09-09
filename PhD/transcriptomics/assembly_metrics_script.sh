#!/bin/zsh

# Script to optionally run TrinityStats.pl script on each assembly in folder and combine outputs

ASSEMBLY_FOLDER=$(pwd)
TRINITY_STATS=/storage/fin/assemblers/trinityrnaseq_r20140717/util/TrinityStats.pl
SUMMARY_FILE=${ASSEMBLY_FOLDER}/assembly_metric_summary.csv
cd $ASSEMBLY_FOLDER

usage() 
{
echo "Usage: $0 [-s re_run stats on all assemblies] [-c compile results of stats]" 1>&2; exit 1;
}

STATS_FLAG=0
COMPILE_FLAG=0


while getopts "sc" opt; do
        case "$opt" in
                s) 
		  STATS_FLAG='1'
		  ;;
		c)
		  COMPILE_FLAG='1'
		  ;;
                \?)
		  usage
		  ;;
        esac
done


if [ $STATS_FLAG -eq 1 ]; then
	for assembly in $(ls -d */ | sed 's/\///'); do
		$TRINITY_STATS ${assembly}/Trinity.fasta > ${assembly}/${assembly}_assembly_stats &
	done
	echo "All stats processes forked: waiting for completion"
	wait
	echo "All stats processes completed"	
fi

if [ $COMPILE_FLAG -eq 1 ]; then
	echo "assembly,gc_percentage,transcripts_number,transcript_n50,median_transcript_length,average_transcript_length,transcript_assembled_bases,gene_number,gene_n50,median_gene_length,average_gene_length,gene_assembled_bases" > $SUMMARY_FILE
	for assembly_stats in $(ls -d */*assembly_stats); do
		name=$(echo "$assembly_stats" | sed 's/\/.*//')
		gc_percentage=$(grep "Percent GC:" ${assembly_stats} | grep -oP "[0-9]+.[0-9]+" | sed 's/\n//')
		transcript_number=$(grep "Total trinity transcripts:" ${assembly_stats} | grep -oP "[0-9]+" | sed 's/\n//')
		gene_number=$(grep "Total trinity 'genes':" ${assembly_stats} | grep -oP "[0-9]+" | sed 's/\n//' )
	
		transcript_stats=$(grep -A12 "Stats based on ALL transcript contigs:" ${assembly_stats})
		transcript_n50=$(echo "$transcript_stats" | grep "Contig N50:" | grep -oP " [0-9]+" | sed 's/ //' | sed 's/\n//' )
		median_transcript_length=$(echo "$transcript_stats" | grep "Median contig length:" | grep -oP "[0-9]+" | sed 's/\n//')
		average_transcript_length=$(echo "$transcript_stats" | grep "Average contig:" | grep -oP "[0-9]+.[0-9]+" | sed 's/\n//' )
		transcript_assembled_bases=$(echo "$transcript_stats" | grep "Total assembled bases:" | grep -oP "[0-9]+" | sed 's/\n//')
		
		gene_stats=$(grep -A12 "Stats based on ONLY LONGEST ISOFORM per 'GENE':" ${assembly_stats})
		gene_n50=$(echo "$gene_stats" | grep "Contig N50:" | grep -oP " [0-9]+" | sed 's/ //' | sed 's/\n//')
		median_gene_length=$(echo "$gene_stats" | grep "Median contig length:" | grep -oP "[0-9]+" | sed 's/\n//')
		average_gene_length=$(echo "$gene_stats" | grep "Average contig:" | grep -oP "[0-9]+.[0-9]+" | sed 's/\n//')
		gene_assembled_bases=$(echo "$gene_stats" | grep "Total assembled bases:" | grep -oP "[0-9]+" | sed 's/\n//')
		
		echo "${name},${gc_percentage},${transcript_number},${transcript_n50},${median_transcript_length},${average_transcript_length},${transcript_assembled_bases},${gene_number},${gene_n50},${median_gene_length},${average_gene_length},${gene_assembled_bases}" >> $SUMMARY_FILE
	done	
fi

