#!/usr/bin/python
import cPickle, sys, csv
import operator, bisect, numpy, math
import matplotlib
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

def import_tsv_values(tsv_to_read, dict_flag, sort_index):
    #function to read in pre-generated tsv files containing data needed to calculate average 3' and 5' end coverage
    assert dict_flag==0 or dict_flag==1
    output_array = []
    output_dict={}
    with open(tsv_to_read,'r') as tsv_fh:
        for line in csv.reader(tsv_fh, delimiter='\t'):
            if dict_flag==1:
                output_dict.update({line[0]: (line[1],line[2])})
            elif dict_flag==0:
                output_array.append(line)
    if dict_flag==1:
        return output_dict
    elif dict_flag==0:
        output_array.sort(key=operator.itemgetter(sort_index))
        return output_array

def import_pickled_values(pickled_file_to_read):
    #function to read pickled lists and dictionaries
    with open(pickled_file_to_read, 'rb') as in_file:
        return cPickle.load(in_file)
    print time.time()-start_time, "seconds"

def lookup_hits_by_contig(contig_name, contig_column_from_alignment, full_alignment_list):
    #lookup all reads that aligned with contig and output as list all read_names
    hit_list=[]
    start_location=(bisect.bisect_left(contig_column_from_alignment, contig_name))
    end_location=(bisect.bisect_right(contig_column_from_alignment, contig_name))
    for read, contig, read1_start, read2_start in full_alignment_list[start_location:end_location]:
        if contig==contig_name:
            hit_list.append((read, contig, read1_start, read2_start))
        else:
            print "fail contig: ", contig
            print "fail contig_name: ", contig_name
            assert False
    return hit_list

def coverage_array_cds_trimming(contig_name, current_contig_array, cds_data_array, just_cds_contigs):
    #trim contig coverage list to cds start and end posistions
    trimmed_coverage_array=[]
    if contig_name not in just_cds_contigs:
        print "No CDS for contig uni or tet encoding"
        return None

    cds_data_location_start=bisect.bisect_left(just_cds_contigs, contig_name)
    cds_data_location_end=bisect.bisect_right(just_cds_contigs, contig_name)
    cds_array=cds_data_array[cds_data_location_start: cds_data_location_end]

    for cds in cds_array:
        if cds[1]!=contig_name:
            assert False
        cds_co_ords_tuple=(int(cds[2])-1, int(cds[3])-1)
        print cds_co_ords_tuple
        trimmed_coverage_list=current_contig_array[cds_co_ords_tuple[0]: cds_co_ords_tuple[1]]
        trimmed_coverage_array.append(trimmed_coverage_list)

    return trimmed_coverage_array

def coverage_array_cds_trimming_with_cds_name(contig_name, current_contig_array, cds_data_array, just_cds_contigs):
    #trim contig coverage list to cds start and end posistions
    trimmed_coverage_array=[]
    if contig_name not in just_cds_contigs:
        print "No CDS for current encoding"
        return None

    cds_data_location_start=bisect.bisect_left(just_cds_contigs, contig_name)
    cds_data_location_end=bisect.bisect_right(just_cds_contigs, contig_name)
    cds_array=cds_data_array[cds_data_location_start: cds_data_location_end]

    for cds in cds_array:
        if cds[1]!=contig_name:
            assert False
        cds_co_ords_tuple=(int(cds[2])-1, int(cds[3])-1)
        trimmed_coverage_list=current_contig_array[cds_co_ords_tuple[0]: cds_co_ords_tuple[1]]
        trimmed_coverage_array.append((cds[0],trimmed_coverage_list))

    return trimmed_coverage_array

def quintiles(coverage_array_list):
    quintile_array=[]
    for coverage_array in coverage_array_list:
        quintile_pos=len(coverage_array)/5
        quintile_length=len(coverage_array[:quintile_pos])
        quint3_length=len(coverage_array[quintile_pos*2:-quintile_pos*2])
        quint1, quint2, quint3, quint4, quint5 = 0, 0, 0, 0, 0
        for pos in coverage_array[:quintile_pos]:
            quint1+=pos
        for pos in coverage_array[quintile_pos:quintile_pos*2]:
            quint2+=pos
        for pos in coverage_array[quintile_pos*2:-quintile_pos*2]:
            quint3+=pos
        for pos in coverage_array[-quintile_pos*2:-quintile_pos]:
            quint4+=pos
        for pos in coverage_array[-quintile_pos:]:
            quint5+=pos
        quintiles=[float(quint1)/quintile_length, float(quint2)/quintile_length, float(quint3)/quint3_length, float(quint4)/quintile_length, float(quint5)/quintile_length]
        quintile_array.append(quintiles)
    return quintile_array

def calculate_and_output_coverage():
    #initialise start_time to use as reference for all time counts
    start_time=time.time()

    ###Import Required Data###
    #read_name          read1_length     read2_length               (paired end reads)
    #key                    val1                  val2
    read_lengths=import_pickled_values('required_data/pickled/read_name_and_lengths_sorted_tsv.dictionary')

    #contig_name        contig_length
    #0                        1
    contig_lengths=import_pickled_values('required_data/pickled/contig_lengths_sorted.tsv_sorted_list.list')
    print "Contig lengths imported to list"
    #list_test(contig_lengths)

    #read_name          aligned_contig_name     start_pos_1_based_of_read1_alignment       start_pos_1_based_of_read2_alignment
    #val0                        val1                            val2                                        val3
    alignment_data=import_tsv_values('required_data/alignment_data_trimmed_sorted_contig.tsv', 0, 1)
    print "Alignment data imported to list"
    #list_test(alignment_data)

    just_alignment_contigs=[item[1] for item in alignment_data] #create list of just contigs for searching bisect location
    print "Contig data filtered into list from alignment data"
    #list_test(just_alignment_contigs)

    #read in files using function into arrays
    #cds_name      contig_name     cds_start_pos_on_contig      cds_end_pos_on_contig
    #0                  1                   2                           3
    universal_cds_data=import_pickled_values('required_data/pickled/uni_cds_data_sorted_list.list')
    print "Universal CDS data imported to list"
    #list_test(universal_cds_data)
    just_universal_cds_contigs=[item[1] for item in universal_cds_data]
    print "Contig data filtered into list from Universal CDS data"
    #list_test(just_universal_cds_contigs)

    tet_cds_data=import_pickled_values('required_data/pickled/tet_cds_data_sorted_list.list')
    print "Tetrahymena CDS data imported to list"
    just_tet_cds_contigs=[item[1] for item in tet_cds_data]
    print "Contig data filtered into list from Tetrahymena CDS data"
    print time.time() - start_time, "seconds\n\n"

    print "\n\n==========DATA_IMPORTED==========\n\n"
    print time.time() - start_time, "seconds"

    ###Get coverage data for CDS portions of contigs###

    # for each contig
    # set up (0,0,0,0,0,0,0,...,0,0,0,0,0,0,0,0) = length of contig

    universal_coverage_array_list=[]

    all_tet_quints, all_uni_quints= [], []
    out_universal_coverage_array_list, out_tet_coverage_array_list = [], []
    tet_coverage_array_list=[]
    uni_coverage_array_list=[]

    out_universal_coverage_array_list
    for contig_name, contig_length in contig_lengths:
        print "\n\n===CONTIG COVERAGE ARRAY GENERATION===\n\n"
        contig_length=int(contig_length)
        current_contig_array=[0]*contig_length
        #search for reads that hit contig (i.e. name)
        #array with read name, contig name, and start of both reads
        alignment_data_for_contig=lookup_hits_by_contig(contig_name, just_alignment_contigs, alignment_data)

        print "Reads aligning to ", contig_name, ": ", len(alignment_data_for_contig)

        for alignment_read_name, alignment_contig_name, read1_start, read2_start in alignment_data_for_contig:

            read_length_tuple=read_lengths.get(alignment_read_name)

            #-1 as the start pos are 1-based rather than the (0,0,0,0) contig list which is 0-based

            zero_based_read1_start=int(read1_start)-1
            zero_based_read2_start=int(read2_start)-1
            #-1 is from being 1-based rather than 0-based locations, second coord = start+read_length
            read1_co_ords=(zero_based_read1_start, zero_based_read1_start+int(read_length_tuple[0]))
            read2_co_ords=(zero_based_read2_start, zero_based_read2_start+int(read_length_tuple[1]))

            #increment those locations in current_contig_ar


            if read1_co_ords[0]>=contig_length:
                read1_co_ords=(contig_length-1, read1_co_ords[1])
            if read1_co_ords[1]>=contig_length:
                read1_co_ords=(read1_co_ords[0],contig_length-1)
            if read2_co_ords[1]>=contig_length:
                read2_co_ords=(read2_co_ords[0],contig_length-1)
            if read2_co_ords[0]>=contig_length:
                read2_co_ords=(contig_length-1,read2_co_ords[1])


            for i in range(read1_co_ords[0],read1_co_ords[1]):
                if i>=contig_length:
                    continue
                else:
                    current_contig_array[i]+=1
            for i in range(read2_co_ords[0],read2_co_ords[1]):
                if i>=contig_length:
                    continue
                else:
                    current_contig_array[i]+=1

        print "\n\n===CONTIG COVERAGE ARRAY_GENERATED===\n\n"
        print "Contig_name: ", contig_name
        print "Contig_length: " , len(current_contig_array)
        print time.time()-start_time
        universal_coverage_array_for_contig=coverage_array_cds_trimming(contig_name, current_contig_array, universal_cds_data, just_universal_cds_contigs)
        tet_coverage_array_for_contig=coverage_array_cds_trimming(contig_name, current_contig_array, tet_cds_data, just_tet_cds_contigs)
        out_universal_coverage_array_for_contig=coverage_array_cds_trimming_with_cds_name(contig_name, current_contig_array, universal_cds_data, just_universal_cds_contigs)
        out_tet_coverage_array_for_contig=coverage_array_cds_trimming_with_cds_name(contig_name, current_contig_array, tet_cds_data, just_tet_cds_contigs)

        print "\n\n===CONTIG ARRAY TRIMMED TO UNI AND TET CDS==="
        if universal_coverage_array_for_contig is not None:
            lengths=[len(x) for x in universal_coverage_array_for_contig]
            print "Contig_length_after_uni_cds_trimmed: " , lengths
        if tet_coverage_array_for_contig is not None:
            lengths=[len(x) for x in tet_coverage_array_for_contig]
            print "Contig length after tet cds trimming: " , lengths
        print time.time() - start_time, "seconds"


        if universal_coverage_array_for_contig is not None:
            for array in universal_coverage_array_for_contig:
                universal_coverage_array_list.append(array)
        if tet_coverage_array_for_contig is not None:
            for array in tet_coverage_array_for_contig:
                tet_coverage_array_list.append(array)


	if out_universal_coverage_array_for_contig is not None:
            for array in out_universal_coverage_array_for_contig:
                out_universal_coverage_array_list.append(array)
        if out_tet_coverage_array_for_contig is not None:
            for array in out_tet_coverage_array_for_contig:
                out_tet_coverage_array_list.append(array)


        print "\n\n===APPENDING TRIMMED ARRAYS TO LIST FOR UNI AND TET==="
        print "Uni coverage progress: ", len(universal_coverage_array_list), "of 10997"
        print "Tet coverage progress: ", len(tet_coverage_array_list), "of 47101"



        if universal_coverage_array_for_contig is not None:
            uni_contig_quints=quintiles(universal_coverage_array_for_contig)
            for quints in uni_contig_quints:
                all_uni_quints.append((contig_name,quints))

        if tet_coverage_array_for_contig is not None:
            tet_contig_quints=quintiles(tet_coverage_array_for_contig)
            for quints in tet_contig_quints:
                all_tet_quints.append((contig_name,quints))


    cPickle.dump(out_universal_coverage_array_list, open('coverage_calc_output/uni_pickled_coverage', 'wb'))

    with open('coverage_calc_output/uni_coverage_quints', 'a') as  uni_out:
        for array in all_uni_quints:
            uni_out.write(str(array)+'\n')

    with open('coverage_calc_output/all_uni_arrays', 'a') as uni_out2:
        for array in universal_coverage_array_list:
            uni_out2.write(str(array)+'\n')

    cPickle.dump(out_tet_coverage_array_list, open('coverage_calc_output/tet_pickled_coverage', 'wb'))

    with open('coverage_calc_output/tet_coverage_quints', 'a') as tet_out:
        for array in all_tet_quints:
            tet_out.write(str(array)+'\n')

    with open('coverage_calc_output/all_tet_arrays', 'a') as tet_out2:
        for array in tet_coverage_array_list:
            tet_out2.write(str(array)+'\n')


def import_bin_accessions(accession_file_name):
    bin_accession_list=[]
    with open(accession_file_name, "r") as accession_fh:
        for accession in accession_fh:
            bin_accession_list.append(accession.rstrip())
    return bin_accession_list

def trim_full_coverage_array_to_bin_accessions(full_coverage_array_for_encoding, bin_accessions):
    bin_trimmed_coverage_array=[]
    just_cds_names_from_full_array=[x[0] for x in full_coverage_array_for_encoding]
    for bin_cds_name in bin_accessions:
        location=bisect.bisect_left(just_cds_names_from_full_array, bin_cds_name)
        if location!=bisect.bisect_right(just_cds_names_from_full_array, bin_cds_name):
            bin_trimmed_coverage_array.append(full_coverage_array_for_encoding[location])
    return bin_trimmed_coverage_array

def calculate_posistional_average(complete_coverage_array):
    five_prime_coverage_array=numpy.array([x[1][:100] for x in complete_coverage_array if len(x[1])>=100])
    three_prime_coverage_array=numpy.array([x[1][-100:] for x in complete_coverage_array if len(x[1])>=100])
    if len(five_prime_coverage_array)!=len(three_prime_coverage_array):
        assert False
    five_prime_coverage_positional_averages=numpy.median(five_prime_coverage_array, 0)
    five_prime_coverage_posistional_std_deviations=numpy.std(five_prime_coverage_array, 0)
    three_prime_coverage_positional_averages=numpy.median(three_prime_coverage_array, 0)
    three_prime_coverage_posistional_std_deviations=numpy.std(three_prime_coverage_array, 0)
    return (five_prime_coverage_positional_averages,
            five_prime_coverage_posistional_std_deviations,
            three_prime_coverage_positional_averages,
            three_prime_coverage_posistional_std_deviations)

def accessions_to_complete_coverage(subset_universal_accessions, subset_tetrahymena_accessions, complete_cds_coverage_uni, complete_cds_coverage_tet):
    uni_accessions_for_subset=import_bin_accessions(subset_universal_accessions)
    tet_accessions_for_subset=import_bin_accessions(subset_tetrahymena_accessions)
    uni_coverage_array_for_subset=trim_full_coverage_array_to_bin_accessions(complete_cds_coverage_uni, uni_accessions_for_subset)
    tet_coverage_array_for_subset=trim_full_coverage_array_to_bin_accessions(complete_cds_coverage_tet, tet_accessions_for_subset)
    complete_coverage_array_for_subset=uni_coverage_array_for_subset+tet_coverage_array_for_subset
    return complete_coverage_array_for_subset

def plotting(figure, y_plot_data, x_plot_data, subplot_co_ords, title, y_label, x_label, labels, color_list, axis_limits):
    ax=figure.add_subplot(subplot_co_ords)
    ax.set_title(title, fontsize='10')
    ax.set_ylabel(y_label, fontsize='10')
    ax.set_xlabel(x_label, fontsize='10')


    for x_data, y_data, colors, labels in zip(x_plot_data, y_plot_data, color_list, labels):

        ax.plot(x_data, y_data, '-', color=colors, label=labels)


    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='best', prop={'size': 8})
    ax.xaxis.set_ticks(numpy.arange(-100,110,10))
    ax.yaxis.set_ticks(numpy.arange(0,130,5))
    ax.set_xticklabels(ax.get_xticks(), fontsize='8')
    ax.set_yticklabels(ax.get_yticks(), fontsize='8')
    #ax.set_yscale('log')
    ax.axis(axis_limits)
    ax.grid(True, linestyle='-', color='0.6')
    return figure


def generate_plots_of_coverage():
    universal_cds_coverage=cPickle.load(open("coverage_calc_output/uni_pickled_coverage","rb"))
    universal_cds_coverage.sort(key=operator.itemgetter(0))
    tet_cds_coverage=cPickle.load(open("coverage_calc_output/tet_pickled_coverage","rb"))
    tet_cds_coverage.sort(key=operator.itemgetter(0))

    complete_e_bin_coverage_array = accessions_to_complete_coverage("plotting/bin_accessions/phy_check_e_bin_universal_accessions",
            "plotting/bin_accessions/phy_check_e_bin_tet_accessions",
            universal_cds_coverage,
            tet_cds_coverage)
    complete_h_bin_coverage_array = accessions_to_complete_coverage("plotting/bin_accessions/phy_check_h_bin_universal_accessions",
            "plotting/bin_accessions/phy_check_h_bin_tet_accessions",
            universal_cds_coverage,
            tet_cds_coverage)
    secretome_pipeline_e_bin_coverage = accessions_to_complete_coverage("plotting/bin_accessions/Secretome_pipeline_e_bin_uni_accessions",
            "plotting/bin_accessions/Secretome_pipeline_e_bin_tet_accessions",
            universal_cds_coverage,
            tet_cds_coverage)
    secretome_pipeline_h_bin_coverage = accessions_to_complete_coverage("plotting/bin_accessions/Secretome_pipeline_h_bin_uni_accessions",
            "plotting/bin_accessions/Secretome_pipeline_h_bin_tet_accessions",
            universal_cds_coverage,
            tet_cds_coverage)
    tmhmm_1_or_more_e_bin_coverage = accessions_to_complete_coverage("plotting/bin_accessions/TMHMM_1ormorepredhel_e_bin_uni_accessions",
            "plotting/bin_accessions/TMHMM_1ormorepredhel_e_bin_tet_accessions",
            universal_cds_coverage,
            tet_cds_coverage)
    tmhmm_1_or_more_h_bin_coverage = accessions_to_complete_coverage("plotting/bin_accessions/TMHMM_1ormorepredhel_h_bin_uni_accessions",
            "plotting/bin_accessions/TMHMM_1ormorepredhel_h_bin_tet_accessions",
            universal_cds_coverage,
            tet_cds_coverage)
    tmhmm_4_or_more_e_bin_coverage = accessions_to_complete_coverage("plotting/bin_accessions/TMHMM_4ormorepredhel_e_bin_uni_accessions",
            "plotting/bin_accessions/TMHMM_4ormorepredhel_e_bin_tet_accessions",
            universal_cds_coverage,
            tet_cds_coverage)
    tmhmm_4_or_more_h_bin_coverage = accessions_to_complete_coverage("plotting/bin_accessions/TMHMM_4ormorepredhel_h_bin_uni_accessions",
            "plotting/bin_accessions/TMHMM_4ormorepredhel_h_bin_tet_accessions",
            universal_cds_coverage,
            tet_cds_coverage)
    sigp_targetp_or_wolfsortp_e_bin_coverage = accessions_to_complete_coverage("plotting/bin_accessions/SigP_TargetP_or_WolfSortP_e_bin_uni_accessions",
            "plotting/bin_accessions/SigP_TargetP_or_WolfSortP_e_bin_tet_accessions",
            universal_cds_coverage,
            tet_cds_coverage)
    sigp_targetp_or_wolfsortp_h_bin_coverage = accessions_to_complete_coverage("plotting/bin_accessions/SigP_TargetP_or_WolfSortP_h_bin_uni_accessions",
            "plotting/bin_accessions/SigP_TargetP_or_WolfSortP_h_bin_tet_accessions",
            universal_cds_coverage,
            tet_cds_coverage)


    (e_bin_five_prime_pos_av_coverages,
        e_bin_five_prime_pos_std_dev,
        e_bin_three_prime_pos_av_coverages,
        e_bin_three_prime_pos_std_dev) = calculate_posistional_average(complete_e_bin_coverage_array)

    (h_bin_five_prime_pos_av_coverage,
        h_bin_five_prime_pos_std_dev,
        h_bin_three_prime_pos_av_coverage,
        h_bin_three_prime_pos_std_dev) = calculate_posistional_average(complete_h_bin_coverage_array)

    (uni_total_five_prime_pos_av_coverage,
        uni_total_five_prime_pos_std_dev,
        uni_total_three_prime_pos_av_coverage,
        uni_total_three_prime_pos_std_dev) = calculate_posistional_average(universal_cds_coverage)

    (tet_total_five_prime_pos_av_coverage,
        tet_total_five_prime_pos_std_dev,
        tet_total_three_prime_pos_av_coverage,
        tet_total_three_prime_pos_std_dev) = calculate_posistional_average(tet_cds_coverage)

    (secretome_pipeline_e_bin_five_prime_pos_av_coverage,
        secretome_pipeline_e_bin_five_prime_pos_std_dev,
        secretome_pipeline_e_bin_three_prime_pos_av_coverage,
        secretome_pipeline_e_bin_three_prime_pos_std_dev) = calculate_posistional_average(secretome_pipeline_e_bin_coverage)

    (secretome_pipeline_h_bin_five_prime_pos_av_coverage,
        secretome_pipeline_h_bin_five_prime_pos_std_dev,
        secretome_pipeline_h_bin_three_prime_pos_av_coverage,
        secretome_pipeline_h_bin_three_prime_pos_std_dev) = calculate_posistional_average(secretome_pipeline_h_bin_coverage)

    (tmhmm_1_or_more_e_bin_five_prime_pos_av_coverage,
            tmhmm_1_or_more_e_bin_five_prime_pos_std_dev,
            tmhmm_1_or_more_e_bin_three_prime_pos_average_coverage,
            tmhmm_1_or_more_e_bin_three_prime_pos_std_dev) = calculate_posistional_average(tmhmm_1_or_more_e_bin_coverage)

    (tmhmm_1_or_more_h_bin_five_prime_pos_av_coverage,
            tmhmm_1_or_more_h_bin_five_prime_pos_std_dev,
            tmhmm_1_or_more_h_bin_three_prime_pos_av_coverage,
            tmhmm_1_or_more_h_bin_three_prime_pos_std_dev) = calculate_posistional_average(tmhmm_1_or_more_h_bin_coverage)

    (tmhmm_4_or_more_e_bin_five_prime_pos_av_coverage,
            tmhmm_4_or_more_e_bin_five_prime_pos_std_dev,
            tmhmm_4_or_more_e_bin_three_prime_pos_av_coverage,
            tmhmm_4_or_more_e_bin_three_prime_pos_std_dev) = calculate_posistional_average(tmhmm_4_or_more_e_bin_coverage)

    (tmhmm_4_or_more_h_bin_five_prime_pos_av_coverage,
            tmhmm_4_or_more_h_bin_five_prime_pos_std_dev,
            tmhmm_4_or_more_h_bin_three_prime_pos_av_coverage,
            tmhmm_4_or_more_h_bin_three_prime_pos_std_dev) = calculate_posistional_average(tmhmm_4_or_more_h_bin_coverage)

    (sigp_targetp_or_wolfsortp_e_bin_five_prime_pos_av_coverage,
            sigp_targetp_or_wolfsortp_e_bin_five_prime_pos_std_dev,
            sigp_targetp_or_wolfsortp_e_bin_three_prime_pos_av_coverage,
            sigp_targetp_or_wolfsortp_e_bin_three_prime_pos_std_dev) = calculate_posistional_average(sigp_targetp_or_wolfsortp_e_bin_coverage)

    (sigp_targetp_or_wolfsortp_h_bin_five_prime_pos_av_coverage,
            sigp_targetp_or_wolfsortp_h_bin_five_prime_pos_std_dev,
            sigp_targetp_or_wolfsortp_h_bin_three_prime_pos_av_coverage,
            sigp_targetp_or_wolfsortp_h_bin_three_prime_pos_std_dev) = calculate_posistional_average(sigp_targetp_or_wolfsortp_h_bin_coverage)


    fig=Figure(figsize=(25,30))
    canvas=FigureCanvas(fig)
    plotting(fig,
            [e_bin_five_prime_pos_av_coverages, h_bin_five_prime_pos_av_coverage, uni_total_five_prime_pos_av_coverage, tet_total_five_prime_pos_av_coverage],
            [range(1,101)]*4,
            321,
            "Average Posistional 5' CDS Coverage\n\n Complete Bins",
            "Average Number of Reads Mapping to CDS Posistion",
            "Relative posistion in CDS from 5'",
            ["Endosymbiont Bin CDS [n=3498]", "Host Bin CDS [n=27482]", "All Universal Encoding CDS [n=10997]", "All Tetrahymena Encoding CDS [n=47101]"],
            ['green','blue','magenta','cyan'],
            [0,100,0,125])


    plotting(fig,
            [e_bin_three_prime_pos_av_coverages, h_bin_three_prime_pos_av_coverage, uni_total_three_prime_pos_av_coverage, tet_total_three_prime_pos_av_coverage],
            [range(-99,1)]*4,
            322,
            "Average Posistional 3' CDS Coverage\n\n Complete Bins",
            " ",
            "Relative Posistion in CDS from 3'",
            ["Endosymbiont Bin CDS [n=3498]", "Host Bin CDS [n=27482]", "All Universal Encoding CDS [n=10997]", "All Tetrahymena Encoding CDS [n=47101]"],
            ['green','blue','magenta','cyan'],
            [-100,0,0,125])

    plotting(fig,
            [secretome_pipeline_e_bin_five_prime_pos_av_coverage, secretome_pipeline_h_bin_five_prime_pos_av_coverage,  sigp_targetp_or_wolfsortp_e_bin_five_prime_pos_av_coverage, sigp_targetp_or_wolfsortp_h_bin_five_prime_pos_av_coverage],
            [range(1,101)]*4,
            323,
            "Putative Secreted Proteins",
            "Average Number of Reads Mapping to CDS Posistion",
            "Relative Posistion in CDS from 5'",
            ["Secretome Pipeline Output Endosymbiont Bin [n=14]", "Secretome Pipeline Output Host Bin [n=221]", "SigP TargetP or WolfSortP secreted sequences Endosymbiont Bin [n=160]", "SigP TargetP or WolfSortP secreted sequences Host Bin [n=1448]"],
            ['yellow', 'red', 'purple', 'cyan'],
            [0,100,0,125])

    plotting(fig,
            [secretome_pipeline_e_bin_three_prime_pos_av_coverage, secretome_pipeline_h_bin_three_prime_pos_av_coverage, sigp_targetp_or_wolfsortp_e_bin_three_prime_pos_av_coverage, sigp_targetp_or_wolfsortp_h_bin_three_prime_pos_av_coverage],
            [range(-99,1)]*4,
            324,
            "Putative Secreted Proteins",
            " ",
            "Rleative Posistion in CDS from 3'",
            ["Secretome Pipeline Output Endosymbiont Bin [n=14]", "Secretome Pipeline Output Host Bin [n=221]", "SigP TargetP or WolfSortP secreted sequences Endosymbiont Bin [n=160]", "SigP TargetP or WolfSortP secreted sequences Host Bin [n=1448]"],
            ["yellow","red", "purple", "cyan"],
            [-100,0,0,125])

    plotting(fig,
            [tmhmm_4_or_more_e_bin_five_prime_pos_av_coverage, tmhmm_4_or_more_h_bin_five_prime_pos_av_coverage, tmhmm_1_or_more_e_bin_five_prime_pos_av_coverage, tmhmm_1_or_more_h_bin_five_prime_pos_av_coverage],
            [range(1,101)]*4,
            325,
            "Putative Transporter Proteins",
            "Average Number of Reads Mapping to CDS Posistion",
            "Relative Posistion in CDS from 5'",
            ["Sequence with 4 or more TM domain Endosymbiont Bin [n=37]", "Sequence with 4 or more TM domain Host Bin [n=2226]", "Sequences with a TM domain Endosymbiont Bin [n=310]", "Sequences with a TM domain Host Bin [n=7238]"],
            ['green', 'blue', 'red', 'black'],
            [0,100,0,125])

    plotting(fig,
            [tmhmm_4_or_more_e_bin_three_prime_pos_av_coverage, tmhmm_4_or_more_h_bin_three_prime_pos_av_coverage, tmhmm_1_or_more_e_bin_three_prime_pos_average_coverage, tmhmm_1_or_more_h_bin_three_prime_pos_av_coverage],
            [range(-99,1)]*4,
            326,
            "Putative Transport Proteins",
            " ",
            "Relative Posistion in CDS from 3'",
            ["Sequence with 4 or more TM domain Endosymbiont Bin [n=37]", "Sequence with 4 or more TM domain Host Bin [n=2226]", "Sequences with a TM domain Endosymbiont Bin [n=310]", "Sequences with a TM domain Host Bin [n=7238]"],
            ['green', 'blue', 'red', 'black'],
            [-100,0,0,125])

    fig.suptitle("Average Terminal Coverage for subsets of Paramecium-Chlorella RNA-Seq Data", fontsize=12)
    canvas.print_figure('plot_output/coverage_plot_median.svg')
    canvas.print_figure('plot_output/coverage_plot_median.png', dpi=500)
#def plotting(figure, y_plot_data, x_plot_data, subplot_co_ords, title, y_label, x_label, labels, color_list, axis_limits):

if __name__=="__main__":

    #calculate_and_output_coverage()
    generate_plots_of_coverage()

