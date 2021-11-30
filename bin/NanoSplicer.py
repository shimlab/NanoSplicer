'''
Input:
    subset of pandas dataframe of JWR to be tested (in HDF5), using key = "data"

Output:
    Results:
        Probabilty Table
        BED file containing corrected JWRs    
    
    Visulisation (for developer only):
        Alignment figure 1:
            * Align the junction squiggle to the base candidate squiggle and 
            another candidate squiggle
            * Highlight distinguishing segment (The definition of distinguishing
            segment is defined in the config file)
    
        Alignment figure 2:
            Figure 2: log LR contribution and accumulative log LR
            * Aligned candidate squiggles same as figure 1.
            * A line short the LR contribution for each segments
            * Highlight distinguishing segment (The definition of distinguishing
            segment is defined in the config file)
    
    Candidate squiggle and junction squiggle (for developer only):
        Output CSV files for candidate squiggle and junction squiggle.

Update to the NanoSplicer.py:
    * Use JWR_checker hdf as input.
    * Default trim model changed to 0.
    * Check the failed JWR.
'''

import importlib
import textwrap
import pandas as pd
import matplotlib 
# matplotlib.rcParams.update(matplotlib.rcParamsDefault)
# output in .svg format
matplotlib.use('svg')
new_rc_params = {'text.usetex': False,
    "svg.fonttype": 'none'
    }
matplotlib.rcParams.update(new_rc_params)
import matplotlib.pyplot as plt
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing as mp
import time
import sys
import h5py
import getopt
import timeit
import os
import numpy as np
import re
import fcntl
import pysam
from tqdm import tqdm
from intervaltree import Interval
from intervaltree import IntervalTree
from pathlib import Path
from ont_fast5_api.fast5_interface import get_fast5_file        
from scipy.stats.mstats import theilslopes

# denosing
#from skimage.restoration import denoise_wavelet
import helper
from junction_identification import find_candidate, canonical_site_finder, \
                                                    candidate_motif_generator
from junc_annotation import GetJuncFromBed, Mapped_junc
from dtw import dtw
import junction_support_bed as jsb

def parse_arg():
    def print_help():
        help_message =\
        '''
        Usage: python3 {} [OPTIONS] <JWR_checker/JWR_subset hdf5 file>
        Options:
            -i      .bam/.sam file (required)
            -f      path to fast5s (required)
            -r      Genome reference file (required)
            -o      Prefix for output files <default: './output'>.
            --threads
                    Number of threads used <default: # of available cpus - 1>.

        More option on the candidate splice junctions (optional):
            User-provided splice junction (from annotation or short-read mapping):
            --junction_BED      
                    User-provided BED file containing known splice junctions. The BED 
                    file can convert from genome annotation or short read mapping result. 
                    If provided, NanoSplicer will include all splice junctions in the BED file
                    as candidates with highest preference level. e.g. --junction_BED='xx.bed'
            --non_canonical_junctions
                    User-provided BED file containing known splice junctions. The BED 
                    file be converted from genome annotation (GTF/GFF file) or short read mapping 
                    result. If provided, NanoSplicer will include non-canonical (i.e., non GT-AT) 
                    splice junctions in the BED file as candidates with highest preference level.
                    e.g. --non_canonical_junctions='xx.bed'
            Long-read mapping:
            --supported_junc_as_candidate
                    Add nearby long read supported splice junctions as candidates (if they have 
                    not been included yet).
            --min_support
                    Minimum support, i.e. number of JWRs (with minimum JAQ) mapped ,for
                     an additional splice junction to be added to candidate list <default: 3>
            --min_JAQ_support
                    Minimum JAQ (range: 0-1) for JWRs that counted as supports. <default: 0>

        For developer only:
            -T      Number of events trimmed from scrappie model <default: 0>
            -t      Number of bases trimmed from raw signam, the raw samples are
                        match to basecalled bases by tombo <default :6>
            -F      Flanking sequence size in each side of candidate searching 
                        window <default: 20>
            -w      Candidate searching window size <default: 10>
            -c      Config filename <default: 'config'>
            --plot_alignment
                    Output figures including the alignments between each junction
                    squiggle to its candidate squiggles
            --plot_LR
                    Output likelihood ratio contribution figures, which show how each
                    segment of the junction squiggle contributes to the likelihood ratio
        '''.format(argv[0])
        print(textwrap.dedent(help_message))
    
    argv = sys.argv
    
    
    try: 
        opts, args = getopt.getopt(argv[1:],"hi:f:r:o:T:t:b:F:w:c:",
                    ["help","input_alignment=","input_fast5_dir=",
                    "genome_ref=","output_path=", "trim_model=",
                    "trim_signal=","dtw_adj","bandwidth=",
                    "flank_size=", "window=", "config_file=", 
                    "annotation_bed=", "plot_alignment", "plot_LR", 
                    "non_canonical_junctions=", 'junction_BED=', 
                    'min_support=', 'supported_junc_as_candidate', 
                    'min_JAQ_support=', 'threads='])
    except getopt.GetoptError:
        helper.err_msg("Error: Invalid argument input") 
        print_help()
        sys.exit(1)


    # load config file
    config_file = 'config'
    for opt, arg in opts:
        if opt in ("-c", "--config_file"):
           config_file = str(arg)

    mdl = importlib.import_module(config_file)
    if "__all__" in mdl.__dict__:
        names = mdl.__dict__["__all__"]
    else:
        names = [x for x in mdl.__dict__ if not x.startswith("_")]
    globals().update({k: getattr(mdl, k) for k in names})

    # DEFAULT VALUE
    alignment_file, fast5_dir, genome_ref = None, None, None
    out_fn = OUTPUT_FILENAME
    plot_alignment = PLOT_ALIGNMENT
    plot_LR = PLOT_LR
    anno_fn = None
    non_cano_anno_fn = None
    trim_signal = 6
    trim_model = 0
    dtw_adj = False
    bandwidth = 0.4
    flank_size = 20
    window = 10
    include_mapped_junc = False
    min_support = 3
    min_JAQ_support = 0
    n_process = mp.cpu_count()-1

    for opt, arg in opts:
        if opt in  ("-h", "--help"):
            print_help()
            sys.exit(0)
        elif opt in ("-i", "--input_alignment"):
            alignment_file = arg
        elif opt in ("-f", "--input_fast5_dir"):
            fast5_dir = arg
        elif opt in ("-r", "--genome_ref"):
            genome_ref = arg
        elif opt in ("-o", "--output_path"):
            out_fn = arg
        elif opt in ("-T", "--trim_model"):
            trim_model = int(arg)
        elif opt in ("-t", "--trim_signal"):
            trim_signal = int(arg)
        elif opt in ("-b", "--bandwidth"):
            bandwidth = float(arg)
        elif opt in ("-F", "--flank_size"):
            flank_size = int(arg)
        elif opt in ("-w", "--window"):
            window = int(arg)
        elif opt in ("-c", "--config_file"):
            config_file = str(arg)
        elif opt == "--junction_BED":
            anno_fn = arg
        elif opt == "--plot_alignment":
            plot_alignment = True
        elif opt == "--plot_LR":
            plot_LR = True
        elif opt == "--non_canonical_junctions":
            non_cano_anno_fn = arg       
        elif opt == "--supported_junc_as_candidate":
            include_mapped_junc = True      
        elif opt == "--min_support":
            min_support = int(arg)      
        elif opt == "--min_JAQ_support":
            min_JAQ_support = int(arg)
        elif opt == "--threads":
            n_process = min(int(arg), n_process)
               

    # import global variable from config file 
    mdl = importlib.import_module(config_file)
    if "__all__" in mdl.__dict__:
        names = mdl.__dict__["__all__"]
    else:
        names = [x for x in mdl.__dict__ if not x.startswith("_")]
    globals().update({k: getattr(mdl, k) for k in names})

    if not args:
        helper.err_msg("Error: Invalid argument input")   
        print_help()# print help doc when no command line args provided
        sys.exit(0)
    pd_file = args[0]
    
    # check input
    if not alignment_file or not fast5_dir or not genome_ref:
        helper.err_msg("Error: Missing input files.") 
        sys.exit(1)        

    if include_mapped_junc and min_JAQ_support > 1:
        helper.err_msg("Error: Invalid value for '--min_JAQ_support', value range: 0-1") 
        sys.exit(1)        

    # choose the version of dtw (sum: minimize the sum of cost for the path)
    def dtw_local_alignment(candidate_squiggle, junction_squiggle, 
                            bandwidth = bandwidth, dist_type = None):
        
        return dtw(candidate_squiggle=candidate_squiggle,
                   junction_squiggle=junction_squiggle, 
                   band_prop = bandwidth,
                   dist_type = dist_type).dtw_local_alignment()
    
    return fast5_dir, out_fn, alignment_file, genome_ref, \
            bandwidth, trim_model, trim_signal, flank_size, \
             window, pd_file, anno_fn, non_cano_anno_fn, \
             plot_alignment, plot_LR, include_mapped_junc, min_support,\
             min_JAQ_support, n_process

def main():
    # parse command line argument
    fast5_dir, out_fn, alignment_file, genome_ref, bandwidth, trim_model, \
        trim_signal, flank_size, window, pd_file, anno_fn, non_cano_anno_fn, \
        plot_alignment, plot_LR, include_mapped_junc, min_support,\
        min_JAQ_support, n_process = parse_arg()

    # output filenames
    error_fn = out_fn + '_error_summary.tsv'
    prob_table_fn = out_fn + '_prob_table.tsv'
    jwr_bed_fn = out_fn + '_jwr.bed'
    support_bed_fn = out_fn + '_junc_support.bed'

    # create directory if not exit
    directory = os.path.dirname(prob_table_fn)
    Path(directory).mkdir(parents=True, exist_ok=True)

    # check file exist
    if os.path.isfile(prob_table_fn):
        helper.err_msg("File '{}' exists, please re-try by specify another output filename or remove the existing file.".format(prob_table_fn)) 
        sys.exit(1)
    if os.path.isfile(jwr_bed_fn):
        helper.err_msg("File '{}' exists, please re-try by specify another output filename or remove the existing file.".format(jwr_bed_fn)) 
        sys.exit(1)
    
    # creat empty file will header
    else:
        f = open(prob_table_fn, "w")
        f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                                            'read_id',
                                            'reference_name',
                                            'inital_junction',
                                            'JAQ',
                                            'candidates',
                                            'candidate_sequence_motif_preference',
                                            'SIQ',
                                            'prob_uniform_prior',
                                            'prob_seq_pattern_prior',
                                            'best_prob'
                                             ))
    f.close()

    # import data
    jwr_df = pd.read_hdf(pd_file, key = 'data')

    # import user-porvided bed file
    if anno_fn:
        anno_df = GetJuncFromBed(anno_fn,
                GTAG_only = False, ref_FastaFile = genome_ref)
    else:
        anno_df = pd.DataFrame([])

    if non_cano_anno_fn:
        anno_df = pd.concat([anno_df, 
                GetJuncFromBed(non_cano_anno_fn, 
                GTAG_only = True, ref_FastaFile = genome_ref)]).drop_duplicates()


    # get all mapped junction with certain hight JAQ support
    if include_mapped_junc:
        map_support_junc_df = Mapped_junc(pd_file, min_support, min_JAQ_support)
    else:
        map_support_junc_df = pd.DataFrame([])

    # get fast5 filenames recursively
    fast5_paths = Path(fast5_dir).rglob('*.fast5') # iterator
    fast5_paths = list(fast5_paths)

    # start to processing the fast5 (multiread format)
    print("running process")
    from itertools import repeat
    #executor = concurrent.futures.ProcessPoolExecutor(mp.cpu_count()-1)
    executor = concurrent.futures.ProcessPoolExecutor(n_process)
    pbar = tqdm(total = len(fast5_paths))

    futures = [executor.submit(run_multifast5, 
                                f, 
                                jwr_df, 
                                alignment_file, 
                                genome_ref, 
                                window, 
                                flank_size, 
                                trim_model, 
                                trim_signal,
                                bandwidth,
                                prob_table_fn, 
                                jwr_bed_fn,
                                anno_df,
                                map_support_junc_df,
                                plot_alignment,
                                plot_LR) for f in fast5_paths]
    
    for future in as_completed(futures):
        pbar.update(1)
    
    pd.concat([x.result() for x in futures]).to_csv(error_fn)


    # Quantification test
    if OUTPUT_JUNC_COUNT:
        jsb.merge_count_bed_df(
        df1 = jsb.pd_hdf5_to_count_bed_df(pd_file, 'skipped'),
        df2 = jsb.NanoSplicer_to_count_bed_df(jwr_bed_fn, 
                            best_p_thresh = BESTQ_THRESH)
    ).to_csv(support_bed_fn, header = None, sep = '\t')

def write_err_msg(d, jwr, error_msg):
    '''
    save failed jwr in d <pd.DataFrame> with specified error message
    '''
    d = \
        pd.concat([d, 
            pd.DataFrame([[jwr.id,
                    jwr.chrID,
                    jwr.loc, 
                    jwr.JAQ,
                    error_msg]],
        columns = ['id', 'chrID' ,'loc', 'JAQ', 'error_msg'])], ignore_index=True)
    return d

# running a single process 
def run_multifast5(fast5_path, jwr_df, AlignmentFile, ref_FastaFile, 
                    window, flank_size, trim_model, trim_signal,
                    bandwidth, prob_table_fn, jwr_bed_fn, anno_df, 
                    map_support_junc_df, plot_alignment, plot_LR):
    '''
    Description:
        Run a single process within a certain multi-read fast5
    Input:
        fast5_path:     path to the target fast5 (multi-read version only)
        jwr_df:  subset of the NanoSplcier output to be plotted(set of splice junctions to be test, it is not the 
                            splice junction candidates but the splice junctions to
                            match the junction within read to be queried)                            
        AlignmentFile:  BAM filename
        ref_fastaFile:  reference genome fasta file
        window:         candidate searching window size
        chrID:          chromosome name, come along with the all_junctions to locate the 
                            junction within reads
        flank_size:     flank size attached to the junction searching window to get candidate
                            motifs
        trim_models:    # of bases to be trimmed on candidate squiggles
        trim_signals:   # of bases to be trimmed on juntion squiggles
        bandwidth:      bandwidth for dtw
        output_file:    outbupt file name (.tsv)    
    Output:
        on each linein the output_file write NanoSplice output in the following columns

        <to be added>
    '''               
    np.random.seed(2021)

    # dtw version (function alias)
    def dtw_local_alignment(candidate_squiggle, junction_squiggle, 
                            bandwidth = bandwidth, dist_type = None):
        return dtw(candidate_squiggle=candidate_squiggle,
                   junction_squiggle=junction_squiggle, 
                   band_prop = bandwidth,
                   dist_type = dist_type, 
                   max_z = MAX_Z).dtw_local_alignment()

    # prepare the inputs
    AlignmentFile = pysam.AlignmentFile(AlignmentFile) 
    ref_FastaFile = pysam.FastaFile(ref_FastaFile)
    multi_fast5 = get_fast5_file(fast5_path, 'r')
    reads_in_file = set(multi_fast5.get_read_ids())         
    failed_jwr = pd.DataFrame(columns = ['id', 'chrID' ,'loc', 'JAQ', 'error_msg'])    

    # loop though each row in output_df
    for jwr in jwr_df.itertuples():
        if jwr.id not in reads_in_file:
            continue
        
        # search read in bam file
        overlap_read = AlignmentFile.fetch(jwr.chrID, jwr.loc[0], jwr.loc[1])
        read = None
        for i in overlap_read:
            if i.qname == jwr.id:
                read = i
                break
        if not read:
            failed_jwr = write_err_msg(failed_jwr, jwr, 'Read not found.')
            continue

        # get candidate motifs
            # find GT-AG pattern in a window nearby
        candidate_tuples = canonical_site_finder(read, jwr.loc, 
                                          ref_FastaFile, 
                                          window, jwr.chrID)

        if len(map_support_junc_df):
            # mapped splice junction in a window nearby
            temp_d = map_support_junc_df[
                    (map_support_junc_df.chrID ==  jwr.chrID) &
                    (np.abs(map_support_junc_df.site1 - jwr.loc[0]) <= window) & 
                    (np.abs(map_support_junc_df.site2 - jwr.loc[1]) <= window)
                        ]
            if len(temp_d):
                map_candidate_tuples = \
                    list(temp_d.apply(
                          lambda row: (row.site1, row.site2), axis = 1))
            else:
                map_candidate_tuples = []
            # get union of GT-AG candidate and annotation candidate
            candidate_tuples =\
                    list(set(candidate_tuples) | set(map_candidate_tuples))

        
        if len(anno_df):
            # annotation in a window nearby
            temp_d = anno_df[
                    (anno_df.chrID ==  jwr.chrID) &
                    (np.abs(anno_df.site1 - jwr.loc[0]) <= window) & 
                    (np.abs(anno_df.site2 - jwr.loc[1]) <= window)
                        ]
            if len(temp_d):
                anno_candidate_tuples = \
                    list(temp_d.apply(
                          lambda row: (row.site1, row.site2), axis = 1))
            else:
                anno_candidate_tuples = []
            # get union of GT-AG candidate and annotation candidate
            candidate_tuples =\
                    list(set(candidate_tuples) | set(anno_candidate_tuples))
        else:
            anno_candidate_tuples = []

        if jwr.loc not in candidate_tuples:
            candidate_tuples.append(jwr.loc)

        candidate_motif, motif_start_ref, motif_end_ref, \
            candidate_preference = candidate_motif_generator(
                                        jwr.chrID, candidate_tuples, 
                                        flank_size, ref_FastaFile,
                                        pattern_preference = PATTERN_PREFERENCE)

        # adjust annotated junction (if provided)
        if len(anno_candidate_tuples):
            for junc in anno_candidate_tuples:
                candidate_preference[candidate_tuples.index(junc)] = 3
        
        candidate_preference = np.array(candidate_preference)
        
        if not candidate_motif:# or len(candidate_motif) == 1:  
            failed_jwr = write_err_msg(failed_jwr, jwr, 'Fail to find candidates.')
            continue
        index_m = candidate_tuples.index(jwr.loc)

        # check reverse
        candidate_motif_rev = [helper.reverse_complement(seq) 
                                for seq in candidate_motif]
        if read.is_reverse:
            candidates_for_read = candidate_motif_rev
        else:
            candidates_for_read = candidate_motif     
 
        # trim signal (use fewer base)
        motif_start_ref += trim_signal
        motif_end_ref -= trim_signal

        #tombo resquiggle
        try:
            tombo_results, tombo_start_clip, tombo_end_clip, std_ref, sd_of_median= \
                tombo_squiggle_to_basecalls(multi_fast5, read)
        except:
            failed_jwr = write_err_msg(failed_jwr, jwr, 'Fail to resquiggle (Tombo).')
            continue
        
        read_length = len(tombo_results.genome_seq) \
                            + tombo_start_clip + tombo_end_clip

        # read coordinate of JWR start and end
        jwr_loc_read = locate_jwr_bases(read_length, std_ref, sd_of_median, 
                            read, motif_start_ref, motif_end_ref)
        
        if not jwr_loc_read:
            failed_jwr = write_err_msg(failed_jwr, jwr, 'Fail to get the JWR bases in read.')
            continue

     
        # get junction squiggle
        junction_squiggle = get_junc_squiggle(tombo_results, 
                                tombo_start_clip, read_length, 
                                read, jwr_loc_read, spike_thres=SPIKE_THRES)
        
        if not len(junction_squiggle):
            failed_jwr = write_err_msg(failed_jwr, jwr, 
                                'Fail to get junction squiggle')
            continue

        # candidate squiggle
        model_dic = helper.expect_squiggle_dict_tombo(seqs=candidates_for_read, 
                                                std_ref = std_ref,
                                                uniform_dwell=UNIFORM_DWELL,
                                                customised_sd = None)

        cum_path = {}
        score_dict = {}
        squiggle_match = {}
        output_prefix = ''
        num_of_cand = len(candidates_for_read)
        aligned_base = {}

        for j, candidate in enumerate(candidates_for_read):
            candidate_squiggle = np.array(model_dic[candidate],float)
            
            if WAVELET_DENOISE:
                wavelet_levels = int(np.ceil(np.log2(len(junction_squiggle))))
                junction_squiggle_wav = \
                    denoise_wavelet(junction_squiggle, method='BayesShrink', 
                              mode='soft', wavelet_levels=wavelet_levels,
                               wavelet='sym8', rescale_sigma='True')
                path , score, cum_matrix = \
                    dtw_local_alignment(candidate_squiggle = candidate_squiggle, 
                                        junction_squiggle = junction_squiggle_wav)
            else:
                path , score, cum_matrix = \
                    dtw_local_alignment(candidate_squiggle = candidate_squiggle, 
                                    junction_squiggle = junction_squiggle)

            # log likelihood (from DTW score to log likelihood)
            score_dict[j] = -1 * score
            
            # candidate squiggle in matched len of junction suqiggle
            squiggle_match[j] = candidate_squiggle[path[:,1] - 1, :]
            
            # change to sd for segment median
            squiggle_match[j][:,1] = [sd_of_median] * squiggle_match[j].shape[0]
            aligned_base[j] = candidate[trim_model + int((path[0,1] - 1)/UNIFORM_DWELL): trim_model+ int((path[-1,1] - 1)/UNIFORM_DWELL) + 1]
            #cum_path[j] = cum_matrix[path[:, 0], path[:, 1]]

            # Alignment plot
            if plot_alignment and j == num_of_cand - 1:
                preference_text_dict = {
                    0: "non GT-AG",
                    1: "GT(-)-AG(-)",
                    2: "GT(+)-AG(-) or GT(-)-AG(+)",
                    3: "GT(+)-AG(+) or annotated"
                }
                # minimap2 candidate as reference
                matched_candidate_ref = squiggle_match[index_m]
                squiggle_match_list = [squiggle_match[x] for x in range(num_of_cand)]
                junction_squiggle_median_list = \
                    [median_denoise(junction_squiggle, [matched_candidate_ref, x]) for x in squiggle_match_list]
                segment_list = \
                    [get_distinguishing_segment(matched_candidate_ref, [x]) for x in squiggle_match_list]

                # LR calculation (distinguishing segment only)
                dist_seg_LR = []

                for cand in range(num_of_cand):
                    junction_squiggle_median = junction_squiggle_median_list[cand]
                    segment, is_dist_seg = segment_list[cand]
                    is_dist_seg = is_dist_seg & (segment[:,1] - segment[:,0] >= MINIMUM_POINT_FOR_DIST_SEG)
                    p_wise_LR = __log_likelihood(junction_squiggle_median, 
                                    squiggle_match[cand][:,0],squiggle_match[cand][:,1],max_z = MAX_Z) -\
                                 __log_likelihood(junction_squiggle_median, 
                                    matched_candidate_ref[:,0],matched_candidate_ref[:,1], max_z = MAX_Z)

                    dist_seg_LR.append(sum(p_wise_LR[segment[is_dist_seg][:,0]]))
                
                
                rel_to_ref_LR = np.exp(dist_seg_LR)/np.exp(dist_seg_LR[index_m])
                post_prob = rel_to_ref_LR/sum(rel_to_ref_LR)
                
                # add sequence pattern prior
                post_prob_prior = post_prob.copy()
                post_prob_prior[candidate_preference==0] =\
                    post_prob[candidate_preference==0]*PRIOR_RATIO_NON_GTAG
                post_prob_prior[candidate_preference==2] =\
                    post_prob[candidate_preference==2]*PRIOR_RATIO
                post_prob_prior[candidate_preference==3] =\
                    post_prob[candidate_preference==3]*(PRIOR_RATIO**2)
                post_prob_prior = post_prob_prior/sum(post_prob_prior)
                
                fig = plt.figure(figsize=(30,5 * num_of_cand), dpi = 400)
                subplot_index = 0
            
                for cand in range(num_of_cand):
                    junction_squiggle_median = junction_squiggle_median_list[cand]
                    segment, is_dist_seg = segment_list[cand]
                    is_dist_seg = is_dist_seg & (segment[:,1] - segment[:,0] >= MINIMUM_POINT_FOR_DIST_SEG)

                    if cand == index_m:
                        continue
                    
                    # frame subplot
                    subplot_index += 1
                    ax = fig.add_subplot(num_of_cand -1,1,subplot_index)
                    #ax.margin
                    for i, seg in enumerate(segment):
                        ax.axvspan(seg[0]-0.5, seg[1]-1+0.5, -0.2,1.2, 
                            clip_on=True, alpha = 0.2 if is_dist_seg[i] else 0.1 , edgecolor = 'black', 
                            lw = 1,  facecolor = 'red' if is_dist_seg[i] else '#E2D3CA')
                    # # data-point wise LR

                    # LR_con = __log_likelihood(junction_squiggle_median, 
                    #             matched_candidate_ref[:,0],matched_candidate_ref[:,1]) - \
                    #         __log_likelihood(junction_squiggle_median, 
                    #             squiggle_match[cand][:,0],squiggle_match[cand][:,1])

                    #ax.plot(distinguish_points_idx,
                    #                np.repeat(max(junction_squiggle_median)+0.5,len(distinguish_points_idx)), 'o')
                    ax.plot(matched_candidate_ref[:,0], '#ee9b00',linewidth=3, alpha = 0.8,
                            label = "Candidate squiggle matching the initial mapping")
                    ax.plot(matched_candidate_ref[:,0] + matched_candidate_ref[:,1], 'y--',linewidth=0.8, 
                            label = "Position of +/- 1 standard deviation")
                    ax.plot(matched_candidate_ref[:,0] - matched_candidate_ref[:,1], 'y--',linewidth=0.8)
                    ax.plot(squiggle_match[cand][:,0], '#0a9396',linewidth=3, alpha = 0.5,
                            label = "Candidate squiggle matching the true splice junction")
                    ax.plot(squiggle_match[cand][:,0] + squiggle_match[cand][:,1], 'y--',linewidth=0.8)
                    ax.plot(squiggle_match[cand][:,0] - squiggle_match[cand][:,1], 'y--',linewidth=0.8)
                    
                    mid_seg_pos = (segment[:,1] + segment[:,0]-1) /2
                    ax.plot(mid_seg_pos,junction_squiggle_median[np.array([int(x) for x in mid_seg_pos])], 'r.',ms=8, 
                            label = "Median of current measurements in segments")
                    ax.plot(junction_squiggle, '#0096c7',linewidth=1, 
                            label = "Junction squiggle")
                    
                    middle_pos = find_candidate_middle_pos(squiggle_match[index_m])
                    for i in range(len(aligned_base[index_m])):
                        ax.annotate(aligned_base[index_m][i],
                            xy=(middle_pos[i], max(junction_squiggle) + 0.05), xycoords='data', ha='center', color = '#ee9b00')
                    middle_pos = find_candidate_middle_pos(squiggle_match[cand])
                    for i in range(len(aligned_base[cand])):
                        ax.annotate(aligned_base[cand][i],
                            xy=(middle_pos[i], min(junction_squiggle) - 0.2), xycoords='data', ha='center', color = '#0a9396')
                    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize = 10)
                    
                    #figure title
                    if True:
                        ax.set_title('Log LR: {:.2f}, average Log-Lik: {:.2f}, post probability: {:.3f}, candidate_preference: {}, post probability(seq prior ratio = 9): {:.3f} '.format(
                            dist_seg_LR[cand], 
                            normalise_loglik(
                                np.mean(
                                even_wise_log_likelihood(junction_squiggle, squiggle_match_list[cand], max_z=MAX_Z)),
                                sd = sd_of_median,
                                C = SPIKE_THRES, 
                                qunatile = QUANTILE
                            ),
                            post_prob[cand], 
                            preference_text_dict[candidate_preference[cand]], 
                            post_prob_prior[cand]), y=1.0, pad=-14)

                fig.suptitle('minimap2 candidate likelihood: {:.2f}, average Log-Lik: {:.2f}, post probability: {:.3f}, candidate_preference: {}, post probability(seq prior ratio = 9): {:.3f}'.format(
                        dist_seg_LR[index_m], 
                        normalise_loglik(
                            np.mean(even_wise_log_likelihood(junction_squiggle, squiggle_match_list[index_m], max_z=MAX_Z)),
                            sd = sd_of_median,
                            C = SPIKE_THRES,
                            qunatile = QUANTILE
                            ),
                        post_prob[index_m],
                        preference_text_dict[candidate_preference[index_m]], 
                        post_prob_prior[index_m]))

                #fig.subplots_adjust(right=0.8)

                fig_dir = os.path.dirname(prob_table_fn)
                fig.savefig(
                    os.path.join(
                        fig_dir, "{}_median_pairwise.png".format(read.qname)), 
                                format = 'png')

                # fig.savefig(
                #     os.path.join(
                #         fig_dir, "{}_median_pairwise.png".format(read.qname)))

                plt.close()
            # LR contribution plot
            if plot_LR and j == num_of_cand - 1:

                # minimap2 candidate as reference
                matched_candidate_ref = squiggle_match[index_m]
                squiggle_match_list = [squiggle_match[x] for x in range(num_of_cand)]
                junction_squiggle_median_list = \
                    [median_denoise(junction_squiggle, [matched_candidate_ref, x]) for x in squiggle_match_list]
                segment_list = \
                    [get_distinguishing_segment(matched_candidate_ref, [x]) for x in squiggle_match_list]

                # LR calculation (distinguishing segment only)
                dist_seg_LR = []
                p_wise_LR_list = []
                
                for cand in range(num_of_cand):
                    junction_squiggle_median = junction_squiggle_median_list[cand]
                    segment, is_dist_seg = segment_list[cand]
                    is_dist_seg = is_dist_seg & (segment[:,1] - segment[:,0] >= MINIMUM_POINT_FOR_DIST_SEG)
                    p_wise_LR = __log_likelihood(junction_squiggle_median, 
                                    squiggle_match[cand][:,0],squiggle_match[cand][:,1],max_z = MAX_Z) -\
                                 __log_likelihood(junction_squiggle_median, 
                                    matched_candidate_ref[:,0],matched_candidate_ref[:,1], max_z = MAX_Z)
                    p_wise_LR_list.append(p_wise_LR)

                    dist_seg_LR.append(sum(p_wise_LR[segment[is_dist_seg][:,0]]))
                
                rel_to_ref_LR = np.exp(dist_seg_LR)/np.exp(dist_seg_LR[index_m])
                post_prob = rel_to_ref_LR/sum(rel_to_ref_LR)
                post_prob_prior = post_prob.copy()
                post_prob_prior[candidate_preference==0] =\
                    post_prob[candidate_preference==0]*PRIOR_RATIO_NON_GTAG
                post_prob_prior[candidate_preference==2] =\
                    post_prob[candidate_preference==2]*PRIOR_RATIO
                post_prob_prior[candidate_preference==3] =\
                    post_prob[candidate_preference==3]*(PRIOR_RATIO**2)
                post_prob_prior = post_prob_prior/sum(post_prob_prior)
                fig = plt.figure(figsize=(20,5 * num_of_cand))
                subplot_index = 0
                for cand in range(num_of_cand):
                    junction_squiggle_median = junction_squiggle_median_list[cand]
                    segment, is_dist_seg = segment_list[cand]
                    is_dist_seg = is_dist_seg & (segment[:,1] - segment[:,0] >= MINIMUM_POINT_FOR_DIST_SEG)

                    if cand == index_m:
                        continue
                    
                    # frame subplot
                    subplot_index += 1
                    ax = fig.add_subplot(num_of_cand -1,1,subplot_index)
                    for i, seg in enumerate(segment):
                        ax.axvspan(seg[0], seg[1]-1, -0.2,1.2, 
                            clip_on=True, alpha = 0.2 if is_dist_seg[i] else 0.1 , edgecolor = 'black', 
                            lw = 0.2,  facecolor = 'red' if is_dist_seg[i] else '#E2D3CA')

                    # # data-point wise LR

                    # LR_con = __log_likelihood(junction_squiggle_median, 
                    #             matched_candidate_ref[:,0],matched_candidate_ref[:,1]) - \
                    #         __log_likelihood(junction_squiggle_median, 
                    #             squiggle_match[cand][:,0],squiggle_match[cand][:,1])

                    #ax.plot(distinguish_points_idx,
                    #                np.repeat(max(junction_squiggle_median)+0.5,len(distinguish_points_idx)), 'o')
                    ax.plot(matched_candidate_ref[:,0], 'r',linewidth=0.6, 
                            label = "minimap2 candidate")
                    ax.plot(matched_candidate_ref[:,0] + matched_candidate_ref[:,1], 'y--',linewidth=0.4, 
                            label = "sd")
                    ax.plot(matched_candidate_ref[:,0] - matched_candidate_ref[:,1], 'y--',linewidth=0.4)
                    ax.plot(squiggle_match[cand][:,0], 'g',linewidth=0.6,
                            label = "Candidate")

                    if True:
                        ax2 = ax.twinx()
                        ax2.plot(p_wise_LR_list[cand], 
                                color='#4b0082', linewidth=0.4, alpha=1, 
                                dash_capstyle='round', label = 'LR contribution')
                        ax2.set_ylabel('log likelihood ratio contribution', color='#4b0082')
                        ax2.axhline(0, linestyle = '--', lw = 0.6)
                    
                    middle_pos = find_candidate_middle_pos(squiggle_match[index_m])
                    for i in range(len(aligned_base[index_m])):
                        ax.annotate(aligned_base[index_m][i],
                            xy=(middle_pos[i], max(junction_squiggle) + 0.05), xycoords='data', ha='center')
                    middle_pos = find_candidate_middle_pos(squiggle_match[cand])
                    for i in range(len(aligned_base[cand])):
                        ax.annotate(aligned_base[cand][i],
                            xy=(middle_pos[i], min(junction_squiggle) - 0.2), xycoords='data', ha='center')
                    ax.legend(frameon=False, loc='best')
                    
                    ax.set_title('Log LR: {:.2f}, post probability: {:.3f}, candidate_preference: {}, post probability(seq prior ratio = 9): {:.3f} '.format(
                        dist_seg_LR[cand], post_prob[cand], candidate_preference[cand], post_prob_prior[cand]), y=1.0, pad=-14)
                fig.suptitle('minimap2 candidate likelihood: {:.2f}, post probability: {:.3f}, candidate_preference: {}, post probability(seq prior ratio = 9): {:.3f}'.format(
                        dist_seg_LR[index_m], post_prob[index_m],candidate_preference[index_m], post_prob_prior[index_m]))

                fig_dir = os.path.dirname(prob_table_fn)
                fig.savefig(
                    os.path.join(fig_dir, "{}_LR_median_pairwise.png".format(
                                read.qname)))

                plt.close()
                #plot cumulative contribution of LR of each data point
                if False:
                    ax2.plot(cum_path[1] - cum_path[0], 
                            color='#4b0082', linewidth=0.2, alpha=1, dash_capstyle='round')
                    #plot diff of candidate instead
                    #plt.plot(squiggle_match[1][:,1] - squiggle_match[0][:,1])
                    plt.savefig("{}_cum_LR_{}.png".format(output_prefix, read.qname))
                    plt.close()
                
                    # med_j, mad_j = medmad(junction_squiggle)
                    # junction_squiggle_norm = (junction_squiggle - med_j)/ mad_j
                    # candidate_squiggle[:,0] = (candidate_squiggle[:,0] - med_c) / mad_c
                    # candidate_squiggle[:,1] = candidate_squiggle[:,1] / mad_c
                
            # Plot DTW 
            if False:
                print(score)
                fig, axes = plt.subplots(nrows=4, figsize=(30,40))
                axes[0].plot(junction_squiggle,linewidth = 10)
                axes[0].tick_params(labelsize=40)
                            #axes[1].figure(figsize=(10,5))
                axes[1].plot(candidate_squiggle[:,0], color = "orange",linewidth =7)
                axes[1].plot(candidate_squiggle[:,0] + candidate_squiggle[:,1], ":",color = "orange",linewidth = 4)
                axes[1].plot(candidate_squiggle[:,0] - candidate_squiggle[:,1], ":",color = "orange",linewidth = 4)
                axes[1].tick_params(labelsize=40)   
                axes[1].set_title("Scrappie model",fontsize=30, pad = 1)
                axes[2].plot(junction_squiggle,linewidth = 10)
                axes[2].tick_params(labelsize=40)
                #path = np.array(path[::-1])
                axes[2].plot(path[:,0]-1, candidate_squiggle[[path[:,1]-1]][:,0],'',color = "orange",linewidth = 7)
                axes[2].plot(path[:,0]-1, candidate_squiggle[[path[:,1]-1]][:,0]\
                                                + candidate_squiggle[[path[:,1]-1]][:,1],':',color = "orange",linewidth = 4)
                axes[2].plot(path[:,0]-1, candidate_squiggle[[path[:,1]-1]][:,0]\
                                                - candidate_squiggle[[path[:,1]-1]][:,1],':',color = "orange",linewidth = 4)
                # axes[2].set_title(
                #     str(i) + \
                #     "\nDist: {:.2f}, path length: {}, Adjusted dist:"\
                #     " {:.2f}, prior: {}".format(
                #             score,len(path),
                #             score/len(path), 
                #             [minimap_priors[i], NanoSplicer_priors[i]][j])
                #         ,fontsize=40, pad = 1)

                pos = axes[3].imshow(cum_matrix, cmap='hot', interpolation='nearest',aspect='auto')
                fig.colorbar(pos,ax = axes[3])
                axes[3].plot(path[:,1], path[:,0])
                axes[3].set_title("Alignment path",fontsize=30, pad = 20)
                if read.is_reverse:
                    strand = 1
                else:
                    strand = 0
                

                fig_dir = os.path.dirname(prob_table_fn)
                fig.savefig(
                    os.path.join(fig_dir,"{}_fig{}_{}.png".format(output_prefix,read.qname, j)))
                plt.close()

            # NanoSplicer result
            if RESULT and j == num_of_cand - 1:
                # minimap2 candidate as reference
                matched_candidate_ref = squiggle_match[index_m]
                squiggle_match_list = [squiggle_match[x] for x in range(num_of_cand)]
                

                junction_squiggle_median_list = \
                    [median_denoise(junction_squiggle, [matched_candidate_ref, x]) for x in squiggle_match_list]
                
                try:
                    segment_list = \
                        [get_distinguishing_segment(matched_candidate_ref, [x]) for x in squiggle_match_list]
                except:
                    failed_jwr = write_err_msg(failed_jwr, jwr, 'Fail to get segments')
                    continue
                # LR calculation (distinguishing segment only)
                dist_seg_logLR = []
                dist_seg_logLR_individual = []

                for cand in range(num_of_cand):
                    junction_squiggle_median = junction_squiggle_median_list[cand]
                    segment, is_dist_seg = segment_list[cand]

                    # remove distinguishing segment with low number of points
                    is_dist_seg = is_dist_seg & (segment[:,1] - segment[:,0] >= MINIMUM_POINT_FOR_DIST_SEG)

                    # point wise likelihood ratio (match the junction squiggle len)
                    p_wise_logLR = __log_likelihood(junction_squiggle_median, 
                                    squiggle_match[cand][:,0],squiggle_match[cand][:,1], max_z = MAX_Z) -\
                                 __log_likelihood(junction_squiggle_median, 
                                    matched_candidate_ref[:,0],matched_candidate_ref[:,1], max_z = MAX_Z)

                    # total log LR (a segment is one data point)
                    dist_seg_logLR.append(sum(p_wise_logLR[segment[is_dist_seg][:,0]]))

                    # list of log LRs for each distinguishing segment
                    dist_seg_logLR_individual.append(p_wise_logLR[segment[is_dist_seg][:,0]])

                # reverse the log and calculate the posterio
                dist_seg_logLR = np.array(dist_seg_logLR)
                rel_to_ref_LR = np.exp(dist_seg_logLR - dist_seg_logLR[index_m])
                post_prob = rel_to_ref_LR/sum(rel_to_ref_LR)
                post_prob_prior = post_prob.copy()
                post_prob_prior[candidate_preference==0] =\
                    post_prob[candidate_preference==0]*PRIOR_RATIO_NON_GTAG
                post_prob_prior[candidate_preference==2] =\
                    post_prob[candidate_preference==2]*PRIOR_RATIO
                post_prob_prior[candidate_preference==3] =\
                    post_prob[candidate_preference==3]*(PRIOR_RATIO**2)
                post_prob_prior = post_prob_prior/sum(post_prob_prior)
                # segment abs diff
                even_logL_list = [even_wise_log_likelihood(junction_squiggle, x, max_z=MAX_Z) for x in squiggle_match_list]
                p_wise_Si = [score_dict[x]/len(junction_squiggle) for x in range(num_of_cand)]
                n_of_aligned_event = [len(x) for x in even_logL_list]
                segment_Si = [normalise_loglik(np.mean(x),sd = sd_of_median, C = SPIKE_THRES, qunatile = QUANTILE) for x in even_logL_list]
                #worst_even_logL_list = [sorted(x)[0] for x in even_logL_list]
                f = open(prob_table_fn, "a")
                #fcntl.flock(f,fcntl.LOCK_EX)
                # f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                #     jwr.id,
                #     str(jwr.loc),
                #     jwr.JAQ,
                #     ','.join([str(x)+','+str(y) for x,y in candidate_tuples]),
                #     index_m,
                #     ','.join([str(x) for x in candidate_preference]),
                #     np.max(segment_Si),
                #     #','.join([str(x) for x in segment_Si]),
                #     ','.join([str(x) for x in dist_seg_logLR]),
                #     ','.join([str(x) for x in post_prob]),
                #     ','.join([str(x) for x in post_prob_prior]),
                #     np.max(post_prob_prior)
                #     ))
                #     # fcntl.flock(f,fcntl.LOCK_UN)
                # f.close()
                f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    jwr.id,
                    jwr.chrID,
                    str(jwr.loc),
                    '{:.3f}'.format(jwr.JAQ),
                    ','.join(['('+str(x)+','+str(y)+')' for x,y in candidate_tuples]),
                    ','.join([str(x) for x in candidate_preference]),
                    '{:.4f}'.format(np.max(segment_Si)),
                    ','.join(['{:.4f}'.format(x) for x in post_prob]),
                    ','.join(['{:.4f}'.format(x) for x in post_prob_prior]),
                    '{:.4f}'.format(np.max(post_prob_prior))
                    ))
                    # fcntl.flock(f,fcntl.LOCK_UN)
                f.close()

                if JWR_BED and np.max(segment_Si) > BED_SIQ:
                    junc_start, junc_end =\
                         candidate_tuples[np.argmax(post_prob_prior)]

                    try:
                        ts = read.get_tag("ts")
                        if read.is_reverse:
                            if ts == '+':
                                ts = '-'
                            elif ts == '-':
                                ts = '+'
                    except:
                        ts = None
                

                    f_bed = open(jwr_bed_fn, "a")
                    
                    best_q = np.max(post_prob_prior)
                    bed_col =\
                         BED_COL[0] if best_q > BESTQ_THRESH else BED_COL[1]
                    
                    f_bed.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        jwr.chrID,
                        str(junc_start - 25),
                        str(junc_end + 25),
                        jwr.id,
                        '{:.4f}'.format(best_q),
                        ts,
                        '0',
                        '0',
                        bed_col,
                        '2',
                        '25,25',
                        '0,{}'.format(junc_end-junc_start + 25),
                        ))
                        # fcntl.flock(f,fcntl.LOCK_UN)
                    f_bed.close()
            
            # save squiggles as csv
            if SAVE_DATA:
                np.savetxt("candidate_{}_{}({}).csv".format(j, jwr.id,jwr.loc[0]), squiggle_match[j], delimiter=",")
                np.savetxt("squiggle_{}({}).csv".format(read.qname,jwr.loc[0]), junction_squiggle, delimiter=",")

    return failed_jwr

def get_gaps_in_read(AlignedSegment):
    blocks = AlignedSegment.get_blocks()
    gap = set([(blocks[i-1][1], blocks[i][0]) for i in range(len(blocks))])
    return gap

def tombo_squiggle_to_basecalls(multi_fast5, AlignedSegment):
    '''
    script from tombo (https://gist.github.com/marcus1487/5cf083644143aafc770dceb213783346)
    Input:
        read_fast5: fast5 filename
    Returns:
        rsqgl_results <class>: Resquiggle result
           # useful rsqgl_results attributes:
                #   - genome_seq: clipped basecalls to model-able positions (can't get k-mers and thus expected levels for first and last few bases)
                #   - raw_signal: median normalized raw signal values (only for portion of raw signal mapped to model-able basecalls)
                #   - read_start_rel_to_raw: start of rsqgl_results.raw_signal and rsqgl_results.segs within original all_raw_signal
                #   - segs: start position for each model-able base (should be one longer than length of genome_seq to include last end point) 
    '''
    from tombo import tombo_helper, tombo_stats, resquiggle
    # extract read info
    #fast5s_f = get_fast5_file(read_fast5_fn, 'r')
    #fast5_data = h5py.File
    seq_samp_type = tombo_helper.seqSampleType('RNA', True)
    #seq_data = resquiggle.get_read_seq(fast5_data, 'Basecall_1D_000', 'BaseCalled_template', seq_samp_type, 0)
    
    if AlignedSegment.is_reverse:
        read_seq = helper.reverse_complement(AlignedSegment.seq)
    else:
        read_seq = AlignedSegment.seq

    if AlignedSegment.query_qualities:
        mean_q_score = AlignedSegment.query_qualities
    else:
        mean_q_score = None
    
    seq_data = tombo_helper.sequenceData(seq=read_seq, 
                                         id=AlignedSegment.qname, 
                                         mean_q_score=mean_q_score)
    # prep tombo objects
    std_ref = tombo_stats.TomboModel(seq_samp_type=seq_samp_type)
    start_clip = std_ref.central_pos
    end_clip = std_ref.kmer_width - std_ref.central_pos - 1
    rsqgl_params = tombo_stats.load_resquiggle_parameters(seq_samp_type)

    # extract raw signal
    #all_raw_signal = tombo_helper.get_raw_read_slot(fast5_data)['Signal'][:]
    all_raw_signal = multi_fast5.get_read(AlignedSegment.qname).get_raw_data()
   
    # if this is a direct RNA read, flip raw signal to process from 5' to 3'
    
    if seq_samp_type.rev_sig:
        all_raw_signal = all_raw_signal[::-1]

    # spoof mapping read to its own basecalls (with pseudo-clipped bases in genome_location for model-able postions)
    map_results = tombo_helper.resquiggleResults(
        align_info=tombo_helper.alignInfo(ID=seq_data.id, Subgroup='BaseCalled_template', ClipStart=0, ClipEnd=0, Insertions=0,
                                        Deletions=0, Matches=len(seq_data.seq), Mismatches=0),
        genome_loc=tombo_helper.genomeLocation(Start=std_ref.central_pos, Strand='+', Chrom='NA'),
        genome_seq=seq_data.seq, raw_signal=all_raw_signal, mean_q_score=seq_data.mean_q_score)

    # align raw signal to basecalls
    try:
        rsqgl_results = resquiggle.resquiggle_read(
        map_results, std_ref, rsqgl_params, all_raw_signal=all_raw_signal, outlier_thresh = np.inf)
    except:
        rsqgl_params = rsqgl_params._replace(bandwidth=2000)
        rsqgl_results = resquiggle.resquiggle_read(
        map_results, std_ref, rsqgl_params, all_raw_signal=all_raw_signal, outlier_thresh = np.inf)

    # estimation sd of segment median
    sd_of_median = sd_from_tombo(rsqgl_results,std_ref,read_seq)

    return rsqgl_results, start_clip, end_clip, std_ref, sd_of_median

def median_denoise(junction_squiggle, squiggle_match_list):
    median_junction_squiggle = np.array([])
    event = junction_squiggle[0:1]
    
    for i in range(1, len(junction_squiggle)):
        if all([(x[i] == x[i-1]).all() for x in squiggle_match_list]):
            event = np.append(event, junction_squiggle[i])
        else:
            median_junction_squiggle = np.append(median_junction_squiggle,[np.median(event)] * len(event))
            event = np.array([junction_squiggle[i]])
    median_junction_squiggle = np.append(median_junction_squiggle,[np.median(event)] * len(event))
    return(median_junction_squiggle)

def even_wise_log_likelihood(junction_squiggle, squiggle_match, max_z):
    even_wise_likelihood = np.array([])
    event = junction_squiggle[0:1]
    
    for i in range(1, len(junction_squiggle)):
        if (squiggle_match[i] == squiggle_match[i-1]).all():
            event = np.append(event, junction_squiggle[i])
        else:
            even_wise_likelihood = np.append(
                even_wise_likelihood,
                __log_likelihood(
                    np.median(event), 
                    squiggle_match[i-1,0], 
                    squiggle_match[i-1,1], max_z = max_z))
            
            event = np.array([junction_squiggle[i]])
    even_wise_likelihood = np.append(
        even_wise_likelihood,
        __log_likelihood(
                    np.median(event), 
                    squiggle_match[-1,0], 
                    squiggle_match[-1,1], max_z = max_z)
        )
    return(even_wise_likelihood)

def sd_from_tombo(tombo_results, std_ref, read_seq):
    '''
    Input: 
        tombo_results：tombo resquiggle results
        std_ref: tombo model
        read_seq: strand sensitive read sequence
    '''
    med_seg = np.array([np.median(
        tombo_results.raw_signal[tombo_results.segs[i]:tombo_results.segs[i+1]])
        for i in range(len(tombo_results.segs) - 1)])
    expected_means, _ = std_ref.get_exp_levels_from_seq(read_seq, rev_strand=False)
    return np.sqrt(np.sum((expected_means - med_seg)**2)/(len(expected_means) - 1))

def find_candidate_middle_pos(squiggle_match):
    '''
    function for plot bases: find the middel point of each aligned k-mer
    '''
    out = [0]
    for i in range(1, len(squiggle_match)):
        if (squiggle_match[i] == squiggle_match[i-1]).all():
            out[-1] += 0.5
        else:
            out.append(i)
    return(out)

def get_distinguishing_segment(matched_candidate_ref, 
                matched_candidate_ls):
    '''
    Input:
        matched_candidate_ref: For the candidate that used as reference,  mean and sd from 
                                                                that matched to junction squiggle
        matched_candidate_ls: a list of queried candidate squiggle
    Output:
        list of distinguishing segment: [(seg1_start, seg1_end), ...]
    '''
    matched_candidate_ls.append(matched_candidate_ref)

    all_means = np.array([x[:,0] for x in matched_candidate_ls])
    is_end_boundary = np.any(all_means[:,1:] != all_means[:,:-1], axis = 0)
    end_index = np.arange(1,len(is_end_boundary) + 1)[is_end_boundary]
    segment = list(zip(np.append(0, end_index[:-1]), end_index))
    seg_start = np.array(segment)[:,0]
    
    is_dist_seg = \
            np.any(
                np.abs(
                    all_means[:,seg_start] - matched_candidate_ref[seg_start,0]) \
                    > DIST_SD * matched_candidate_ref[seg_start,1], axis = 0)
    
    return np.array(segment), is_dist_seg
        
def get_junc_map_quality(cigar):
    '''
    The junc map quality is simply as the proportion of 'M's within the cigar string
    '''
    if not cigar:
        return np.nan
    else:
        return cigar.count('M')/len(cigar)

def __log_likelihood(a, b_mean, b_std, 
                        max_z):
    '''
     log likelihood by assuming normal distribution
    retrun:
        log density in standard normal distribution
    '''
    def norm_log_density(x, mean, sd, max_diff = None):                
        diff = np.abs(mean - x)
        if max_diff is not None:
            diff = np.minimum(diff, max_diff)
        z = diff/sd
        return -np.log(sd)-0.9189385 - z**2/2 #norm

    if max_z:
        max_diff = max_z * b_std
        return norm_log_density(a, b_mean, b_std, max_diff)
    else:
        return norm_log_density(a, b_mean, b_std) 

def locate_jwr_bases(read_length, std_ref, sd_of_median, 
                            read, motif_start_ref, motif_end_ref):
    '''
    locating the JWR bases within a read:
    read: AlignmentSegment class from pysam

    return:
        0-based read coordinate of JWR start and end
    '''

    # genome pos to read pos mapping vector
    g_r_mapping = genome_to_read_pos_conversion(read.cigarstring)

    # convert to read relative pos (forward direction)
    start_pos_rel_to_mapped_start = motif_start_ref - read.reference_start
    end_pos_rel_to_mapped_start = motif_end_ref - read.reference_start

    # make sure the read covers the genome coordinate
    if  start_pos_rel_to_mapped_start >= 0 \
                and start_pos_rel_to_mapped_start < len(g_r_mapping):
        motif_start_read = g_r_mapping[start_pos_rel_to_mapped_start]

        if motif_start_read == -1: #queried motif start/end mapped to intron 'N'
            return None
        elif g_r_mapping[start_pos_rel_to_mapped_start] == \
            g_r_mapping[start_pos_rel_to_mapped_start - 1]:
            motif_start_read += 1 
    else:
        return None
    
    if end_pos_rel_to_mapped_start < len(g_r_mapping):
        motif_end_read = g_r_mapping[end_pos_rel_to_mapped_start - 1] + 1
        if motif_end_read == -1 + 1:
            return None
        else:
            return motif_start_read, motif_end_read
    else:
        return None

def get_junc_squiggle(tombo_results, tombo_start_clip, read_length, 
                        read, jwr_loc_read, spike_thres):
    '''
    Return junction squiggle with spike removed
    '''
    motif_start_read, motif_end_read = jwr_loc_read

    if not read.is_reverse:
        seg_start = max(motif_start_read - tombo_start_clip, 0)
        seg_end = motif_end_read - tombo_start_clip
    else:
        seg_start = \
            max(read_length - motif_end_read -1 - tombo_start_clip, 0)
        seg_end = read_length - motif_start_read - 1 - tombo_start_clip
    
    # take into account the end clip
    seg_end = min(seg_end, len(tombo_results.segs) - 1)
    
    # getting junction squiggle
    signal = tombo_results.raw_signal[
        tombo_results.segs[seg_start]:tombo_results.segs[seg_end]]

    # check output
    if not len(signal):
        return []
    else:
        return np.array(signal, float)
    


    # spike removal
    junction_squiggle = junction_squiggle[abs(junction_squiggle) < spike_thres]   

def genome_to_read_pos_conversion(cigar):
    '''
    Input:
        g_pos: 0-base position within genome
        cigar: CIGAR string from BAM/SAM file
    
    Returns:
        r_pos: 0-base position within a read, '-1' if the corresponding genome 
            coordinate is skipped ('N' in cigar) in the alignment 
    '''
    cigar_long = []
    for count, type in re.findall('(\d+)([A-Za-z=])', cigar):
        cigar_long += int(count) * [type]
    r_index = -1
    g_r_mapping = []
    for i in cigar_long:
        if i in "MIS=X":
            r_index += 1
        if i in "MD=X":
            g_r_mapping.append(r_index)
        if i == "N":
            g_r_mapping.append(-1)
    return g_r_mapping

def normalise_loglik(loglike, sd, C, qunatile):
    quantile_value = scipy.stats.norm.ppf(qunatile)*sd
    density = scipy.stats.norm.pdf(quantile_value, scale = sd)
    
    # area under the modified normal distribution
    auc_tail = (C-quantile_value)*density*2
    auc_body = abs(0.5-qunatile) * 2
    return loglike - np.log(auc_tail + auc_body)
if __name__ == "__main__":
    main()