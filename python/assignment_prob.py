import sys
import h5py
import getopt
import timeit
import os
import numpy as np
import re
import fcntl


# import NanoSplicer module and function
from . import helper
from helper import expect_squiggle_dict
from helper import parse_candidate_file
from . import junction_squiggle_selection
from .dtw import dtw


class __parse_arg(object):
    def __init__(self):
        argv = sys.argv
        if len(argv) <= 2:     
            self.__print_help()
            sys.exit(1)

        try: 
            opts, args = getopt.getopt(argv[1:], "ho:c:T:t:ab:s:",
                                       ["help=", "output_csv=",
                                        "candidate_file=", "trim_model=",
                                        "trim_signal=", "dtw_adj", 
                                        "band_prop=", "SAM="])
        except getopt.GetoptError:
            self.__print_help()
            sys.exit(1)

        args.output_file = "Untiled"
        
        try:
            self.fast5_filename = args[0]
        except:
            print("InputError: missing fast5 file name!")
            sys.exit(1)

        # DEFAULT VALUE
        self.trim_signal = 0
        self.trim_model = 4
        self.dtw_adj = False
        self.band_prop = False

        for opt, arg in opts:
            if opt in ('-h', "--help"):
                self.__print_help()
                sys.exit(0)
            elif opt in ("-o", "--output_csv"):
                self.output_file = arg
            elif opt in ("-c", "--candidate_file"):
                self.candidate_file = arg
            elif opt in ("-s", "--SAM"):
                self.sam_line = arg
            elif opt in ("-T", "--trim_model"):
                self.trim_model = int(arg)
            elif opt in ("-t", "--trim_signal"):
                self.trim_signal = int(arg)
            elif opt in ("-a", "--dtw_adj"):
                self.dtw_adj = True
            elif opt in ("-b", "--band_prop"):
                self.band_prop = float(arg)
        
    def __print_help():
        print("\n\nUsage: python {} [OPTIONS] <fast5 filename>".format(argv[0]))
        print("Options:\n\tIndexing:")
        print('\t\t-c INT\tcandidate file name (required)')
        print("\t\t-o INT\toutput csv name, default: 'Untitled'")
        print('\t\t-T \tNumber of events trimmed from scrappie model')
        print('\t\t-t \tNumber of samples trimmed from raw signal')
        return None


def main():
    os.system("mkdir -p .tmp")
    fcount_pass = ".tmp/pass"
    ftombo_fail = ".tmp/tombo_fail"
    fabnormal_sam_flag = ".tmp/abnormal_sam_flag"
    fbad_junction_mapping = ".tmp/bad_junction_mapping"

    os.system("touch {} {} {} {}".format(
        fcount_pass,
        ftombo_fail,
        fabnormal_sam_flag,
        fbad_junction_mapping
        ))
    
    def count_in_tmp(filename):
        f = open(filename, "r+")
        fcntl.flock(f, fcntl.LOCK_EX)
        count = f.read()
        count = int(count) + 1 if count else 1
        f.seek(0)
        f.write(str(count))
        f.truncate()
        f.close()

    args = __parse_arg()
    outf = open(args.output_file, 'w')

    candidates = parse_candidate_file(args.candidate_file)

    # read line in sam file
    sam_flag, mapped_pos0, cigar = [args.sam_line.strip().split('\t')[index] 
                                    for index in [1, 3, 5]]
    mapped_pos0 = int(mapped_pos0)
    # 1-based to 0-based
    mapped_pos0 -= 1

    # determine strand
    if sam_flag == "0":
        strand = "+"
    elif sam_flag == "16":
        strand = "-"
    else:
        print("Error: Abnormal mapped reads.")
        count_in_tmp(fabnormal_sam_flag)
        sys.exit(0)

    # tombo resquiggle
    try:
        tombo_results, tombo_start_clip, tombo_end_clip = \
            junction_squiggle_selection.tombo_squiggle_to_basecalls(
                args.fast5_filename)
    except:
        print("tombo resquiggle failed!!")
        count_in_tmp(ftombo_fail)
        sys.exit(0)
    
    read_length = len(tombo_results.genome_seq) + \
        tombo_start_clip + tombo_end_clip
    normalised_raw_signal = tombo_results.raw_signal/1.4826
    
    # genome pos to read pos mapping vector
    g_r_mapping = \
        junction_squiggle_selection.genome_to_read_pos_conversion(cigar)

    for candidate in candidates:
        outf.write(args.fast5_filename + ',' + strand)

        # take reverse compliment seq if nesseccary
        if strand == "-":
            candidate.sequences = [helper.reverse_complement(s)
                                   for s in candidate.sequences]

        # trim signal
        candidate.start += args.trim_signal
        candidate.end -= args.trim_signal

        # convert to read relative pos (forward direction)
        start_pos_rel_to_mapped_start = candidate.start - mapped_pos0
        end_pos_rel_to_mapped_start = candidate.end - mapped_pos0
        
        if start_pos_rel_to_mapped_start >= 0 and \
                start_pos_rel_to_mapped_start < len(g_r_mapping):
            candidate.start = g_r_mapping[start_pos_rel_to_mapped_start]
            
            # discard junction squiggle with the queried motif start/end 
            # mapped to gaps
            if candidate.start == -1:
                print("Warning: Junction squiggle start index point to mapped"
                      " intron, junction squiggle skipped.")
                count_in_tmp(fbad_junction_mapping)
                outf.write("\n")
                continue

            elif g_r_mapping[start_pos_rel_to_mapped_start] == \
                    g_r_mapping[start_pos_rel_to_mapped_start - 1]:
                candidate.start += 1 
        else:
            print("candidate start pos out of bound.")
            outf.write("\n")
            continue

        if end_pos_rel_to_mapped_start < len(g_r_mapping):
            candidate.end = g_r_mapping[end_pos_rel_to_mapped_start - 1] + 1
            
            # discard junction squiggle with the queried motif start/end mapped 
            # to gaps
            if candidate.end == -1 + 1:
                print("Warning: Junction squiggle end index point to mapped "
                      "intron, junction squiggle skipped.")
                count_in_tmp(fbad_junction_mapping)
                outf.write("\n")
                continue
        else:
            print("candidate end pos out of bound.")
#            print(candidate.start, candidate.end, mapped_pos0, g_r_mapping[-1])
            outf.write("\n")
            continue

        # get signal
        if strand == "+":
            seg_start = max(candidate.start - tombo_start_clip, 0)
            seg_end = candidate.end - tombo_start_clip
            
        elif strand == "-":
            seg_start = max(read_length - candidate.end - 1 - tombo_start_clip,
                            0)
            seg_end = read_length - candidate.start - 1 - tombo_start_clip
        
        # take into account the end clip
        seg_end = min(seg_end, len(tombo_results.segs) - 1)
        junction_squiggle = normalised_raw_signal[
            tombo_results.segs[seg_start]:tombo_results.segs[seg_end]]

        if not len(junction_squiggle):
            for i in range(len(candidate.sequences)):
                outf.write(",NA")
            outf.write("\n")
            print("read discarded")
            continue
        
        else:
            print("pass")
            count_in_tmp(fcount_pass)
        # DTW
        model_dic = expect_squiggle_dict(candidate.sequences, 
                                         trim=args.trim_model)
        junction_squiggle = junction_squiggle[abs(junction_squiggle)-3 < 0]
        for key in candidate.sequences:

            candidate_squiggle = np.array(model_dic[key], float)
            candidate_squiggle = np.array(candidate_squiggle)

            dtw_analysis = dtw(candidate_squiggle, junction_squiggle,
                               band_prop=self.band_prop,
                               dist_type="log_likelihood")
            dtw_analysis.__dtw_local_alignment()
            outf.write(',{}'.format(dtw_analysis.best_score))

        outf.write('\n')


if __name__ == "__main__":
    main()
