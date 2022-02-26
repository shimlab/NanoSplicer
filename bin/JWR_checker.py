'''
    Finding junctions within reads (JWRs) from a bam file
    Input:
        bam_filename <str>:
            filename of the BAM
        output_name <str> <default: "JWR_full_list.hd5">:
            the filename of the output table.
    Output:
        A HDF5 with given output_name will be generated with the following 
        columns:
            1. read_id: <str> the id of original read of the JWR
            2. mapQ: <Int> the mapping quality in the input BAM
            3. transcript_strand: '+'/'-'
            4. chrID: mapped chromosome name
            5. <tuple> location of the splice junciton in the reference genome
            6. JAQ: <float> junction alignment quality
'''
import h5py
import pandas as pd
import textwrap
import pysam
import os
import re
import numpy as np
import sys
import getopt
from tqdm import tqdm
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing as mp
import helper


class JWR_to_hdf(h5py.File):
    '''
    modify the output HDF5 fil
    '''  
    def __init__(self, Filename):
        super().__init__(Filename, 'r+')
    def add_setting_info(self, name, value):
        name = 'setting/' + name
        self[name] = value

class JWR_from_reads:
    '''
    Get information from bam
    '''
    def __init__(self, read_iter):
        '''
        read_iter: read iterator from pysam. e.g. from .fetch function
        '''
        self.read_iter = read_iter
    def get_JWR(self):
        """
        Input:
            read: AlignedSegment object from pysam
        Return:
            list of introns [(start, stop),...]
            Listing the intronic sites in the reads (identified by 
            'N' in the cigar strings).
        """
        for read in self.read_iter:
            match_or_deletion = {0, 2, 7, 8} # only M/=/X (0/7/8) and D (2) are related to genome position
            ref_skip = 3
            base_position = read.pos
            for op, nt in read.cigartuples:
                if op in match_or_deletion:
                    base_position += nt
                elif op == ref_skip:
                    junc_start = base_position
                    base_position += nt
                    yield JWR_class(read, (junc_start, base_position),read.reference_name)

class JWR_class:
    def __init__(self, read, loc, chrID):
        '''
        Define a junction within read (JWR)
        Input:
            read: pysam AlignedSegment
            loc: tuple of intron location corresponding to the JWR
        '''
        self.qname = read.qname
        self.cigarstring = read.cigarstring
        self.reference_start = read.reference_start
        self.reference_end =read.reference_end
        self.loc = loc
        self.chrID = chrID
        self.mapq = read.mapping_quality
        try:
            self.ts = read.get_tag("ts")
            if read.is_reverse:
                if self.ts == '+':
                    self.ts = '-'
                elif self.ts == '-':
                    self.ts = '+'
        except:
            self.ts = None
    def get_JAQ(self, half_motif_size=25):
        '''
        Report the junction alignment quality
        '''
        junction_cigar = \
            self.get_junction_cigar(half_motif_size=25)
        junction_alignment_quality =\
            self.get_junc_map_quality(junction_cigar)
        return junction_alignment_quality
    def get_junction_cigar(self, half_motif_size):
        '''
        Return the cigar letter of around the JWR
        '''
        cigar_long = []
        for count, op in re.findall('(\d+)([A-Za-z=])', self.cigarstring):
            cigar_long += int(count) * [op]
        junc1_rel_read_start = self.loc[0] - self.reference_start
        junc2_rel_read_start = self.loc[1] - self.reference_start
        junction_cigar1 = ''
        junction_cigar2 = ''
        ref_index = -1
        for i in cigar_long:
            if i in 'MND=X':
                ref_index += 1
            if ref_index >= max(0, junc1_rel_read_start - half_motif_size) \
                    and ref_index < junc1_rel_read_start \
                    and i != 'N':
                junction_cigar1 +=  i
            if  ref_index >= junc1_rel_read_start \
                    and ref_index < min(junc2_rel_read_start + half_motif_size,
                        self.reference_end - self.reference_start) \
                    and i != 'N':
                junction_cigar2 +=  i
            if ref_index >= min(junc2_rel_read_start + half_motif_size,
                                    self.reference_end - self.reference_start):
                break   
        return junction_cigar1[-25:] + junction_cigar2[:25]
    def get_junc_map_quality(self, cigar):
        '''
        The junc map quality is simply as the proportion of 'M's within the cigar string
        '''
        if not cigar:
            return np.nan
        elif 'M' in cigar:
            return cigar.count('M')/len(cigar)
        elif '=' in cigar:
            return cigar.count('=')/len(cigar)
        else:
            return 0

class JWR_checker_param:
    def __init__(self, arg = sys.argv):
        self.arg = arg
        try: 
            opts, args = getopt.getopt(arg[1:],"hw:",
                        ["help", 
                         "window=",
                         "chrID=",
                         "genome-loc=",
                         "threads=",
                         "output_csv"])
        except getopt.GetoptError:
            helper.err_msg("ERROR:Invalid input.")
            self.print_help()
            sys.exit(1)
        
        # DEFAULT VALUE
        self.junc_cigar_win = 25
        self.chrID, self.g_loc = None, (None,None)
        self.threads = 32
        self.output_csv = False

        for opt, arg in opts:
            if opt in ("-h", "--help"):
                self.print_help()
                sys.exit(0)
            elif opt in ("-w", "--window"):
                self.junc_cigar_win = int(arg)
            elif opt == "--chrID":
                self.chrID = arg
            elif opt == "--genome-loc":
                self.g_loc =\
                    tuple([int(i) for i in arg.split('-')])
                if len(self.g_loc) != 2:
                    helper.err_msg("ERROR:Invalid genome-loc.")
                    sys.exit(1)
            elif opt == "--threads":
                self.threads = int(arg)
            elif opt == "--output_csv":
                self.output_csv = True
        
        if len(args) <2:   
            helper.err_msg("ERROR:Invalid input.")
            self.print_help()
            sys.exit(1)
        self.bamfile, self.outfile = args            

    def print_help(self):
        help_message =\
        '''
        Finding junctions within reads (JWRs) from a spliced mapping result (BAM). 

        Usage: python3 {} [OPTIONS] <BAM file> <output hdf5 file>
        
        Options:
            -h/--help        Print this help text
            -w/--window      Candidate search window size (nt) <default: 25>
            --chrID          Target a specific chromosome, chrID should match
                                the chromosome name in the BAM. All chromosomes
                                in the BAM will be searched if not specified
            --genome-loc     Target a specific genomic region, e.g. --genome-loc=0-10000
                                Use in conjunction with --chrID option. Entire
                                chromosome will be searched if location not specified
            --threads        Number of threads used <default: 32>.
            --output_csv     With this option, a csv file will be output along with the hdf5 file
        '''.format(sys.argv[0])

        print(textwrap.dedent(help_message))

def tqdm_parallel_map(executor, fn, *iterables, **kwargs):
    """
    Equivalent to executor.map(fn, *iterables),
    but displays a tqdm-based progress bar.
    
    Does not support timeout or chunksize as executor.submit is used internally
    
    **kwargs is passed to tqdm.
    """
    futures_list = []
    for iterable in iterables:
        futures_list += [executor.submit(fn, i) for i in iterable]
    for f in tqdm(concurrent.futures.as_completed(futures_list), total=len(futures_list), **kwargs):
        yield f.result()

def get_row(jwr_class_list, junc_cigar_win):
    d = pd.DataFrame(
        {'id': [jwr.qname for jwr in jwr_class_list],
             'mapQ': [jwr.mapq for jwr in jwr_class_list],
             'transcript_strand': [jwr.ts for jwr in jwr_class_list],
             'chrID': [jwr.chrID for jwr in jwr_class_list],
             'loc': [jwr.loc for jwr in jwr_class_list], 
             'JAQ': [jwr.get_JAQ(junc_cigar_win) for jwr in jwr_class_list]
            })
    return d

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def main():

    # read command line arguments
    param = JWR_checker_param()

    # check mismatch recorded in CIGAR
    algn_file = pysam.AlignmentFile(param.bamfile)
    for read in algn_file.fetch():
        if 'M' in read.cigarstring:
            warning_text =\
                '''
                    Warning:  Mismatches were not recorded in the CIGAR from the input BAM file. The JAQ
                    can still be calculated but mismatched bases will be treated as matched base. Please 
                    update the mapping setting (e.g., use '--eqx' option in minimap2) to take into account 
                    the mismatches in the JAQ calculation (recommanded).
                '''
            helper.warning_msg(textwrap.dedent(warning_text))

            break
        if '=' in read.cigarstring or 'X' in read.cigarstring:
            break

    print("Searching for JWRs ...\n\n")

    # get splice junctions from BAM
    algn_file = pysam.AlignmentFile(param.bamfile)
    reads_fetch = algn_file.fetch(param.chrID,
                                    param.g_loc[0], 
                                    param.g_loc[1])
    JWR_fetch = JWR_from_reads(reads_fetch)

    JWRs = list(JWR_fetch.get_JWR())
    # JWRs = [x for x in JWRs if x.loc[0] > param.g_loc[0] - 50 and 
    #              x.loc[1] < param.g_loc[1] + 50 ]
    
    print("Calculating JAQ for {} JWRs found".format(len(JWRs)))
    executor = concurrent.futures.ProcessPoolExecutor(param.threads)    
    futures = [executor.submit(get_row, jwr, param.junc_cigar_win) for jwr in chunks(JWRs, 500)]
    
    pbar = tqdm(total = len(JWRs))
    for future in as_completed(futures):   
        pbar.update(500)

    d = pd.concat([x.result() for x in futures])

    d.to_hdf(param.outfile, 'data')
    pd.DataFrame([]).to_hdf(param.outfile, 'skipped')
    
    # output csv file if required
    if param.output_csv:
        d.to_csv(param.outfile + ".csv")

if __name__ == "__main__":
    main()
