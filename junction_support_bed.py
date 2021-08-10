'''
Script for generate a BED file containing all the junctions identified with the number of supports
Input:
    NanoSPlicer_out.tsv

'''
import pandas as pd
import numpy as np
import sys
import getopt
import helper
import textwrap

def pd_hdf5_to_count_bed_df(f,key):
    '''
    Read data from pandas hdf file and convert to bed (junction with read count)
    return:
        count df with columns
        'chrID': chr name
        'loc' : splice junction
        'count': number of support in the hdf5 table
    '''
    d = pd.read_hdf(f, key)
    d['score'] =  np.repeat(1,len(d))
    d = d[['id', 'transcript_strand', 'chrID','loc', 'score']].groupby(
        ['chrID', 'loc','transcript_strand'], as_index=False).agg( {'id': np.random.choice, 'score':'count'})
    
    # reorganise to bed
    bed_df = pd.DataFrame({
        0:  d.chrID,
        1:  d['loc'].apply(lambda x: x[0] - 25),
        2:  d['loc'].apply(lambda x: x[1] + 25),
        3:  d.id,
        4:  np.repeat(1,len(d)),
        5:  d.transcript_strand,
        6:  d['loc'].apply(lambda x: x[0] - 25),
        7:  d['loc'].apply(lambda x: x[1] + 25),
        8:  np.repeat('0,0,255',len(d)),
        9:  np.repeat(2,len(d)),
        10: np.repeat('25,25',len(d)),
        11: ['0,{}'.format(i) for i in d['loc'].apply(lambda x: x[1]-x[0] + 25)]
    })

    return bed_df

def NanoSplicer_to_count_bed_df(f, best_p_thresh):
    '''
    Convert the bed output from NanoSplicer.py to a pd.DataFrame in python

    best_p_thresh: threshold to define confident identification
    '''
    
    df = pd.read_csv(f, sep = '\t', header = None)
    columns = df.columns
    df = df[df[4] > best_p_thresh]
    #df['count'] = np.repeat(1,len(d))
    df = df.groupby(
        [0,1,2,5,6,7,8,9,10,11],as_index=False).agg(
                                {3: np.random.choice, 4:'count'})
    
    # restore the order of columns
    df = df[columns]
    return df

def merge_count_bed_df(df1, df2):
    df = pd.concat([df1, df2])
    columns = df.columns
    
    # color columns
    df[8] = np.repeat('0,0,0', len(df))
    df = df.groupby(
        [0,1,2,5,6,7,8,9,10,11],as_index=False).agg(
                                {3: np.random.choice, 4:'sum'})

    return df[columns]


class convertor_param:
    def __init__(self, arg = sys.argv):
        self.arg = arg
        try: 
            opts, args = getopt.getopt(arg[1:],"ho:",
                        ["help", 
                         "input_BED=",
                         "input_HDF5=",
                         "best_p_thres="])
        except getopt.GetoptError:
            helper.err_msg("ERROR:Invalid input.")
            self.print_help()
            sys.exit(1)
        
        # DEFAULT VALUE
        self.input_BED,  self.input_HDF5 = None, None
        self.output = "junction_count.bed"
        self.best_p_thres = 0.8

        for opt, arg in opts:
            if opt in ("-h", "--help"):
                self.print_help()
                sys.exit(0)
            elif opt == "--input_BED":
                self.input_BED = arg
            elif opt == "--input_HDF5":
                self.input_HDF5 = arg
            elif opt == "-o":
                self.output = arg
            elif opt == "--best_p_thres":
                self.best_p_thres = arg
        
        if not self.input_BED and not self.input_HDF5:   
            helper.err_msg("ERROR:Invalid input. Please specify --input_BED and --input_HDF5")
            self.print_help()
            sys.exit(1)      

    def print_help(self):
        help_message =\
        '''
        Generate a BED file containing counts for each splice junction. For JWR didn't 
        higer than the JAQ threshold in JWR_subset (i.e. JWRs excluded in the NanoSplicer),
        the identification from initial mapping is retained. 

        Usage: python {} [OPTIONS]
        
        Options:
            -h/--help       Print this help text
            -o              Output BED filename <default: junction_count.bed>
            --input_BED     BED file from NanoSplicer.py
            --input_HDF5    HDF5 file from JWR.subset.py
            --best_p_thres  Threshold to define confident identification in NanoSplicer. 
        '''.format(sys.argv[0])
        print(textwrap.dedent(help_message))

def main():
    param = convertor_param()

    merge_count_bed_df(
        df1 = pd_hdf5_to_count_bed_df(param.input_HDF5, 'skipped'),
        df2 = NanoSplicer_to_count_bed_df(param.input_BED, 
                            best_p_thresh = param.best_p_thres)
    ).to_csv(param.output, header = None, sep = '\t')
   
   
   # pd_hdf5_to_count_bed_df(param.input_HDF5, 'skipped').to_csv(param.output, header = None, sep = '\t')

if __name__ == '__main__':
    main()
                    

