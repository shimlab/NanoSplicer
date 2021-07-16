'''
Subset the output of JWR_checker by:
    --genome-loc: e.g. 1-100000
    --chrID: reference name matches in the bam
    --best_JAQ: a number from 0-1, JWRs with JAQ 
                above the threshold will not be included
'''
import pandas as pd
import sys
import getopt
import textwrap
import helper

class JWR_subset_param:
    def __init__(self, arg = sys.argv):
        self.arg = arg
        try: 
            opts, args = getopt.getopt(arg[1:],"h",
                        ["help", 
                         "best_JAQ=",
                         "chrID=",
                         "genome-loc="])
        except getopt.GetoptError:
            print("ERROR:Invalid input.")
            print_help()
            sys.exit(1)
        

        # DEFAULT VALUE
        self.best_jaq, self.chrID, self.g_loc = None, None, None

        for opt, arg in opts:
            if opt in ("-h", "--help"):
                self.print_help()
                sys.exit(0)
            elif opt in ("--best_JAQ"):
                self.junc_cigar_win = float(arg)
            elif opt == "--chrID":
                self.chrID = arg
            elif opt == "--genome-loc":
                self.g_loc =\
                    tuple([int(i) for i in arg.split('-')])
                if len(self.g_loc) != 2:
                    helper.err_msg("ERROR:Invalid genome-loc.")
                    sys.exit(1)


        if len(args) <2:   
            helper.err_msg("ERROR:Invalid input.")
            self.print_help()
            sys.exit(1)

        self.input_h5, self.output_h5 = args

    def print_help(self):
        help_message =\
        '''
        Usage: python {} [OPTIONS] <hdf5 file> <output file>
        Options:
            -h/--help       Print this help text
            --bset_JAQ      A number from 0-1, JWRs with JAQ above
                             the threshold will not be included <defalt: no filter>
            --chrID         Target on specific chromosome, chrID should match
                                the chromosome name in the BAM
            --genome-loc    Target on specific genome region, chrID should be 
                                specified. e.g. --genome-loc=0-10000
        '''.format(sys.argv[0])

        print(textwrap.dedent(help_message))

def main():
    param = JWR_subset_param()
    d = pd.read_hdf(param.input_h5, 'data')
    
    if d.chrID:
        d  = d[d.chrID == param.chrID]
    
    if param.g_loc:
        d = d[d.loc.apply(lambda x: 
            x[0] >= param.g_loc[0] and x[1] <= param.g_loc[1])]
    
    if param.best_jaq:
        d = d[d.JAQ <= param.best_jaq]

    d.to_hdf(param.output_h5, 'data')
if __name__ == '__main__':
    main()