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
                         "genome-loc=",
                         "output_csv"])
        except getopt.GetoptError:
            print("ERROR:Invalid input.")
            print_help()
            sys.exit(1)
        

        # DEFAULT VALUE
        self.best_jaq, self.chrID, self.g_loc, \
            self.output_csv = None, None, None, False

        for opt, arg in opts:
            if opt in ("-h", "--help"):
                self.print_help()
                sys.exit(0)
            elif opt in ("--best_JAQ"):
                self.best_jaq = float(arg)
            elif opt == "--chrID":
                self.chrID = arg
            elif opt == "--genome-loc":
                self.g_loc =\
                    tuple([int(i) for i in arg.split('-')])
                if len(self.g_loc) != 2:
                    helper.err_msg("ERROR:Invalid genome-loc.")
                    sys.exit(1)
            elif opt == "--output_csv":
                self.output_csv = True


        if len(args) <2:   
            helper.err_msg("ERROR:Invalid input.")
            self.print_help()
            sys.exit(1)

        self.input_h5, self.output_h5 = args

    def print_help(self):
        help_message =\
        '''
        Usage: python {} [OPTIONS] <input file: hdf5 output from JWR_checker> <output file: hdf5>
        Options:
            -h/--help       Print this help text
            --bset_JAQ      A number from 0-1, JWRs with junction alignment quality (JAQ) above
                             the threshold will not be included <defalt: no filter>
            --chrID         Target on specific chromosome, chrID should match
                                the chromosome name in the BAM
            --genome-loc    Target on specific genome region, chrID should be 
                                specified. e.g. --genome-loc=0-10000
            --output_csv    With this option, a csv file will be output with
                                 the hdf5 file
        '''.format(sys.argv[0])

        print(textwrap.dedent(help_message))

def main():
    param = JWR_subset_param()
    d = pd.read_hdf(param.input_h5, 'data')
    
    if param.chrID:
        d = d[d.chrID == param.chrID]
    
    if param.g_loc:
        d = d[d.loc.apply(lambda x: 
            x[0] >= param.g_loc[0] and x[1] <= param.g_loc[1])]
    
    if param.best_jaq:
        d = d[d.JAQ <= param.best_jaq]

    d.to_hdf(param.output_h5, 'data')
    
    if param.output_csv:
        d.to_csv(param.output_h5 + ".csv")


if __name__ == '__main__':
    main()