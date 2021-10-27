'''
Read junction annotation in BED format (probability merge this script into the 
junction_identification.py after testing).
'''
import pandas as pd
import sys
import os
from tqdm import tqdm
import pysam

import helper

def Mapped_junc(jwr_hdf, min_support, min_JAQ = 0):
    '''
    Input:
        jwr_hdf: HDF5 file from either JWR_checker or JWR_subset. 
        min_support: minimum number of supporting JWR required
        min_JAQ: minimum JAQ require for defining a "supporting JWR"
    Output:
        pd.DataFrame with columns:
            1. site1: 5' splice site
            2. site2: 3' splice site
        Note: the 5' and 3' is based on the direction of the reference strand,
        so we always have site1 < site2 
    '''
    
    # get all JWR
    d = pd.concat([pd.read_hdf(jwr_hdf, key = 'data'), 
                    pd.read_hdf(jwr_hdf, key = 'skipped')])
    
    d = d[d.JAQ >= min_JAQ].groupby(['loc', 'chrID'], as_index = False).count()
    d = d[d.id >= min_support]
    return(
        pd.DataFrame(
            {   'chrID': d['chrID'],
                'site1': d['loc'].apply(lambda x: x[0]),
                'site2': d['loc'].apply(lambda x: x[1])
            }
        )
    )



def GetJuncFromBed(anno_fn, GTAG_only = False, ref_FastaFile = ''):
    '''
    Read junction annotation in BED format (probability merge this script into the 
    junction_identification.py after testing).

    Read the bed line and create pd.DataFrame with following columns:
    'trans_id' : transcript ID in the BED file,
    'strand': strand recoreded in the BED file,
    'site1': splice site towards the 5' end in the forward genome strand,
    'site2': splice site towards the 3' end in the forward gemone strand
    '''
    df = pd.DataFrame(columns=['chrID','strand','site1','site2'])

    with open(anno_fn, 'r') as bed_f:
        chrID_ls, strand_ls, site1_ls, site2_ls = [], [], [], []
        for line in tqdm(bed_f):
            chrom, strand, splice_sites, *_ = ReadBedLine(line)
            for x,y in splice_sites:
                chrID_ls.append(chrom)
                strand_ls.append(strand)
                site1_ls.append(x)
                site2_ls.append(y)

    df = pd.DataFrame({'chrID' : chrID_ls,
                        'strand': strand_ls,
                        'site1': site1_ls,
                        'site2': site2_ls}).drop_duplicates()
    if GTAG_only:
        GTAG = df.apply(lambda row: check_GTAG(
                    row.chrID, row.site1, row.site2, 
                    row.strand, ref_FastaFile), axis = 1)
        df = df[~GTAG]
        
    return df



def ReadBedLine(bedline):
    """
    Take a bed12 format line as input and return the corresponding splicing
    junction position (0-based and genome related)

    return:
        name: transcript name
        strand
        splice_sites list of tuple (genome 0-based position)
        junction_pos list (transript 0-based position)
        transcript_length
    """
    bedline = bedline.strip().split()
    chrom, chromStart, chromEnd, name, score, strand, thickStart,\
     thickEnd, itemRgb, blockCount, blockSizes,\
      blockStarts, *_ = bedline
    
    blockSizes = [int(x) for x in blockSizes.strip(',').split(',')]
    blockStarts =[int(chromStart) + int(x) for x in blockStarts.strip(',').split(',')]
    blockEnds = [x + y for x,y in zip(blockSizes, blockStarts)]

    splice_sites = \
        [(blockEnds[i], blockStarts[i + 1]) for i in range(int(blockCount) - 1)]

    return chrom, strand, splice_sites


def check_GTAG(chrID, start, end, strand, ref_FastaFile):
    '''
        check the splice junction pattern
        Retrun:
            True: if GT-AG
            False: if not GT-AG
    '''
    if not ref_FastaFile:
        helper.err_msg("Error when checking GT-AG in annotation file: ref_FastaFile is required when GTAG_only=True")
        sys.exit(1)
    if not os.path.isfile(ref_FastaFile):
        helper.err_msg("Error when checking GT-AG in annotation file: {} doesn't exist!".format(ref_FastaFile))
        sys.exit(1)
    ref_FastaFile = pysam.FastaFile(ref_FastaFile)
    splice_pattern =  ref_FastaFile.fetch(
        chrID, start, start + 2).upper() + \
        ref_FastaFile.fetch(chrID, end -2, end).upper()

    if strand == '+' and splice_pattern == "GTAG":
        return True
    elif strand == '-' and splice_pattern == "CTAC":
        return True
    else:
        return False

