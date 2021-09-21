'''
Read junction annotation in BED format (probability merge this script into the 
junction_identification.py after testing).
'''
import pandas as pd
import sys
from tqdm import tqdm
def GetAnnoFromBed(anno_fn):
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
                        'site2': site2_ls})
    return df.drop_duplicates()



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