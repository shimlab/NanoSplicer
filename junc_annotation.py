'''
Read junction annotation in BED format (probability merge this script into the 
junction_identification.py after testing).
'''
import pandas as pd
import sys

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
        for line in bed_f:
            chrom, strand, splice_sites, *_ = ReadBedLine(line)
            for x,y in splice_sites:
                df = df.append(pd.Series({
                    'chrID' : chrom,
                    'strand': strand,
                    'site1': x,
                    'site2': y}), ignore_index = True)
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

    junction_pos = []
    for start_pos in blockStarts:
        exon_transcript_pos = start_pos - int(chromStart)
        for intron_pos in splice_sites:
            if start_pos >=  intron_pos[1]:
                exon_transcript_pos -= intron_pos[1] - intron_pos[0]
        junction_pos.append(exon_transcript_pos)
    
    junction_pos = junction_pos[1:]

    transcript_length = sum(blockSizes)

    return chrom, strand, splice_sites, junction_pos, transcript_length