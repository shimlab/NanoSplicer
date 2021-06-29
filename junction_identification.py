import pysam
import intervaltree
from intervaltree import IntervalTree
import re
import itertools
import numpy as np

# read .bam file
def get_intron_tree(pysamAlignment, chrID):
    '''
        Read bam file and generate IntervalTrees with format: 
        (intron start, intron end, num of support)
        Parameter:
            bam_filename:
                <string> filename for the bam file.
    '''
    #pysamAlignment = pysam.AlignmentFile(bam_filename)
    f_fetch = pysamAlignmentClass.fetch(chrID) # reads iterator
    f_introns = pysamAlignmentClass.find_introns(f_fetch) # dictionary

    # converty dictionary to IntervalTree
    intron_tree = IntervalTree()
    for (begin, end), data in f_introns.items():
        intron_tree.addi(begin, end, data)

    #built non overlapped range
    intron_tree_non_overlaps = intron_tree.copy()
    intron_tree_non_overlaps.merge_overlaps()

    count = 0
    single_support = 0
    total_candidate = 0
    total_junc_with_in_read = 0
    for interval in intron_tree_non_overlaps:
        candidates = find_candidate(interval.begin, interval.end, intron_tree, 10)
        
        # some statistics
        total_candidate += len(candidates)
        total_junc_with_in_read += sum([x.data for x in candidates])
        single_support += sum([x.data for x in candidates if x.data < 3])
        count += 1 
        if count < 0:
            break

# get alignment given genome pos
def find_candidate(Interval_list, window=10, min_primary = 0, 
                   min_support=0, secondary_thres=0.0, primary_thres=1.0):
    '''
    Find candidate exon boundary (i.e. intron boundary) within a given range.
    Parameter:
        begin:
            start (left-most) position of the range to be searched (0-based)
        end:
            end (right-most) possition of the range to be searched (0-based)
        tree:
            IntervalTree containing all boundary pairs 
        window: 
            window size for group surrounding boundaries (difference 
            of boundary in each size of the intron will be grouped together if 
            the absolute difference < window size)
        min_support:
            The best supported boundary need will be included only when the num
            of support reaches the minimum
        secondary_thres:
            only the junctions with multiple well supported boundary will
            be included. Well supported junction is defined as 
            secondary_thres * support num of the most supported boundary.
    '''
    # get boundaries with in searching window, sorted by the number of support
    intervals_tree = IntervalTree()
    for interval in Interval_list:
        intervals_tree.addi(interval.begin, interval.end, interval.data)    
        
    candidate_boundaries = []
    while intervals_tree:
        interval = max(intervals_tree, key = lambda x: x.data)
        best_support = interval.data
        if interval.data < min_primary: # lower bound of the support
            return candidate_boundaries
            
        #candidate_boundaries.append(interval)
        intervals_tree.remove(interval)
        
        # include surrounding boundaries
        enveloped_interval = intervals_tree.envelop(interval.begin - window, 
                                          interval.end + window)
        neighbour_found = []
        for i in enveloped_interval:
            if i.begin <= interval.begin + window and \
                    i.end >= interval.end - window:
                if i.data > secondary_thres * best_support:        
                    neighbour_found.append((interval, i))
                intervals_tree.remove(i)
        if neighbour_found:
            neighbour_found.append((interval, interval))
            count = sum([x.data for y, x in neighbour_found])
            if count >= min_support and best_support/count <= primary_thres:
                candidate_boundaries += neighbour_found
    return candidate_boundaries

def canonical_site_finder(aligned_seg, candidate_Interval, ref_FastaFile,
                              AlignmentFile, window, chrID):
    '''
        Generate candidate motif given the position of intron boundary.
        Input:
            candidate_Interval:
                tuple of intron boundary position
            ref_FastaFile:
                pysam.Fastafile class read from reference genom
            AlignmentFile:
                pysam.AlignmentFile class read from bam/sam file
            window:
                size of window in which NanoSplicer search for the candidate
            chrID: 
                Chromosome ID name.
    '''
    start, end = candidate_Interval
    
    # check determined transcript direction (by minimap2)
    try:
        ts = read.get_tag("ts")
    except:
        ts = None

    donor_pattern = ref_FastaFile.fetch(chrID, start - window, \
                                        start + window).upper()
    acceptor_pattern = ref_FastaFile.fetch(chrID, end - window, \
                                        end + window).upper()

    if ts == '+':
        intron_start_candidate = [start - window + m.start() 
                                    for m in re.finditer("GT",donor_pattern)]
        intron_end_candidate = [end - window + m.start() + 2
                    for m in re.finditer("AG",acceptor_pattern)]
        if not intron_start_candidate or not intron_end_candidate:
            return []
        else:
            return list(itertools.product(
                intron_start_candidate, intron_end_candidate))
    
    elif ts == '-':
        intron_start_candidate = [start - window + m.start() 
                                    for m in re.finditer("CT",donor_pattern)]
        intron_end_candidate = [end - window + m.start() + 2
                    for m in re.finditer("AC",acceptor_pattern)]
        if not intron_start_candidate or not intron_end_candidate:
            return []
        else:
            return list(itertools.product(
                intron_start_candidate, intron_end_candidate))
        
    else:
        # no ts tag in the bam
        intron_start_candidate1 = [start - window + m.start() 
                                    for m in re.finditer("GT",donor_pattern)]
        intron_end_candidate1 = [end - window + m.start() + 2
                    for m in re.finditer("AG",acceptor_pattern)]
        intron_start_candidate2 = [start - window + m.start() 
                                    for m in re.finditer("CT",donor_pattern)]
        intron_end_candidate2 = [end - window + m.start() + 2
                    for m in re.finditer("AC",acceptor_pattern)]

        candidate_tuples = []
        if intron_start_candidate1 and intron_end_candidate1:
            candidate_tuples += list(itertools.product(
                intron_start_candidate1, intron_end_candidate1))
        if intron_start_candidate2 and intron_end_candidate2:
            candidate_tuples += list(itertools.product(
                intron_start_candidate2, intron_end_candidate2))
        return candidate_tuples

def candidate_motif_generator(chrID, candidates_tuple, 
                              flank_size, ref_FastaFile, pattern_preference):
    '''
        generate candidate motif (genome forward strand) with flanking seq
    '''

    if not candidates_tuple:
        return None, None, None, None
    
    # candidate preference
    prefer_pattern = [1] * len(candidates_tuple)
    if pattern_preference:
        for i in range(len(candidates_tuple)):
            start, end = candidates_tuple[i]
            start_3 = ref_FastaFile.fetch(chrID, start, start + 3).upper()
            end_3 = ref_FastaFile.fetch(chrID, end -3, end).upper()
            if start_3[:2] +  end_3[-2:] != "GTAG" and \
                start_3[:2] +  end_3[-2:] != "CTAC":
                prefer_pattern[i] -= 1
                continue
            if start_3 in ['GTA', 'GTG', "CTA", 'CTG']:
                prefer_pattern[i] += 1
            if end_3 in ['CAG', 'TAG', "CAC", 'TAC']:
                prefer_pattern[i] += 1

    motif_start_pos = min([x[0] for x in candidates_tuple]) - flank_size
    motif_end_pos = max([x[1] for x in candidates_tuple]) + flank_size

    motif_list = []
    for donor, acceptor in candidates_tuple:
        motif_list.append(ref_FastaFile.fetch( 
                                    chrID, motif_start_pos, donor).upper() + 
                          ref_FastaFile.fetch( 
                                    chrID, acceptor, motif_end_pos).upper())

    return  motif_list, motif_start_pos, motif_end_pos, prefer_pattern

def find_junctions(intron_tree, window = 10):
	def split_by_begin(list_of_intervals, window = window):
		'''
		Input:
			1-D list of interval
		return:
			2-D list: splited input list by window.
		'''
		list_of_intervals = sorted(list_of_intervals, key = lambda x: x.begin)
		adj_diff = [list_of_intervals[i].begin - list_of_intervals[i-1].begin for i in range(1, len(list_of_intervals))]
		split_index = np.where(np.array(adj_diff) > window)[0] + 1
		split_by_begin = []
		pre_index = 0
		for index in np.where(np.array(adj_diff) > 10)[0] + 1:
			split_by_begin.append(list_of_intervals[pre_index:index])
			pre_index = index
		return split_by_begin
	def split_by_end(list_of_intervals, window = window):
		'''
		Input:
			1-D list of interval
		return:
			2-D list: splited input list by window.
		'''
		list_of_intervals = sorted(list_of_intervals, key = lambda x: x.end)
		adj_diff = [list_of_intervals[i].end - list_of_intervals[i-1].end for i in range(1, len(list_of_intervals))]
		split_index = np.where(np.array(adj_diff) > window)[0] + 1
		split_by_end = []
		pre_index = 0
		for index in np.where(np.array(adj_diff) > 10)[0] + 1:
			split_by_end.append(list_of_intervals[pre_index:index])
			pre_index = index
		return split_by_end
	output_list = []
	split_by_begin = split_by_begin(intron_tree)
	for l in split_by_begin:
		output_list += split_by_end(l)
	return output_list




def main():
    return None

if __name__ == "__main__":
    main()