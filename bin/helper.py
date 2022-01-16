'''
collections of frequently used functions
'''

import numpy as np
import h5py
import os
from collections import defaultdict
import matplotlib.pyplot as plt
import sys

def reverse_complement(seq):
	'''
	Args: <str>
		queried seq
	Returns: <str>
		reverse_complement seq
	'''
	comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
					'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
	letters = \
		[comp[base] if base in comp.keys() else base for base in seq]
	return ''.join(letters)[::-1]

# tombo expected squiggle
def sequence_to_squiggle_tombo(seq, std_ref, customised_sd = None):
	kmer_expect = std_ref.get_exp_levels_from_seq(seq, rev_strand=False)
	if customised_sd:
		kmer_expect[1] = np.repeat(customised_sd, len(kmer_expect[1]))
	# reorganise the output into two columns which are means and sds
	kmer_expect = np.array(kmer_expect).T
	return kmer_expect

def expect_squiggle_dict_tombo(seqs, std_ref, uniform_dwell=4, customised_sd = None):
	
	if seqs:
		expect_squiggle_dic = defaultdict(list)
		for seq in seqs:
			squiggle = sequence_to_squiggle_tombo(seq = seq, std_ref = std_ref, customised_sd = customised_sd)
			if uniform_dwell > 1:			
				for mean, std in squiggle:
					expect_squiggle_dic[seq] += \
								[[mean, std]] * uniform_dwell
			else:
				expect_squiggle_dic[seq] = squiggle
	else:
		print("No valid input detected when generating expect squiggle")
		sys.exit(0)
	return expect_squiggle_dic

# scrappie squiggle
def sequence_to_squiggle(seq, trim = 0, model = 'squiggle_r94'):
	'''input:
		seq:
			<str> sequence to be converted to squiggle
		trim:
			<int> the number of events that will be ignored in each side.
		model:
			<str> scrappy model name:	{'squiggle_r94',
										'squiggle_r94_rna',
										'squiggle_r10'}
		
		output:
			numpy array: [[mean, std, dwell]...]
	'''
	simulated_seq = scrappy.sequence_to_squiggle(seq,model = model, \
		rescale =True).data(as_numpy = True, sloika = False)
	if trim:
		return simulated_seq[trim:-trim]
	else:
		return simulated_seq

def expect_squiggle_dict(seqs, trim=0, model='squiggle_r94',
						 uniform_dwell=4):
	'''
	read squiggle data from scrappie, ready for dtw
	Args:
		seqs: list of sequence motifs <list of strings>
		trim: number of bases trimmed off from both sides
				 of the candidate quiggle
		uniform_dwell: choosing between scrappie dwell time (when 
				uniform_dwell = False) and uniform dwell time ( when
				uniform_dwell = <int>, in which <int> is a integer 
				indicating the number of data points for each kmer)
	Returns: 
		python dictionary of simulated squiggle by scrappie squiggle
	Raises:

	'''
	
	if seqs:
		expect_squiggle_dic = defaultdict(list)
		for seq in seqs:
			squiggle = sequence_to_squiggle(seq = seq, trim = trim, model = model)
			for mean, std, dwell_time in squiggle:
				if uniform_dwell:
					expect_squiggle_dic[seq] += \
								[[mean, std]] * uniform_dwell
				else:
					expect_squiggle_dic[seq] += \
										[[mean, std]] * int(round(dwell_time))
	else:
		print("No valid input detected when generating expect squiggle")
		sys.exit(0)
	return expect_squiggle_dic

def parse_candidate_file(filename):
	'''
	Args:
		filename: Candidate file generated by Candidates_from_gtf.py
	Returns:
		candidate_list: list of <candidate class>
	'''
	class candidate(object):
		def __init__(self, sequences, start, end):
			self.sequences = sequences
			self.start = int(start)
			self.end = int(end)
			
	with open(filename, 'r') as f:
		candidate_list = []
		for line in f:
			line = line.strip().split(',')
			try:
				candidate_list.append(candidate(line[:-3], line[-3], line[-2]))
			except:
				print("parse_candidate_file failed")
				sys.exit(0)
	return candidate_list

def plot_dtw_alignment( long_seq, short_seq, dtw_path, dtw_score = None, \
	show_sd = True, figure_name = "Untitled", \
	figure_title = "Untitled",**plot_args):
	'''
	Args:
		figure_name: <string>
			the figure name that will be saved as
		figure_title: "string"
			figure title
		long_seq: <list or np.array>
			the long sequence in the dtw alignment
		short_seq: 2D np.array [[mean, sd]..]
			the short sequence in the dtw alignment
		dtw_path: <list or np.array>
			the best path of the alignment
		dtw_score: <INT>
			alignment score
		show_sd: True or False
			whether or not plot the sd in dash line
		**plot_args:
			arguments for matplotlib.pyplot
	Returns:
		DTW alignment visulisation.
	'''

	plt.figure(**plot_args)
	plt.plot(long_seq)
	path = np.array(dtw_path)
	plt.plot(path[:,1]-1, short_seq[[path[:,0]-1]][:,0],'g')
	if show_sd:
		plt.plot(path[:,1]-1, short_seq[[path[:,0]-1]][:,0]
						 + short_seq[[path[:,0]-1]][:,1],'g--')
		plt.plot(path[:,1]-1, short_seq[[path[:,0]-1]][:,0]
						 - short_seq[[path[:,0]-1]][:,1],'g--')

	add_info = "\nDist: {:.2f}, path length: {},  \
	Adjusted dist: {:.2f}".format(dtw_score,len(path), \
		dtw_score/len(path)) if dtw_score else ""
	
	plt.title(figure_title + add_info, fontsize=20)
	
	
	i = -1
	while os.path.exists(figure_name + ".png"):
		i += 1
	if i >-1:
		figure_name += str(i)
		
	plt.savefig(figure_name + ".png")
	plt.close()

class Fast5Class(object):
	def __init__(self, filename):
		self.filename = filename
		

	def get_read_id(self):
		with h5py.File(self.filename, 'r') as h5_f: 
			read_key = list(h5_f["Raw/Reads/"].keys())[0]
			read_id = h5_f["Raw/Reads/"][read_key].attrs['read_id']
		return read_id
		
	def get_signal(self, start = None, end = None, normalization = True,\
	 rm_outlier = False):

		start = int(start) if start else 0
		end = int(end) if end else -1

		if end <= start:
			print("InputError: Invalid start and end position when fetching \
				when fetching the raw signals")
			sys.exit(0)
		


		with h5py.File(self.filename, 'r') as h5_f: 
			read_key = list(h5_f["Raw/Reads/"].keys())[0]
			signal = list(h5_f["Raw/Reads/"][read_key]["Signal"]) 

		if normalization:
			self.get_alignment()
			signal = (signal - self.norm_shift)/self.norm_scale
		if rm_outlier:
			signal = signal[start:end]
			signal = self.remove_outlier(signal, thresh = 3)
			return(signal)

		return(signal[start:end])

	def get_alignment(self, output=None):
		with h5py.File(self.filename, 'r') as h5_f:
			path = "Analyses/"
			subpath = list(h5_f[path].keys())
			
			for i in subpath:
				if "RawGenomeCorrected" in i:
					path = path + i +'/BaseCalled_template/'

			try:
				self.mapped = dict(h5_f[path]["Alignment"].attrs)
				self.events = np.array(h5_f[path]["Events"])
				self.read_start_rel_to_raw = \
					h5_f[path]["Events"].attrs['read_start_rel_to_raw']
				self.norm_shift = h5_f[path].attrs['shift']
				 # rescale MAD to sd (https://en.wikipedia.org/wiki/Median_absolute_deviation)
				self.norm_scale = h5_f[path].attrs['scale']*1.4826
			except:
				print("Alignment doesn't exist in fast5!!")
				self.mapped = False
			
			if output == "mapping_info":
				return self.mapped
			elif output == "events":
				return self.events
			elif output == "norm_params":
				return self.norm_shift, self.norm_scale
			else:
				return None
	def remove_outlier(self, normalised_signal, thresh = 3):
		'''
		remove the data points in signal whose obsolute value reaches the thresh.

		Args:
			<list>/<np.array>Normalised signal
			<int> threshold for outliers
		Returns:
			<np.array>Normalised signal with outlier removed
		
		'''
		normalised_signal = np.array(normalised_signal)
		return normalised_signal[np.abs(normalised_signal < thresh)]	

def err_msg(msg):
	CRED = '\033[91m'
	CEND = '\033[0m'
	print(CRED + msg + CEND)	

def warning_msg(msg):
	CRED = '\033[93m'
	CEND = '\033[0m'
	print(CRED + msg + CEND)	