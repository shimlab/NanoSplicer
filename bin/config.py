import scipy.stats

# CHROMOSOME_NAME (temp)
#CHROMOSOME_NAME = "NC_000001.11"
# CHROMOSOME_NAME = "chrIS"
# IS_SEQUINS_DATA = True

# truncated quantile
QUANTILE = 0.99

# max Z value in likelihood calculation
MAX_Z = scipy.stats.norm.ppf(QUANTILE).round(3)

# Distinguishing segment definition
DIST_SD = 0

# wavelet denoising level
WAVELET_DENOISE = False

# spike threshold
SPIKE_THRES = 4

# Control the level of Wavelet denoise 
#WAVELET_LEVEL = np.ceil(np.log(len(squiggle)))

# minimum data point assigned to a k-mer in dtw
UNIFORM_DWELL = 4
MINIMUM_POINT_FOR_DIST_SEG = 4

# candidate generate
FIND_GCAG = False
FIND_ATAC = False

# prior ratio for sequence pattern
PATTERN_PREFERENCE = True
PRIOR_RATIO = 9
PRIOR_RATIO_NON_GTAG = 1

# output
SAVE_DATA = False
PLOT_ALIGNMENT = False
PLOT_LR = False
RESULT = True
OUTPUT_FILENAME = './output'

# BED output 
JWR_BED = True
BED_SIQ = -0.8
BED_COL = ['0,0,0','200,200,200']
OUTPUT_JUNC_COUNT = False
BESTQ_THRESH = 0.8
#JUNC_BED_FN = 'junction_support.bed'
