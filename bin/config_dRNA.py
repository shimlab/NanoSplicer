import scipy.stats

# DNA or RNA
SAMPLE_TYPE = 'RNA'

# truncated quantile
QUANTILE = 0.95

# max Z value in likelihood calculation
MAX_Z = scipy.stats.norm.ppf(QUANTILE).round(3)

# Distinguishing segment definition
DIST_SD = 0

# denoise and segmentation
SQUIGGLE_PREPROCESSING=True
WAVELET_DENOISE = False
MEDIAN_DENOISE = True
# Control the level of Wavelet denoise 
#WAVELET_LEVEL = np.ceil(np.log(len(squiggle)))
# segmentation parameter
SEGMENT_NUM = 100
MIN_BASE_OBS = 12
SEGMENT_WIN_SIZE = 8

# spike threshold
SPIKE_THRES = 4

# DTW
BANDWIDTH_PROP = 0.8
DIST_TYPE = 'manhattan'


# minimum data point assigned to a k-mer in dtw
UNIFORM_DWELL = 4
MINIMUM_POINT_FOR_DIST_SEG = 4

# prior ratio for sequence pattern
PATTERN_PREFERENCE = True
PRIOR_RATIO = 9
PRIOR_RATIO_NON_GTAG = 1

# output
SAVE_DATA = False
PLOT_ALIGNMENT = False
PLOT_SUFFIX = 'alignment'

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
