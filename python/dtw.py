import numpy as np
from numpy import exp, pi, sqrt, log
import matplotlib.pyplot as plt


class dtw(object):
    '''
        Dynamic time warping (local alignment):
            calculating the dtw distance and path between a junction squiggle and a 
            candidate squiggle. 
            Setting:
                Compare a candidate squiggle (long) to a junction squiggle (short)
                A band with given width (band_prop * long squiggle length)
            Input:
                long:
                    Candidate squiggle <np.array> with format [[mean1,sd1],...]
                short:
                    Junction squiggle <np.array (1D)>
                band_prop:
                    proportion of bandwidth relative to the long squiggle length
                dist_type:
                    Cost function
                    
        '''
    def __init__(self, candidate_squiggle, junction_squiggle,
                 band_prop, dist_type=None):
        self.candidate_squiggle = candidate_squiggle
        self.junction_squiggle = junction_squiggle
        self.band_prop = band_prop
        self.dist_type = dist_type

    def __cost(x, *y, dist_type=None):
        '''
        Define the cost function in DTW
        input:
            x,y when dist_type = "manhattan"
            x, y_mean, y_std when dist_type = "z_score" or "log_likelihood"
        '''
        def __manhattan(a, b):
            return abs(a - b)

        def __z_score(a, b_mean, b_std):
            diff = min(abs(a-b_mean))
            return diff/b_std

        def __log_likelihood(a, b_mean, b_std):
            '''
            negative log likelihood by assuming normal distribution
            retrun:
                log density in standard normal distribution
            '''
            diff = min(abs(a - b_mean))
            z = diff/b_std
            return 0.9189385 + z**2/2
        
        if len(y) not in (1,2):
            sys.stderr.write("Error: unexpected input in distance matrics.\n")
            sys.stderr.flush()
            sys.exit(1)
        
        if dist_type == None:
            dist_type = "manhattan" if len(y) == 1 else "log_likelihood" 

        if dist_type == "manhattan":
            return __manhattan(x, y[0])
        elif dist_type == "z_score":
            return __z_score(x, y[0], y[1])
        elif dist_type == "log_likelihood":
            return __log_likelihood(x, y[0], y[1])
        else:
            sys.stderr.write("Error: can't recognise the distance matrics"\
                            "['manhattan', 'z_score', 'log_likelihood']\n")
            sys.stderr.flush()
            sys.exit(1)

    #_max_sum_band_flipped
    def __dtw_local_alignment(self):

        short_len = len(self.junction_squiggle)
        long_len = len(self.candidate_squiggle)
        band_len = int(np.ceil(long_len* self.band_prop))

        # cummulated matrix for DTW
        cum_matrix = np.full((short_len + 1, long_len + 1), np.inf)
        cum_matrix[0, 0:band_len+1] = 0
        
        pre_step_matrix = np.zeros((short_len + 1, long_len + 1), dtype = int)
        '''
        matrix for recording the best path. Each of the cells contains one of three
        possible integer: 0 ,1 ,2 indecating the corresponding value in cum_matrix
        (cum_matrix[i,j]) came from the cum_matrix[i-1,j], cum_matrix[i - 1,j - 1],
        and cum_matrix[i, j - 1] respectively.
        '''
        
        for i in range(1, short_len + 1):
            for j in range(1, long_len + 1):
                
                pre_values = (cum_matrix[i-1,j], 
                            cum_matrix[i - 1,j - 1])

                pre_step_matrix[i, j] = np.random.choice(\
                                np.where(pre_values==min(pre_values))[0])

                cum_matrix[i, j] = min(pre_values) +\
                    self.__cost(self.junction_squiggle[i -1], 
                                *self.candidate_squiggle[j - 1], 
                                dist_type = self.dist_type)
        
        best_score = min(cum_matrix[-1,-band_len:])
        best_path = []
    
        traced_short_index = short_len
        traced_long_index = long_len + 1 - band_len + np.random.choice(\
                    np.where(cum_matrix[-1,-band_len:]==best_score)[0])

        while traced_short_index > 0:

            best_path.append([traced_short_index, traced_long_index])
            pre_step = pre_step_matrix[traced_short_index, traced_long_index]

            if pre_step in (0, 1):
                traced_short_index -= 1

            if pre_step in (1, 2):
                traced_long_index -= 1
        
        # best_path: 0-based coordinate on the (i+1)*(j+1) matrix
        best_path = np.array(best_path)
        
        self.best_path = best_path[::-1]
        self.best_score = best_score
        self.cum_matrix = cum_matrix
        return best_path[::-1], best_score ,cum_matrix


def main():
    short = np.array([[1,1],[1,1],[1,1],[3,1],[4,1],[5,1],[5,1],[5,1],[7,1],[5,1],[5,1],[4,1],[5,1]])
    long = np.array([1,1,1,3,4,5,7,5,4,5])
    path , score = dtw_local_alignment_max_sum(long, short)
    path , score = dtw_local_alignment_max_sum_band_flipped(long, short)
    plt.plot(long)
    plt.plot(np.array(path)[:,1]-1, short[[np.array(path)[:,0]-1]][:,0])
    plt.savefig('test.png')
    print(path)
if __name__ == '__main__':
    main()

    
            
