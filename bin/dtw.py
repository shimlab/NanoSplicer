import numpy as np
from numpy import exp, pi, sqrt, log
import matplotlib.pyplot as plt
import sys

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
                 band_prop=1, dist_type=None, max_z = None):
        self.candidate_squiggle = candidate_squiggle
        self.junction_squiggle = junction_squiggle
        self.band_prop = band_prop
        self.dist_type = dist_type
        self.max_z = max_z

    def _cost(self, x, *y, dist_type=None):
        from scipy.stats import norm
        '''
        Define the cost function in DTW
        input:
            x,y when dist_type = "manhattan"
            x, y_mean, y_std when dist_type = "z_score" or "log_likelihood"
        '''
        def _manhattan(a, b):
            return abs(a - b)

        def _z_score(a, b_mean, b_std):
            diff = min(abs(a-b_mean))
            return diff/b_std


        # normal version
        def _log_likelihood(a, b_mean, b_std, 
                                max_z = self.max_z):
            '''
            negative log likelihood by assuming normal distribution
            retrun:
                log density in standard normal distribution
            '''
            def norm_log_density(x, mean, sd, max_diff = None):
                diff = np.abs(mean - x)
                if max_diff:
                    diff = np.minimum(diff, max_diff)
                z = diff/sd
                return -np.log(sd)-0.9189385 - z**2/2 #norm

            if max_z:
                max_diff = max_z * b_std
                return -1 * norm_log_density(a, b_mean, b_std, max_diff)
            else:
                return -1 * norm_log_density(a, b_mean, b_std)           
            # diff = abs(a - b_mean)
            # z = diff/b_std
            # laplacc_b = b_std/np.sqrt(2)
            # #return 0.9189385 + z**2/2 #norm
            # #return 1.14473 + log(1+z**2) # t with df = 1
            # return np.log(2*laplacc_b) + diff/laplacc_b


        if len(y) not in (1,2):
            sys.stderr.write("Error: unexpected input in distance matrics.\n")
            sys.stderr.flush()
            sys.exit(1)
        
        if dist_type == None:
            dist_type = "manhattan" if len(y) == 1 else "log_likelihood" 

        if dist_type == "manhattan":
            return _manhattan(x, y[0])
        elif dist_type == "z_score":
            return _z_score(x, y[0], y[1])
        elif dist_type == "log_likelihood":
            return _log_likelihood(x, y[0], y[1])
        else:
            sys.stderr.write("Error: can't recognise the distance matrics"\
                            "['manhattan', 'z_score', 'log_likelihood']\n")
            sys.stderr.flush()
            sys.exit(1)
   
    def dtw_alignment(self):
        pass

class dtw_local_alignment(dtw):
    #_max_sum_band_flipped
    def dtw_alignment(self, **kwargs):
        # seed for tie-breaker when trace back the DTW path
        np.random.seed(seed=2021)

        short_len = len(self.junction_squiggle)
        long_len = len(self.candidate_squiggle)
        band_len = int(np.ceil(long_len* self.band_prop))
        band_move = (long_len - band_len)/short_len

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
        #kmer_bound = [False, False] + \
        #             [not np.array_equal(long[i],long[i-1])
        #              for i in range(1, len(long))]
        
        for i in range(1, short_len + 1):
            j_start = int((i - 1) * band_move)
            for j in range(j_start, j_start + band_len):
                
                pre_values = (cum_matrix[i-1,j], 
                            cum_matrix[i - 1,j - 1],)
                            #cum_matrix[i, j - 1],

                pre_step_matrix[i, j] = np.random.choice(\
                                np.where(pre_values==min(pre_values))[0])

                cum_matrix[i, j] = min(pre_values) +\
                    self._cost(self.junction_squiggle[i -1], 
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
        
        return dtw_res(best_score,best_path[::-1],cum_matrix, 
                       self.candidate_squiggle,
                        self.junction_squiggle,
                        self.band_prop,
                        self.max_z)

class dtw_global_alignment_seg_mean(dtw):
    
    def dtw_alignment(self, min_dwell = 1):
        def _update_mean(x , pre_mu, pre_n):
            return (pre_mu*pre_n+x)/(pre_n+1),pre_n+1
        # seed for tie-breaker when trace back the DTW path
        np.random.seed(seed=2021)

        short_len = len(self.junction_squiggle)
        long_len = len(self.candidate_squiggle)

        # cummulated matrix for DTW
        cum_matrix = np.full((short_len + 1, long_len + 1,3), np.inf)
        # init the count (for mean calculation) to 1
        cum_matrix[:,:,1] = 0
        cum_matrix[:,:,2] = 1
        cum_matrix[0, 0, 0] = 0
        cum_matrix[0, 0, 2] = np.inf
        pre_step_matrix = np.zeros((short_len + 1, long_len + 1), dtype = int)
        '''
        matrix for recording the best path. Each of the cells contains one of three
        possible integer: 0 ,1 ,2 indecating the corresponding value in cum_matrix
        (cum_matrix[i,j]) came from the cum_matrix[i-1,j], cum_matrix[i - 1,j - 1],
        and cum_matrix[i, j - 1] respectively.
        '''
        #kmer_bound = [False, False] + \
        #             [not np.array_equal(long[i],long[i-1])
        #              for i in range(1, len(long))]
        
        for i in range(1, short_len + 1):
            #j_start = int((i - 1) * band_move)
            #for j in range(j_start, j_start + band_len):
            for j in range(1, long_len + 1):    
                # pre_values = (cum_matrix[i-1,j,:], 
                #             cum_matrix[i - 1,j - 1,:])
                            #cum_matrix[i, j - 1],
                current_cost = self._cost(self.junction_squiggle[i-1], 
                                *self.candidate_squiggle[j-1], 
                                dist_type = self.dist_type)

                if cum_matrix[i-1,j-1,2] < min_dwell:
                    new_mu, new_n = _update_mean(current_cost,cum_matrix[i-1,j,1], cum_matrix[i-1,j,2])        
                    
                    pre_step_matrix[i, j] = 0
                    cum_matrix[i, j, 0] = cum_matrix[i-1, j, 0]   
                    cum_matrix[i, j, 1] = new_mu
                    cum_matrix[i, j, 2] =  new_n     
                    
                else:
                    new_mu, new_n = _update_mean(current_cost,cum_matrix[i-1,j,1], cum_matrix[i-1,j,2])
                    new_values = [new_mu+cum_matrix[i-1,j,0],
                                cum_matrix[i-1,j-1,0] + cum_matrix[i-1,j-1,1] + current_cost]
                    
                    
                    best = np.random.choice(\
                                    np.where(new_values==min(new_values))[0])
                    
                    pre_step_matrix[i, j] = best
                    cum_matrix[i, j, 0] = cum_matrix[i-1, j, 0] if best == 0 else cum_matrix[i-1,j-1,0] + cum_matrix[i-1,j-1,1]         
                    cum_matrix[i, j, 1] = new_mu if best == 0 else current_cost     
                    cum_matrix[i, j, 2] =  new_n if best == 0 else 1               
        
        #best_score = min(cum_matrix[-1,-band_len:])
        best_score = cum_matrix[-1,-1, 0] + cum_matrix[-1,-1, 1]
        best_path = []
    
        traced_short_index = short_len
        # traced_long_index = long_len + 1 - band_len + np.random.choice(\
        #             np.where(cum_matrix[-1,-band_len:]==best_score)[0])
        traced_long_index = long_len

        while traced_short_index > 0:

            best_path.append([traced_short_index, traced_long_index])
            pre_step = pre_step_matrix[traced_short_index, traced_long_index]

            if pre_step in (0, 1):
                traced_short_index -= 1

            if pre_step in (1, 2):
                traced_long_index -= 1
        
        # best_path: 0-based coordinate on the (i+1)*(j+1) matrix
        best_path = np.array(best_path)
        res_matrix = cum_matrix[:,:,0] + cum_matrix[:,:,1]
        return dtw_res(best_score,best_path[::-1],res_matrix, 
                       self.candidate_squiggle,
                        self.junction_squiggle,
                        self.band_prop,
                        self.max_z)

class dtw_local_alignment_seg_mean(dtw):
    
    def dtw_alignment(self,min_dwell = 1):
        def _update_mean(x , pre_mu, pre_n):
            return (pre_mu*pre_n+x)/(pre_n+1),pre_n+1
        # seed for tie-breaker when trace back the DTW path
        np.random.seed(seed=2021)

        short_len = len(self.junction_squiggle)
        long_len = len(self.candidate_squiggle)
        band_len = int(np.ceil(long_len* self.band_prop))
        band_move = (long_len - band_len)/short_len
        # cummulated matrix for DTW
        cum_matrix = np.full((short_len + 1, long_len + 1,3), np.inf)
        # init the count (for mean calculation) to 1
        cum_matrix[:,:,1] = 0
        cum_matrix[:,:,2] = 1
        cum_matrix[0, 0:band_len+1, 0] = 0
        cum_matrix[0, 0:band_len+1, 2] = 0
        pre_step_matrix = np.zeros((short_len + 1, long_len + 1), dtype = int)
        '''
        matrix for recording the best path. Each of the cells contains one of three
        possible integer: 0 ,1 ,2 indecating the corresponding value in cum_matrix
        (cum_matrix[i,j]) came from the cum_matrix[i-1,j], cum_matrix[i - 1,j - 1],
        and cum_matrix[i, j - 1] respectively.
        '''
        #kmer_bound = [False, False] + \
        #             [not np.array_equal(long[i],long[i-1])
        #              for i in range(1, len(long))]
        
        for i in range(1, short_len + 1):
            #j_start = int((i - 1) * band_move)
            for j in range(1, long_len + 1):    
            #for j in range(j_start, j_start + band_len):
                current_cost = self._cost(self.junction_squiggle[i-1], 
                                *self.candidate_squiggle[j-1], 
                                dist_type = self.dist_type)
                
                
                if cum_matrix[i-1,j-1,2] < min_dwell:
                    new_mu, new_n = _update_mean(current_cost,cum_matrix[i-1,j,1], cum_matrix[i-1,j,2])        
                    
                    pre_step_matrix[i, j] = 0
                    cum_matrix[i, j, 0] = cum_matrix[i-1, j, 0]   
                    cum_matrix[i, j, 1] = new_mu
                    cum_matrix[i, j, 2] =  new_n     
                else:
                    new_mu, new_n = _update_mean(current_cost,cum_matrix[i-1,j,1], cum_matrix[i-1,j,2])
                    new_values = [new_mu+cum_matrix[i-1,j,0],
                                cum_matrix[i-1,j-1,0] + cum_matrix[i-1,j-1,1] + current_cost]
                    
                    
                    best = np.random.choice(\
                                    np.where(new_values==min(new_values))[0])
                    
                    pre_step_matrix[i, j] = best
                    cum_matrix[i, j, 0] = cum_matrix[i-1, j, 0] if best == 0 else cum_matrix[i-1,j-1,0] + cum_matrix[i-1,j-1,1]         
                    cum_matrix[i, j, 1] = new_mu if best == 0 else current_cost     
                    cum_matrix[i, j, 2] =  new_n if best == 0 else 1               
        
        cum_matrix = cum_matrix[:,:,0] + cum_matrix[:,:,1]
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

        return dtw_res(best_score,best_path[::-1],cum_matrix, 
                       self.candidate_squiggle,
                        self.junction_squiggle,
                        self.band_prop,
                        self.max_z)

class dtw_local_alignment_mean_of_seg_mean(dtw):
    
    def dtw_alignment(self, min_dwell=1):
        def _update_mean(x , pre_mu, pre_n):
            return (pre_mu*pre_n+x)/(pre_n+1),pre_n+1

        def _move_to_new_col(pre_mean, pre_n, pre_col_mean,pre_col_n, current_cost):
            '''
            return cost and new cell when path move to a new col
            '''
            cur_n = pre_n+1
            cur_col_n = 1
            if pre_n == 0:
                cur_mean=current_cost+pre_mean
            else:
                cur_mean = _update_mean(current_cost, pre_mean, pre_n)[0]
            cur_col_mean = current_cost
            return [cur_mean,cur_n,cur_col_mean,cur_col_n]


        def _stay_col(pre_mean, pre_n, pre_col_mean,pre_col_n, current_cost):
            '''
            return cost and new cell when path stay col
            '''
            cur_n = pre_n
            cur_col_n = pre_col_n + 1
            cur_col_mean = _update_mean(current_cost, pre_col_mean, pre_col_n)[0] 
            if pre_n == 0:
                cur_mean=current_cost+pre_mean
                cur_n=1
                cur_col_n=1
            else:
                cur_mean = (pre_mean*pre_n-pre_col_mean+cur_col_mean)/pre_n
            return [cur_mean,cur_n,cur_col_mean,cur_col_n]

        # seed for tie-breaker when trace back the DTW path
        np.random.seed(seed=2021)

        short_len = len(self.junction_squiggle)
        long_len = len(self.candidate_squiggle)
        band_len = int(np.ceil(long_len* self.band_prop))
        band_move = (long_len - band_len)/short_len
        # cummulated matrix for DTW
        cum_matrix = np.full((short_len + 1, long_len + 1,4), np.inf)
        '''
        each cell in cum_matrix:
            [mean of current col means, number of current cols, mean of currect col, number of cell in correct col]
        '''
        # init the count (for mean calculation) to 1
        cum_matrix[:,:,1] = 0
        cum_matrix[:,:,2] = 0
        cum_matrix[:,:,3] = 0
        cum_matrix[0, 0:band_len+1, 0] = 0
        #cum_matrix[0, 0:band_len+1, 1] = 0
        pre_step_matrix = np.zeros((short_len + 1, long_len + 1), dtype = int)
        '''
        matrix for recording the best path. Each of the cells contains one of three
        possible integer: 0 ,1 ,2 indecating the corresponding value in cum_matrix
        (cum_matrix[i,j]) came from the cum_matrix[i-1,j], cum_matrix[i - 1,j - 1],
        and cum_matrix[i, j - 1] respectively.
        '''
        #kmer_bound = [False, False] + \
        #             [not np.array_equal(long[i],long[i-1])
        #              for i in range(1, len(long))]
        
        for i in range(1, short_len + 1):
            #j_start = int((i - 1) * band_move)
            for j in range(1, long_len + 1):    
            #for j in range(j_start, j_start + band_len):
                current_cost = self._cost(self.junction_squiggle[i-1], 
                                *self.candidate_squiggle[j-1], 
                                dist_type = self.dist_type)
                
                if cum_matrix[i-1,j-1,3] < min_dwell:
                    new_cell_stay = _stay_col(cum_matrix[i-1,j,0], cum_matrix[i-1,j,1], cum_matrix[i-1,j,2],cum_matrix[i-1,j,3], current_cost)                    
                    pre_step_matrix[i, j] = 0
                    cum_matrix[i, j, :] = new_cell_stay
                else:
                    new_cell_move = _move_to_new_col(cum_matrix[i-1,j-1,0], cum_matrix[i-1,j-1,1], cum_matrix[i-1,j-1,2],cum_matrix[i-1,j-1,3], current_cost)
                    new_cell_stay = _stay_col(cum_matrix[i-1,j,0], cum_matrix[i-1,j,1], cum_matrix[i-1,j,2],cum_matrix[i-1,j,3], current_cost)
                    new_values = [new_cell_stay[0],new_cell_move[0]]
                    best = np.random.choice(\
                                    np.where(new_values==min(new_values))[0])
                    
                    pre_step_matrix[i, j] = best
                    cum_matrix[i, j, :] = new_cell_stay if best==0 else  new_cell_move

        cum_matrix = cum_matrix[:,:,0]
        
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

        return dtw_res(best_score,best_path[::-1],cum_matrix, 
                       self.candidate_squiggle,
                        self.junction_squiggle,
                        self.band_prop,
                        self.max_z)

class dtw_res(object):
    '''
    Result from dtw:
        Attr:
            score
            path
            cum_matrix
        Methods:
            plot alignment
            plot warping path
    '''
    def __init__(self, score, path, cum_matrix,candidate_squiggle,junction_squiggle,band_prop,max_z):
        self.dtw_score = score
        self.dtw_path = path
        self.cum_matrix = cum_matrix
        self.candidate_squiggle = candidate_squiggle
        self.junction_squiggle = junction_squiggle
        self.band_prop = band_prop
        self.max_z = max_z
    
    def plot_alignment(self, save_fig='', **kwargs):        
        if 'figsize' in kwargs.keys():
            plt.figure(figsize=kwargs['figsize'])
        else:
            plt.figure(figsize=(30,5))
            
        plt.plot(self.junction_squiggle)
        plt.plot(self.candidate_squiggle[[np.array(self.dtw_path[:,1])-1]][:,0])
        if 'xlabel' in kwargs.keys():
            plt.xlabel(kwargs['xlabel'])
        if 'ylabel' in kwargs.keys():
            plt.ylabel(kwargs['ylabel'])
        if 'title' in kwargs.keys():
            plt.title(kwargs['title'])
            
        if kwargs.get('legend'):
            plt.legend()
            
        if save_fig:
            plt.savefig(save_fig)
        else:
            plt.show()
    def plot_cum_matrix(self, save_fig='', **kwargs):
        fig, axes = plt.subplots(nrows=1, figsize=(10,8))
        
        plot_mat = self.cum_matrix
        if len(plot_mat.shape) > 2:
            plot_mat = plot_mat[:,:,0]
        #plot_mat[plot_mat>10] = 10
        pos = axes.imshow(plot_mat, cmap='hot', interpolation='nearest',aspect='auto')
        fig.colorbar(pos,ax = axes)
        axes.plot(self.dtw_path[:,1], self.dtw_path[:,0])
        axes.set_title("Alignment path",fontsize=30, pad = 20)
        if kwargs.get('show_value'):
            for x_index, x in enumerate(range(plot_mat.shape[0])):
                for y_index, y in enumerate(range(plot_mat.shape[1])):
                    label = '{:.2f}'.format(plot_mat[x_index, y_index])
                    axes.text(y_index, x_index, label, color='white', ha='center', va='center')            
        if 'xlabel' in kwargs.keys():
            plt.xlabel(kwargs['xlabel'])
        if 'ylabel' in kwargs.keys():
            plt.ylabel(kwargs['ylabel'])
        if 'title' in kwargs.keys():
            plt.title(kwargs['title'])
        
        
        axes.set_xticklabels(np.concatenate(([0,0], self.candidate_squiggle[:,0])))
        axes.set_yticklabels(np.concatenate(([0,0], self.junction_squiggle)))        
        if not save_fig:
            plt.show()
        else:
            fig.savefig(save_fig)
        plt.close()


def main():
    # print('testing...\n')
    # short = np.array([[1,1],[1,1],[1,1],[3,1],[4,1],[5,1],[5,1],[5,1],[7,1],[5,1],[5,1],[4,1],[5,1]])
    # long = np.array([1,1,1,3,4,5,7,5,4,5])
    # print(short)
    # print(long)
    # path , score = dtw_local_alignment_max_sum(long, short)
    # path , score = dtw_local_alignment_max_sum_band_flipped(long, short)
    # plt.plot(long)
    # plt.plot(np.array(path)[:,1]-1, short[[np.array(path)[:,0]-1]][:,0])
    # plt.savefig('test.png')
    # print(path)
    
    short = np.array([[1,1],[1,1],[1,1],[3,1],[4,1],[5,1],[5,1],[5,1],[7,1],[5,1],[5,1],[4,1],[5,1]])
    long = np.array([1,1,1,3,4,5,7,5,4,5])
    print('sequence 1 (short)\n')
    print(short)
    print('sequence 2 (long)\n')
    print(long)

    print('test: dtw local alignment, max sum of cost\n')
    path , score = dtw_local_alignment_max_sum(long, short)
    plt.plot(long)
    plt.plot(np.array(path)[:,1]-1, short[[np.array(path)[:,0]-1]][:,0])
    plt.savefig('test_local_max_sum.png')

    print('test: dtw global alignment, max sum of mean of each segment\n')
    path , score = dtw_local_alignment_max_sum(long, short)
    plt.plot(long)
    plt.plot(np.array(path)[:,1]-1, short[[np.array(path)[:,0]-1]][:,0])
    plt.savefig('test_local_max_sum.png')
    print(path)


    print('test: dtw local alignment, max sum of mean of each segment\n')
    path , score = dtw_local_alignment_max_sum(long, short)
    plt.plot(long)
    plt.plot(np.array(path)[:,1]-1, short[[np.array(path)[:,0]-1]][:,0])
    plt.savefig('test_local_max_sum.png')
    print(path)

if __name__ == '__main__':
    main()

    
            
