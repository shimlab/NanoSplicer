# -*- coding: utf-8 -*-
"""
Created on Sun May  8 16:55:12 2022

@author: yupeiy
"""

import numpy as np


# segmentation function
class junction_squiggle_preprocessing:
    def __init__(self, junction_squiggle):
        self.junction_squiggle_org = junction_squiggle.copy()
        self.segment_break_points = np.array([])
        self.junction_squiggle = junction_squiggle.copy()
    
    def reset_raw(self):
        """
        Reset the processed junction squiggle back to raw
        """
        self.junction_squiggle = self.junction_squiggle_org
        
    def get_segment_break_point(self, min_base_obs, running_stat_width, inplace = False):
        raw_cumsum = np.cumsum(np.concatenate([[0.0], self.junction_squiggle]))
        
        # get difference between all neighboring running_stat_width regions
        candidate_poss = np.argsort(np.abs(
            (2 * raw_cumsum[running_stat_width:-running_stat_width]) -
            raw_cumsum[:-2*running_stat_width] -
            raw_cumsum[2*running_stat_width:]))[::-1]
    
        cpts = [candidate_poss[0]]
        blacklist_pos = set()
    
        for pos in candidate_poss[1:]:
            if pos not in blacklist_pos:
                cpts.append(pos)
                blacklist_pos.update(range(
                    pos-min_base_obs+1, pos+min_base_obs+1))
        
        if inplace:
            self.segment_break_points = np.array(cpts) + running_stat_width
        else:    
            return np.array(cpts) + running_stat_width
    def segment_summary(self, n_seg, min_base_obs, 
                            running_stat_width,  method = 'median', inplace=False):
        '''
        Get an array of summary value in each segment:
        Parameters
        ----------
        n_seg : INT
            number of segment seleted, the number shouldn't be more than 
            the number of break_points stored in self.segment_break_points.
        method : STR, optional
            how to summarise the values within a segment, choose from 'median' or 'mean'. The default is 'median'.

        Returns
        -------
        np.array of segment summaries
        '''
        if not len(self.segment_break_points):
            self.get_segment_break_point(min_base_obs, 
                                         running_stat_width, 
                                         inplace=True)

        if method == 'median':
            rst = np.array(
                [[np.median(x), len(x)] for x in np.split(self.junction_squiggle,np.sort(self.segment_break_points[:n_seg]))])
        elif method == 'mean':
            rst = np.array(
                [[np.mean(x), len(x)] for x in np.split(self.junction_squiggle,np.sort(self.segment_break_points[:n_seg]))])
        
        if inplace:
            self.junction_squiggle_segment = rst
        else:
            return rst
        
    def median_smoother(self, window_size = 3, inplace = False):
        '''
        Denoise the original squiggle before DTW

        Parameters
        ----------
        junction_squiggle : NP.ARRAY
            junction squiggle
        window_size : INT, optional
            Sliding window size to do the median denoise. The default is 3.

        Returns
        -------
        np.array of median denoiced squiggle
        '''
        # add padding
        pad_start = np.repeat(self.junction_squiggle[0], np.ceil((window_size - 1) /2))
        pad_end = np.repeat(self.junction_squiggle[-1], np.floor((window_size - 1) /2))
        self.junction_squiggle = np.hstack([pad_start,self.junction_squiggle,pad_end])
        
        if inplace:
            self.junction_squiggle = np.array(
                [np.median(self.junction_squiggle[i:i+window_size]) 
                    for  i in range(len(self.junction_squiggle) - window_size + 1)])
        else:
            return np.array(
                [np.median(self.junction_squiggle[i:i+window_size]) 
                    for  i in range(len(self.junction_squiggle) - window_size + 1)])
        