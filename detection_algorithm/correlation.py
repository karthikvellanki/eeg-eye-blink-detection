import numpy as np
from scipy.signal import *

from scipy.interpolate import interp1d

def compute_correlation(p_blinks_t, data_sig, chan_id, fs):
    total_p_blinks = len(p_blinks_t)
    corr_matrix = np.ones([total_p_blinks, total_p_blinks])
    pow_matrix = np.ones([total_p_blinks, total_p_blinks])
    for idx_i in range(total_p_blinks):
        for idx_j in range(idx_i+1,total_p_blinks):

            blink_i_left = data_sig[int(fs*p_blinks_t[idx_i,0]):int(fs*p_blinks_t[idx_i,1]), chan_id]
            blink_i_right = data_sig[int(fs*p_blinks_t[idx_i,1]):int(fs*p_blinks_t[idx_i,2]), chan_id]

            blink_j_left = data_sig[int(fs*p_blinks_t[idx_j,0]):int(fs*p_blinks_t[idx_j,1]), chan_id]
            blink_j_right = data_sig[int(fs*p_blinks_t[idx_j,1]):int(fs*p_blinks_t[idx_j,2]), chan_id]

            left_interp = interp1d(np.arange(blink_i_left.size), blink_i_left)
            compress_left = left_interp(np.linspace(0,blink_i_left.size-1, blink_j_left.size))
            right_interp = interp1d(np.arange(blink_i_right.size), blink_i_right)
            compress_right = right_interp(np.linspace(0,blink_i_right.size-1, blink_j_right.size))

            sigA = np.concatenate((compress_left, compress_right))
            sigB = np.concatenate((blink_j_left, blink_j_right))
            
            corr = np.corrcoef(sigA, sigB)[0,1]
            corr_matrix[idx_i, idx_j] = corr
            corr_matrix[idx_j, idx_i] = corr
            
            if np.std(sigA) > np.std(sigB):
                pow_ratio = np.std(sigA)/np.std(sigB)
            else:
                pow_ratio = np.std(sigB)/np.std(sigA)
            
            pow_matrix[idx_i, idx_j] = pow_ratio
            pow_matrix[idx_j, idx_i] = pow_ratio
            

    return corr_matrix, pow_matrix