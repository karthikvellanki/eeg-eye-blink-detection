import os
import numpy as np
from scipy.signal import *
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster

from init_args import args_init
from correlation import compute_correlation
from identify_peaks import identify_signal_peaks


data_path = '../sample_data'
file_idx = 7

fs = 250.0
chan_id = 1

blink_len_max = 2.0
blink_len_min = 0.3


delta_init = 100

correlation_threshold_1 = 0.2
correlation_threshold_2 = 0.7

std_threshold_window = int(5*fs)


def lowpass(sig, fc, fs, butter_filt_order):
    B,A = butter(butter_filt_order, np.array(fc)/(fs/2), btype='low')
    return lfilter(B, A, sig, axis=0)

sample_data_file = 'sample_eyeblink_data.csv'

data_sig = np.loadtxt(open(os.path.join(data_path,sample_data_file), "rb"), delimiter=";", skiprows=1, usecols=(0,1,2))

data_sig[:,1] = lowpass(data_sig[:,1], 10, fs, 4)
data_sig[:,2] = lowpass(data_sig[:,2], 10, fs, 4)

time_min = data_sig[0,0]
time_max = data_sig[-1,0]

data_len = len(data_sig)

def moving_standard_deviation(data_sig, chan_id, fs):
    std_length = int(0.5*fs)
    
    running_std = np.zeros([data_len,1])
    idx = 0
    while(idx < len(data_sig) - std_length):
        running_std[idx] = np.std(data_sig[idx:(idx + std_length), chan_id])
        idx = idx + 1
    running_std[idx:-1] = running_std[idx-1]
            
    return running_std
    

def refine_blink_events(stat_min2, data_sig, chan_id):
    offset_t = 0.00
    win_size = 25
    win_offset = 10
    search_maxlen_t = 1.5


    offset_f = int(offset_t*fs)
    search_maxlen_f = int(search_maxlen_t*fs)
    iters = int(search_maxlen_f/win_offset)

    data_len = len(data_sig)
    p_blinks_t, p_blinks_val = [], []
    for idx in range(len(stat_min2)):
        x_indR = int(fs*stat_min2[idx,0]) + offset_f
        x_indL = int(fs*stat_min2[idx,0]) - offset_f
        start_index = max(0, int(fs*stat_min2[idx,0]) - std_threshold_window)
        end_index = min( int(fs*stat_min2[idx,0]) + std_threshold_window, data_len)
        stable_threshold = 2*min(running_std[start_index:end_index])
        min_val = stat_min2[idx,1];
        max_val = min_val;
        found1, found2 = 0, 0
        state1, state2 = 0, 0

        for iter in range(iters):
            if(x_indR + win_size > data_len):
                x_indR = x_indR - (x_indR + win_size - data_len)
            if(x_indL < 0):
                x_indL = 0
            if (np.std(data_sig[x_indR:x_indR+win_size, chan_id]) < stable_threshold) and state1==1 and data_sig[x_indR, chan_id]>min_val:
                found1 = 1
                max_val = max(data_sig[x_indR, chan_id],max_val)
            if (np.std(data_sig[x_indL:x_indL+win_size, chan_id]) < stable_threshold) and state2==1 and data_sig[x_indL + win_size, chan_id]>min_val:
                found2 = 1
                max_val = max(data_sig[x_indL + win_size, chan_id],max_val)
            if (np.std(data_sig[x_indR:x_indR+win_size, chan_id]) > 2.5*stable_threshold) and state1==0:
                state1 = 1
            if (np.std(data_sig[x_indL:x_indL+win_size, chan_id]) > 2.5*stable_threshold) and state2==0:
                state2 = 1
            if (found1==1) and data_sig[x_indR, chan_id] < (max_val + 2*min_val)/3:
                found1=0
            if (found2==1) and data_sig[x_indL + win_size, chan_id] < (max_val + 2*min_val)/3:
                found2=0
            if (found1==0):
                x_indR = x_indR + win_offset
            if (found2==0):
                x_indL = x_indL - win_offset;
            if found1==1 and found2==1:
                break
        if found1==1 and found2==1:
            if (x_indL + win_size)/fs > stat_min2[idx,0]:
                p_blinks_t.append([(x_indL)/fs, stat_min2[idx,0], x_indR/fs])
                p_blinks_val.append([data_sig[x_indL, chan_id], stat_min2[idx,1], data_sig[x_indR,chan_id]])         
            else:
                p_blinks_t.append([(x_indL + win_size)/fs, stat_min2[idx,0], x_indR/fs])
                p_blinks_val.append([data_sig[x_indL + win_size, chan_id], stat_min2[idx,1], data_sig[x_indR,chan_id]])
            

    p_blinks_t = np.array(p_blinks_t)        
    p_blinks_val = np.array(p_blinks_val)  
    
    return p_blinks_t, p_blinks_val

running_std = moving_standard_deviation(data_sig, chan_id, fs)

def main():
    args_chan1 = args_init(delta_init)

    for idx in range(len(data_sig[:,0])):
        identify_signal_peaks(data_sig[idx,0], data_sig[idx, chan_id], args_chan1)

    min_pts = np.array(args_chan1['mintab'])
    p_blinks_t, p_blinks_val = refine_blink_events(min_pts, data_sig, chan_id)
    corr_matrix, pow_matrix = compute_correlation(p_blinks_t, data_sig, chan_id, fs)

    blink_fp_idx = np.argmax(sum(corr_matrix))
    t = corr_matrix[blink_fp_idx,:] > correlation_threshold_1
    blink_index = [i for i, x in enumerate(t) if x]

    blink_template_corrmat = corr_matrix[np.ix_(blink_index,blink_index)]
    blink_template_powmat = pow_matrix[np.ix_(blink_index,blink_index)]
    blink_templates_corrWpower = blink_template_corrmat/blink_template_powmat

    blink_var = []
    for idx in blink_index:
        blink_var.append(np.var(data_sig[int(fs*p_blinks_t[idx,0]):int(fs*p_blinks_t[idx,2]), chan_id]))



    Z = linkage(blink_templates_corrWpower, 'complete', 'correlation')
    groups = fcluster(Z,2,'maxclust')

    grp_1_blinks_var = [blink_var[i] for i, x in enumerate(groups==1) if x]
    grp_2_blinks_var = [blink_var[i] for i, x in enumerate(groups==2) if x]
    if np.mean(grp_1_blinks_var) > np.mean(grp_2_blinks_var):
        selected_group = 1
    else:
        selected_group = 2
    template_blink_idx = [blink_index[i] for i, x in enumerate(groups==selected_group) if x]

    delta_new = 0
    for idx in template_blink_idx:
        delta_new = delta_new + min(p_blinks_val[idx,0], p_blinks_val[idx,2]) - p_blinks_val[idx,1]
    delta_new = delta_new/len(template_blink_idx)


    args_chan1 = args_init(delta_new/3.0)

    for idx in range(len(data_sig[:,0])):
        identify_signal_peaks(data_sig[idx,0], data_sig[idx, chan_id], args_chan1)

    min_pts = np.array(args_chan1['mintab'])
    p_blinks_t, p_blinks_val = refine_blink_events(min_pts, data_sig, chan_id)
    corr_matrix, pow_matrix = compute_correlation(p_blinks_t, data_sig, chan_id, fs)

        
    s_fc = (sum(corr_matrix))
    sort_idx = sorted(range(len(s_fc)), key=lambda k: s_fc[k])

    t = corr_matrix[sort_idx[-1],:] > correlation_threshold_2        
    blink_index1 = set([i for i, x in enumerate(t) if x])
    t = corr_matrix[sort_idx[-2],:] > correlation_threshold_2        
    blink_index2 = set([i for i, x in enumerate(t) if x])
    t = corr_matrix[sort_idx[-3],:] > correlation_threshold_2        
    blink_index3 = set([i for i, x in enumerate(t) if x])

    blink_index = list(blink_index1.union(blink_index2).union(blink_index3))

    blink_template_corrmat = corr_matrix[np.ix_(blink_index,blink_index)]
    blink_template_powmat = pow_matrix[np.ix_(blink_index,blink_index)]
    blink_templates_corrWpower = blink_template_corrmat/blink_template_powmat

    blink_var = []
    for idx in blink_index:
        blink_var.append(np.var(data_sig[int(fs*p_blinks_t[idx,0]):int(fs*p_blinks_t[idx,2]), chan_id]))


    Z = linkage(blink_templates_corrWpower, 'complete', 'correlation')
    groups = fcluster(Z,2,'maxclust')

    grp_1_blinks_var = [blink_var[i] for i, x in enumerate(groups==1) if x]
    grp_2_blinks_var = [blink_var[i] for i, x in enumerate(groups==2) if x]

    if np.mean(grp_1_blinks_var) > np.mean(grp_2_blinks_var) and np.mean(grp_1_blinks_var)/np.mean(grp_2_blinks_var) > 10:
        blink_index = [blink_index[i] for i, x in enumerate(groups==1) if x]
    elif np.mean(grp_2_blinks_var) > np.mean(grp_1_blinks_var) and np.mean(grp_2_blinks_var)/np.mean(grp_1_blinks_var) > 10:
        blink_index = [blink_index[i] for i, x in enumerate(groups==2) if x]

    detected_blinks = p_blinks_t[blink_index,:]

    print("Detected Blinks:")
    print(detected_blinks)

if __name__ == "__main__":
    main()
