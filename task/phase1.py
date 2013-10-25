#! /usr/bin/python

"""
======================================
Linear analysis of task data - Phase1
======================================

This script will process data beginning with the first step after the generation of the forward solution and ending with the calculation of PSD values.

"""
# Authors: Vincent Rupp Jr. <<ruppjr@hawaii.edu>>; Morgan Hough <<morgan@gazzaleylab.ucsf.edu>>

import mne
import os
import time
import csv
from mne import fiff, write_cov
from mne.fiff import Raw, pick_types
from mne.minimum_norm import read_inverse_operator, compute_source_psd_epochs, apply_inverse_epochs, write_inverse_operator, make_inverse_operator

###############################################################################
# Set global parameters
data_path = os.getcwd() + '/' 
subjects_dir = os.environ['SUBJECTS_DIR']
age = raw_input('YA or OA?\n')
fmin = float(raw_input('fmin:'))
fmax = float(raw_input('fmax:')) 
list_num = raw_input('Which list?\n')

event_id, tmin, tmax = 1, 0.0, 4.0
snr = 1.0 
lambda2 = 1.0 / snr ** 2
method = "dSPM" 

def intra(subj):
    '''
    Performs main process, including generation of inverse solution and PSD computation.
    '''
    print('Now beginning intra processing on ' + subj + '...\n') * 5

    # Set function parameters
    fname_raw = data_path + subj[:6] + '/' + subj 
    fname_fwd = data_path + subj[:6] + '/' + subj[:-4] + '-fwd.fif'

    # Load data 
    raw = fiff.Raw(fname_raw) 
    forward_meg = mne.read_forward_solution(fname_fwd) 

    # Estimate noise covariance from the raw data 
    precov = mne.compute_raw_data_covariance(raw, reject=dict(eog=150e-6)) 
    write_cov(data_path + subj + '/' + subj + '-cov.fif', precov) 

    # Find events from raw file
    events = mne.find_events(raw, stim_channel='STI 014')

    # Write events to file
    mne.write_events(subj + '-eve.fif', events)

    # Set up pick list: (MEG minus bad channels)
    include []
    exclude = raw.info['bads']
    picks = fiff.pick_types(raw.info, meg=True, eeg=False, stime=True, eog=True, include=include, exclude=exclude)

    # Read epochs and remove bad epochs
    epochs = mne.Epochs(raw, events, event_id, tmin, tmax, proj=True, picks=picks, baseline=(None, 0), preload=True, reject=dict(grad=4000e-13, mag=4e-12, eog=150e-6))    

    # Average epochs to produce an evoked dataset, then write to disk
    evoked = epochs.average()
    evoked.save(data_path + subj[:6] + '/' = subj + '-ave.fif')

    # Regularize noise cov
    cov = mne.cov.regularize(precov, evoked.info, grad=4000e-13, mag=4e-12, eog=150e-6, proj=True)

    # Restrict forward solution as necessary for MEG
    restricted_fwd = mne.fiff.pick_types_forward(forward_meg, meg=True, eeg=False) 

    # Make inverse operator
    info = evoked.info
    inverse_operator = make_inverse_operator(info, restricted_fwd, cov, loose=None, depth=0.8)

    # Pull data for averaging later
    epc_array = epochs.get_data()

    # Compute the inverse solution
    inv = apply_inverse(evoked, inverse_operator, lambda2, "dSPM", pick_normal=False)

    

#################################################################################
#################################################################################
	
    # define frequencies of interest
    bandwidth = 4.  # bandwidth of the windows in Hz
    
    # compute source space psd in label
    
    # Note: By using "return_generator=True" stcs will be a generator object
    # instead of a list. This allows us so to iterate without having to
    # keep everything in memory.
    
    psd = compute_source_psd_epochs(epochs, inverse_operator, lambda2=lambda2,
                                     method=method, fmin=fmin, fmax=fmax,
                                     bandwidth=bandwidth, label=label, return_generator=False)
    
    epoch_num = 1
    epoch_num_str = str(epoch_num)
    for i in psd:
    	i.save(data_path + subj + '/' + 'tmp' + '/' + label_name[3:] + '_dspm_snr-1_PSD'+ epoch_num_str)
    	epoch_num = epoch_num + 1
        epoch_num_str = str(epoch_num)
    
    if hemi == 'left':
        filelist = [ f for f in os.listdir(data_path + subj + '/tmp') if f.endswith("-rh.stc") ]
        for f in filelist:
            os.remove(data_path + subj + '/tmp/' + f)
    	keepers = [ f for f in os.listdir(data_path + subj + '/tmp') if f.endswith("-lh.stc") ]
        for f in keepers:
            src = f
            os.rename(data_path + subj + '/tmp/' + src,data_path + subj + '/psd/' + src)
    
    elif hemi == 'right':
        filelist = [ f for f in os.listdir(data_path + subj + '/tmp') if f.endswith("-lh.stc") ]
        for f in filelist:
            os.remove(data_path + subj + '/tmp/' + f)
    	keepers = [ f for f in os.listdir(data_path + subj + '/tmp') if f.endswith("-rh.stc") ]
        for f in keepers:
            src = f
            os.rename(data_path + subj + '/tmp/' + src,data_path + subj + '/psd/' + src)
   
 
    # This code computes the average PSDs of each epoch. Each PSD file is an array of shape N_vertices*N_frequencies. This code averages the PSD value of each vertex together and outputs the average PSD value of each frequency. Then, it averages the PSD values of each epoch, outputting one average PSD value per frequency value, i.e., this is the average across epochs.
    
    n_epochs = len(epc_array)
    for i, stc in enumerate(psd):
        if i >= n_epochs:
            break
    
        if i == 0:
            psd_avg = np.mean(stc.data, axis=0)
        else:
            psd_avg += np.mean(stc.data, axis=0)
    
    print('Length of psd for subject ' + subj + ' is ' + str(len(psd)) + '.')
    print('Number of epochs for subject ' + subj + ' is ' + str(n_epochs) + '.')
   
    if len(psd) != 0:
        psd_avg /= n_epochs
    
    # Compute variance for each epoch and then variance across epochs 
    
    n_epochs = len(epc_array)
    for i, stc in enumerate(psd):
        if i >= n_epochs:
            psd_var = np.array()
	    break
        
        if i == 0:
            psd_var = np.var(stc.data, axis=0)
        else:
            psd_var = np.vstack((psd_var,np.var(stc.data, axis=0)))
    
    if len(psd) >= 2:
        tot_var = np.var(psd_var, axis=0)

    if len(psd) <= 1:
	failed_subj = subj
	print(failed_subj + ' failed. No PSD values calculated, likely because all epochs were rejected.')
	return failed_subj, failed_subj, failed_subj

    if len(psd) >= 2:
        return (psd_avg, tot_var, len(psd_avg))


failed_list = []
# List subjects in cwd
for subject in os.listdir(os.getcwd()):
# Operate only on specified age group
    if subject.startswith(age):
# Check if preprocessing has been done (fwd solutions are last step of preprocessing) and process subjects with intra()
        if os.path.isfile(str(os.getcwd()) + '/' + subject + '/' + subject + '_list' + list_num + '_raw_sss-ico-4-fwd.fif'):
            for i in os.listdir(os.getcwd()):
                if i.endswith(list_num + '_raw_sss.fif'):
                    rawfile = i
                    individual_avg_psd, individual_var_var, n_freqs = intra(rawfile)
# Verify subject did not fail PSD calculations and calculate combined PSD averages and variances. Set up empty array for combined_avg_psd and stack combined_var_var into one vertical array
                    if type(individual_avg_psd) == np.ndarray: 
                        if not 'combined_avg_psd' in locals():
                            combined_avg_psd = np.zeros(n_freqs)
                        else:
                            combined_avg_psd = individual_avg_psd + combined_avg_psd
                        if not 'combined_var_var' in locals():
                            combined_var_var = individual_var_var
                        else:
                            combined_var_var = np.vstack((combined_var_var, individual_var_var))		
                    else:
                        failed_list.append(subject)
        else:
            failed_list.append(subject)
            print('Subject ' + subject + ' has not yet been preprocessed. (No forward solution was found)')

print('The following subjects failed PSD calculations:')
print(failed_list)

# Average across subjects
n_subj = len(os.listdir(os.getcwd())) - len(failed_list)
avg_psd = combined_avg_psd / n_subj

# Average across frequencies
final_avg_psd = np.mean(avg_psd) 

# Compute variance across subjects
var_var = np.var(combined_var_var, axis=0)

# Compute variance across frequencies
final_var = np.var(var_var)

# Compute linear regression
step = 0.25
xi = np.arange(fmin, fmax, step)   
print('xi length:' + str(len(xi)))
A = np.array([ xi, np.ones(len(xi))]) 
print('A.T length:' + str(len(A.T)))
y = avg_psd
print('y length:' + str(len(y)))
w = np.linalg.lstsq(A.T, y)[0] #obtaining the parameters
print('w:')
print(w)

# plotting the line
line = w[0]*xi+w[1] #regression line
pl.plot(xi,line,'r-',xi,y,'o')
pl.show()

print('avg_psd (across subjects) follows:')
print(avg_psd)
print('var_var (across subjects) follows:')
print(var_var)
print('final_avg_psd (across freqs) follows:')
print(final_avg_psd)
print('final_var (across freqs) follows:')
print(final_var)

# Save averages and variances across subjects in csv document
if age == 'YA':
    np.savetxt('psd_avg_YA.csv', avg_psd, delimiter=',')
    np.savetxt('var_var_YA.csv', var_var, delimiter=',')
if age == 'OA':
    np.savetxt('psd_avg_OA.csv', avg_psd, delimiter=',')
    np.savetxt('var_var_OA.csv', var_var, delimiter=',')
