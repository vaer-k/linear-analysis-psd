#! /usr/bin/python

"""
==========================
   linear-analysis-psd
==========================

This program will compute average PSD and variance of across subjects in specified frequency range, the perform analysis using Welche's t-test. Run this from the directory of your MEG subjects.

"""
# Authors: Vincent Rupp Jr. <<ruppjr@hawaii.edu>>; Morgan Hough <<morgan@gazzaleylab.ucsf.edu>>

print __doc__

import numpy as np
import pylab as pl
import scipy
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
age = 'YA'
#age = raw_input('YA or OA?\n')
label_name = 'lh.BA45'
#label_name = raw_input('Which region label would you like to compute PSD for?\n')
fmin = 80.0
fmax = 100.0
#fmin = float(raw_input('fmin:'))
#fmax = float(raw_input('fmax:')) 

event_id, tmin, tmax = 1, 0.0, 4.0
snr = 1.0 
lambda2 = 1.0 / snr ** 2
method = "dSPM" 

	
def intra(subj):
    '''
    Performs initial computations within subject and returns average PSD and variance of all epochs.
    '''
    print('Now beginning intra processing on ' + subj + '...\n') * 5

    # Set function parameters
    fname_label = '/home/vrupp/data/search_fMRI/' + subj + '/' + 'label/%s.label' % label_name
    fname_raw = data_path + subj + '/' + subj + '_rest_raw_sss.fif'
    if os.path.isfile(data_path + subj + '/' + subj + '_rest_raw_sss-ico-4-fwd.fif'): 
	fname_fwd = data_path + subj + '/' + subj + '_rest_raw_sss-ico-4-fwd.fif'
    else: 
        print('Subject ' + subj + ' does not have a ico-4-fwd.fif on file.')	

    if label_name.startswith('lh.'):
    	hemi = 'left'
    elif label_name.startswith('rh.'):
    	hemi = 'right'	
    
    # Load data
    label = mne.read_label(fname_label)
    raw = fiff.Raw(fname_raw)
    forward_meg = mne.read_forward_solution(fname_fwd)
    
    # Estimate noise covariance from teh raw data
    cov = mne.compute_raw_data_covariance(raw, reject=dict(eog=150e-6))
    write_cov(data_path + subj + '/' + subj + '-cov.fif', cov)
    
    # Make inverse operator
    info = raw.info
    inverse_operator = make_inverse_operator(info, forward_meg, cov, loose=None, depth=0.8)
    
    # Epoch data into 4s intervals
    events = mne.make_fixed_length_events(raw, 1, start=0, stop=None, 
    		duration=4.)
    
    # Set up pick list: (MEG minus bad channels)
    include = []
    exclude = raw.info['bads']
    picks = fiff.pick_types(raw.info, meg=True, eeg=False, stim=False, eog=True, include=include, exclude=exclude)
    
    # Read epochs and remove bad epochs
    epochs = mne.Epochs(raw, events, event_id, tmin, tmax, proj=True, picks=picks, baseline=(None, 0), preload=True, reject=dict(grad=4000e-13, mag=4e-12, eog=150e-6))
    
    # Pull data for averaging later
    epc_array = epochs.get_data()
    
    # Compute the inverse solution
    inv = apply_inverse_epochs(epochs, inverse_operator, lambda2, method, label=label)
    
    #Need to add a line here to automatically create stc directory within subj
    
    epoch_num = 1
    epoch_num_str = str(epoch_num)
    for i in inv:
#    	i.save(data_path + subj + '/tmp/' + label_name[3:] + '_rest_raw_sss-oct-6-inv' + epoch_num_str)
	i.save(data_path + subj + '/tmp/' + label_name[3:] + '_rest_raw_sss-ico-4-inv' + epoch_num_str)
    	epoch_num = epoch_num + 1
    	epoch_num_str = str(epoch_num)
    
    # The following is used to remove the empty opposing hemisphere files
    # and then move the files to save into the appropriate directory
    
    if hemi == 'left':
    	filelist = [ f for f in os.listdir("/data/restMEG/" + subj + '/tmp') if f.endswith("-rh.stc") ]	
    	for f in filelist:
            os.remove("/data/restMEG/" + subj + '/tmp/' + f)
    	keepers = [ f for f in os.listdir("/data/restMEG/" + subj + '/tmp') if f.endswith("-lh.stc") ]
    	for f in keepers:
    	    src = f 
            os.rename("/data/restMEG/" + subj + '/tmp/' + src,"/data/restMEG/" + subj + '/inv/' + src)
    
    elif hemi == 'right':
    	filelist = [ f for f in os.listdir("/data/restMEG/" + subj + '/tmp') if f.endswith("-lh.stc") ]
        for f in filelist:
            os.remove("/data/restMEG/" + subj + '/tmp/' + f)
    	keepers = [ f for f in os.listdir("/data/restMEG/" + subj + '/tmp') if f.endswith("-rh.stc") ]
        for f in keepers:
            src = f 
            os.rename("/data/restMEG/" + subj + '/tmp/' + src,"/data/restMEG/" + subj + '/inv/' + src)
    
    
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
    	i.save('/data/restMEG/' + subj + '/' + 'tmp' + '/' + label_name[3:] + '_dspm_snr-1_PSD'+ epoch_num_str)
    	epoch_num = epoch_num + 1
        epoch_num_str = str(epoch_num)
    
    if hemi == 'left':
        filelist = [ f for f in os.listdir("/data/restMEG/" + subj + '/tmp') if f.endswith("-rh.stc") ]
        for f in filelist:
            os.remove("/data/restMEG/" + subj + '/tmp/' + f)
    	keepers = [ f for f in os.listdir("/data/restMEG/" + subj + '/tmp') if f.endswith("-lh.stc") ]
        for f in keepers:
            src = f
            os.rename("/data/restMEG/" + subj + '/tmp/' + src,"/data/restMEG/" + subj + '/psd/' + src)
    
    elif hemi == 'right':
        filelist = [ f for f in os.listdir("/data/restMEG/" + subj + '/tmp') if f.endswith("-lh.stc") ]
        for f in filelist:
            os.remove("/data/restMEG/" + subj + '/tmp/' + f)
    	keepers = [ f for f in os.listdir("/data/restMEG/" + subj + '/tmp') if f.endswith("-rh.stc") ]
        for f in keepers:
            src = f
            os.rename("/data/restMEG/" + subj + '/tmp/' + src,"/data/restMEG/" + subj + '/psd/' + src)
   
 
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
    
    # Compute variance for each epoch and then variance of those results
    
    n_epochs = len(epc_array)
    for i, stc in enumerate(psd):
        if i >= n_epochs:
            psd_var = np.array()
	    break
        
        if i == 0:
            psd_var = np.var(stc.data, axis=0)
        else:
            psd_var = np.vstack((psd_var,np.var(stc.data, axis=0)))
    
    if len(psd) != 0:
        tot_var = np.var(psd_var, axis=0)

    if len(psd) == 0:
	failed_subj = subj
	print(failed_subj + ' is a no go. No PSD values calculated, likely because all epochs were rejected.')
	time.sleep(30)
	return failed_subj, failed_subj, failed_subj

    if len(psd) != 0:
        return (psd_avg, tot_var, len(psd_avg))


failed_list = []
# List subjects in cwd
for subject in os.listdir(os.getcwd()):
# Operate only on specified age group
        if subject.startswith(age):
	    print(subject) * 5
# Check if preprocessing has been done (fwd solutions are last step of preprocessing) and process subjects with intra()
            if os.path.isfile(str(os.getcwd()) + '/' + subject + '/' + subject + '_rest_raw_sss-ico-4-fwd.fif'):
		individual_avg_psd, individual_var_var, n_freqs = intra(subject)
# Verify subject did not fail PSD calculations and calculate combined PSD averages and variances. Set up empty array for combined_avg_psd and stack combined_var_var into one vertical array
		print(type(individual_avg_psd))
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
                print('Subject ' + subject + ' has not yet been preprocessed.')

print('The following subjects failed PSD calculations:')
print(failed_list)

n_subj = len(os.listdir(os.getcwd())) - len(failed_list)
avg_psd = combined_avg_psd / n_subj
var_var = np.var(combined_var_var, axis=0)
print('avg_psd follows:')
print(avg_psd)
print('var_var follows:')
print(var_var)
np.savetxt('psd_avg.csv', avg_psd, delimiter=',')
np.savetxt('var_var.csv', var_var, delimiter=',')
