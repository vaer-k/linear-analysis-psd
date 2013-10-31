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
    evoked.save(data_path + subj[:6] + '/' = subj[:-4] + '-ave.fif')

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
    inv.save(data_path + subj[:6] + '/' + subj[:-4] + '-inv.fif')

    # picks MEG gradiometers
    picks = fiff.pick_types(raw.info, meg=True, eeg=False, eog=True, stim=False, exclude=exclude)

    # Compute source power spectral density and save to file
    psd = compute_source_psd(raw, inverse_operator, method='dSPM', fmin=fmin, fmax=fmax, pick_normal=True)
    psd.save(data_path + subj[:6] + '/' + subj[:-4] + '-psd.fif')
