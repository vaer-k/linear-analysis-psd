#! /usr/bin/python

"""
=====================================
Linear analysis of task data - Phase0
=====================================

Phase0 is used to generate a list of raw files to be processed.

"""
# Authors: Vincent Rupp Jr. <<ruppjr@hawaii.edu>>; Morgan Hough <<morgan@gazzaleylab.ucsf.edu>>

import os
import re

##############################################################

# Initialize 
subj_list = []
raw_lists = []

def filegen():
    # Get subjects in cwd (should be MEG subjects dir) corresponding to desired age group
    for subject in os.listdir(os.getcwd()):
        if subject.startswith(age):
            subj_list.append(subject)
    
    # Generate a list of raw list files to be processed
    m = re.compile('list..?_raw_sss.fif')
    for subject in subj_list:
        l = m.findall(str(os.listdir(os.getcwd() + '/' + subject)))
        raw_lists.append(l)
