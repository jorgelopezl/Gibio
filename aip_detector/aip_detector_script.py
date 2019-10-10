#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 15:51:56 2019

@author: Jorge López Luna da Silva. UTN FRBA
"""

# Python Code "Versatile Detector of Pseudo-periodic Patterns" (A. Santini, M. Llamedo Soria, E. Diez)

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import wfdb
import scipy
from scipy.signal import lfilter
from scipy.signal import find_peaks

# Load the ECG
da , info = wfdb.rdsamp( '100', pb_dir='mitdb/', sampfrom=0 , channels=[0])
data = pd.DataFrame({'hart':da[:,0]})

# N is the pattern width
N = 40

# g(n) : Gaussian function
sigma = (N-1)/5
g = scipy.signal.gaussian( N , ((N-1)/5) )

# Pattern:  p = dg * g
dg = np.diff(g)
g = g[1:]    # Delete last element
p = dg*g

# Normalize_pattern = true
maxim = np.max(np.abs(p))   
p = p/maxim

# Filtering
s = np.array(data['hart'])
rise_det = lfilter( p , 1 , np.flip(s))
rise_det = lfilter( p , 1 , np.flip(rise_det))

# Low Pass filter
example = 120e-3
pattern_size = 2*np.round(example/2*info['fs'])+1
lp_size = np.round(1.2*pattern_size)  
vector = np.ones(int(lp_size))
vector = vector/lp_size
rise_det = lfilter( vector, 1, np.flip(np.abs(rise_det)) )
rise_det = lfilter( vector, 1, np.flip(rise_det) )



# Here we obtain the maximum values above 30% of the observations
initial_thr = 30
actual_thr = np.percentile(rise_det, initial_thr)  
peaks_loc, _ = find_peaks(rise_det, height=actual_thr)

plt.figure()
plt.title("Results")
plt.plot(rise_det, alpha=0.5, color='green')
plt.plot(s, alpha=0.5, color='black')
plt.scatter(peaks_loc,rise_det[peaks_loc], alpha=0.9, color='blue', marker="X")  
plt.show()
    
# Histogram
q = list(range(1,100))
prctile_grid = np.percentile(rise_det[peaks_loc],q)
grid_step = np.median(np.diff(prctile_grid))
max_values = np.max(rise_det[peaks_loc])
thr_grid = np.arange(actual_thr,max_values,grid_step)
hist_max_values,a = np.histogram(rise_det[peaks_loc],thr_grid)

plt.figure()
plt.hist(rise_det[peaks_loc],bins = len(thr_grid))
plt.show()


# Statistical threshold
first_bin_idx = 1
thr_idx = find_peaks(hist_max_values)
thr_max = hist_max_values[thr_idx[0]]
thr_idx_expected = np.floor(np.dot(thr_idx[0],thr_max) *(1/np.sum(thr_max)))
aux_seq = np.array(range(1,len(thr_grid)))
hist_max_values = hist_max_values[aux_seq < thr_idx_expected]
min_hist_max_values = np.min(hist_max_values)
aux_seq = aux_seq[ aux_seq < thr_idx_expected ]
aux_seq = aux_seq[ aux_seq >= first_bin_idx ]
thr_min_idx = np.round(np.mean(np.nonzero(aux_seq & hist_max_values == min_hist_max_values)))
actual_thr = thr_grid[int(thr_min_idx)]
thr_grid = np.arange(actual_thr,max_values,grid_step)

# Here we plot the thresholds on the histogram
plt.axvline(x=actual_thr, color = 'red', label = 'Threshold_2')
plt.axvline(x=thr_grid[int(thr_idx_expected)], color = 'orange', label = 'Threshold_1')
plt.legend(framealpha=1, frameon=True);
plt.figure()
plt.hist(rise_det[peaks_loc],bins = thr_grid)
plt.show()

# Final results
peaks_loc, _ = find_peaks(rise_det, height=actual_thr)
plt.figure()
plt.title("Results")
plt.plot(rise_det, alpha=0.5, color='green')
plt.plot(s, alpha=0.5, color='black')
plt.scatter(peaks_loc,rise_det[peaks_loc], alpha=0.9, color='blue', marker="X")
plt.show()



##### Matlab #####
#Poner el test en la carpeta del ecg-kit
#Correr el test, con un breakpoint en el aip_detector
#
#Comandos:
#thr_grid = actual_thr : grid_step : max(max_values)
#histogram(max_values,thr_grid)