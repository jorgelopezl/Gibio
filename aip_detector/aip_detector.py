#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 15:51:56 2019

@author: Jorge López Luna da Silva. UTN FRBA
"""

# Adaptación a Python del "Versatile Detector of Pseudo-periodic Patterns" (A. Santini, M. Llamedo Soria, E. Diez)


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import wfdb
import scipy
from scipy.signal import lfilter
from scipy.signal import find_peaks

# Adaptación a Python del "Versatile Detector of Pseudo-periodic Patterns" (A. Santini, M. Llamedo Soria, E. Diez)


# Cargamos el archivo que contiene el ECG
da , info = wfdb.rdsamp( '100', pb_dir='mitdb/', sampfrom=0, sampto = 30000 , channels=[0])

datos = pd.DataFrame({'hart':da[:,0]})


# N es el largo del patron y deberìa obtenerse del diccionario que pasan como param
# se llama pattern_width 
N = 40

# Creamos la gaussiana g(n)
sigma = (N-1)/5
g = scipy.signal.gaussian( N , ((N-1)/5) )

# Creamos el patrón:  p = dg * g

dg = np.diff(g)
g = g[1:]    # Saco el ùltimo elemento para poder multiplicar con misma dimensión
p = dg*g

# Normalizo al patrón      normalize_pattern = true
maxim = np.max(np.abs(p))
p = p/maxim

# Filtramos 
s = np.array(datos['hart'])
rise_det = lfilter( p , 1 , np.flip(s))
rise_det = lfilter( p , 1 , np.flip(rise_det))

# Pasabajos
ejemplo = 120e-3
pattern_size = 2*np.round(ejemplo/2*info['fs'])+1
lp_size = np.round(1.2*pattern_size)  #Parametro de entrada
vector = np.ones(int(lp_size))
vector = vector/lp_size
rise_det = lfilter( vector, 1, np.flip(np.abs(rise_det)) )
rise_det = lfilter( vector, 1, np.flip(rise_det) )



# Obtenemos valores maximos de acuerdo a un unmbral

#maxim_rise_det = np.max(rise_det)
#print(maxim_rise_det)
#percent = 0.3 * maxim_rise_det
#print(percent)

initial_thr = 30
actual_thr = np.percentile(rise_det, initial_thr)  # Valor para el cual se encuentra el 30% de las observaciones
peaks, _ = find_peaks(rise_det, height=actual_thr)    
    
plt.figure()
plt.title("Resultado")
plt.plot(rise_det, alpha=0.5, color='green')
plt.plot(s, alpha=0.5, color='black')
plt.scatter(peaks,rise_det[peaks], alpha=0.9, color='blue', marker="X")    ###################
plt.show()
    
# Segundo threshold
q = list(range(1,100))
prctile_grid = np.percentile(peaks,q)
print(prctile_grid)
grid_step = np.median(np.diff(prctile_grid))
print("Step")
print(grid_step)
max_values = np.max(peaks)
print(max_values)
thr_grid = np.arange(actual_thr,max_values,grid_step)
print("Thr_grid")
print(len(thr_grid))
print(thr_grid)
hist_max_values,a = np.histogram(peaks,thr_grid)
print("%%%%%%%%%%%%%%%%%%%%%")
print(hist_max_values)
print(a)

plt.figure()
plt.hist(hist_max_values,bins = len(thr_grid))
plt.show()
    
    
    
    
    
    
    