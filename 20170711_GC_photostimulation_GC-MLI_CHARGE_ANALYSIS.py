# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 14:01:32 2017
@author: Kevin.Dorgans

>> load individual WCP voltage-clamp traces for each GC-MLI photostimulation ePSCs
>> gives a index name to each trace
>> sorts WT and Synpasin2KO data

OUTPUT GRAPHS
    mean WT and S2KO ePSCs following unitary GC photostimulation
    cumulative plot WT vs Syn2KO
    
IMPORTANT VARIABLES
    eEPSC_HISTOGRAM                      individual traces per unitary GC photostimulation
    eEPSC_HISTOGRAM_ID_LIST__            index for eEPSC_HISTOGRAM

"""


import os
import numpy as np
import pandas as pd
import scipy.stats as spst
from matplotlib import pyplot as plt


bin_size = 5 #bin size for charge calculation - ms
count = 0

eEPSC_HISTOGRAM = []
eEPSC_HISTOGRAM_ID_LIST__ = []

path = os.listdir(r'D:/DATA/GC-MLI SeqPatch/MLI PATCH/INDIVIDUALS_2')

for k in range(len(path)):

    A = np.array(pd.read_csv(r'D:/DATA\GC-MLI SeqPatch\MLI PATCH/INDIVIDUALS_2/'+str(path[k]), sep=' '), dtype=np.float64)
    count += len(A)
    if (len(A)) > 3:
        eEPSC_HISTOGRAM_ID_LIST__.append(str(path[k]))
        LIST__ = np.mean(np.array(A)/bin_size, axis=0)
        eEPSC_HISTOGRAM.append(LIST__-np.mean(LIST__[0:50]))


plt.figure(figsize=(3, 3))
list__ = np.mean([eEPSC_HISTOGRAM[i] for i in range(len(eEPSC_HISTOGRAM)) if ('S2KO' in eEPSC_HISTOGRAM_ID_LIST__[i]) == True], axis=0)
list_se = spst.sem([eEPSC_HISTOGRAM[i] for i in range(len(eEPSC_HISTOGRAM)) if ('S2KO' in eEPSC_HISTOGRAM_ID_LIST__[i]) == True], axis=0)

plt.plot(np.linspace(-500, 1500, len(list__)), list__, color='red')
plt.fill_between(np.linspace(-500, 1500, len(list__)), list__+list_se, list__-list_se, color='red', alpha=0.5)

list__ = np.mean([eEPSC_HISTOGRAM[i] for i in range(len(eEPSC_HISTOGRAM)) if ('S2KO' in eEPSC_HISTOGRAM_ID_LIST__[i]) == False], axis=0)
list_se = spst.sem([eEPSC_HISTOGRAM[i] for i in range(len(eEPSC_HISTOGRAM)) if ('S2KO' in eEPSC_HISTOGRAM_ID_LIST__[i]) == True], axis=0)

plt.plot(np.linspace(-500, 1500, len(list__)), list__, color='black')
plt.fill_between(np.linspace(-500, 1500, len(list__)), list__+list_se, list__-list_se, color='black', alpha=0.5)

plt.xlim(0, 150)
plt.ylabel("eEPSC Charge (fC)")
plt.xlabel('Time Post-Stimulus (ms)')
plt.tight_layout()


plt.figure(figsize=(3, 3))
list_cumsum = np.cumsum([eEPSC_HISTOGRAM[i] for i in range(len(eEPSC_HISTOGRAM)) if ('S2KO' in eEPSC_HISTOGRAM_ID_LIST__[i]) == True], axis=1)
list_cumsum_se = spst.sem(list_cumsum, axis=0)
list_cumsum = np.mean(list_cumsum, axis=0)
plt.plot(np.linspace(-500, 1500, len(list_cumsum)), list_cumsum, color='red')
plt.fill_between(np.linspace(-500, 1500, len(list_cumsum)), list_cumsum+list_cumsum_se, list_cumsum-list_cumsum_se, color='red', alpha=0.5)


plt.xlim(0, 150)
plt.ylabel("eEPSC Cumulative Charge (fC)")
plt.xlabel('Time Post-Stimulus (ms)')

list_cumsum = np.cumsum([eEPSC_HISTOGRAM[i] for i in range(len(eEPSC_HISTOGRAM)) if ('S2KO' in eEPSC_HISTOGRAM_ID_LIST__[i]) == False], axis=1)
list_cumsum_se = spst.sem(list_cumsum, axis=0)
list_cumsum = np.mean(list_cumsum, axis=0)
plt.plot(np.linspace(-500, 1500, len(list_cumsum)), list_cumsum, color='black')
plt.fill_between(np.linspace(-500, 1500, len(list_cumsum)), list_cumsum+list_cumsum_se, list_cumsum-list_cumsum_se, color='black', alpha=0.5)
plt.xlim(0, 150)
plt.ylabel("eEPSC Cumulative Charge (fC)")
plt.xlabel('Time Post-Stimulus (ms)')
plt.ylim(0, -8000)
plt.tight_layout()
