# -*- coding: utf-8 -*-
"""
Created on Thu May 11 19:50:58 2017
@author: Kevin.Dorgans

>> loads MLI photostimulation-evoked firing acceleration data from *.csv files
>> sorts data between WT data and Synapsin2 Knockout data
>> computes ISI histograms and various parameters (CV, CV2, LvR)

OUTPUT FIGURES
    WT vs Syn2KO MLI instantaneous firing frequency after unitary GC photostimulation
    WT vs Syn2KO firing probability locked to previous spontaneous spike (spiking MLI) or stimulation onset (silent MLI)

OUTPUT DATA
    HISTOGRAM_PER_CELL__                     all ISIs for all MLI following unitary GC photostimulation
    ID_LIST__                                index for HISTOGRAM_PER_CELL__
"""
import os
import numpy as np
import scipy as sp
import pandas as pd
from matplotlib import pyplot as plt


DIR_CONTENT = os.listdir(r'D:/DATA/GC-MLI SeqPatch/MLI FIRING/INDIVIDUALS')

MED_ISI__ = []
CV__ = []
CV2__ = []
LVR__ = []
_R_ = 5
ALL = []
DISPLAYED___ = []
DISPLAYED_ISI = []
DISPLAYED_FREQ = []
PSTH = []
TOTAL_SPIKES = []
HISTOGRAM_PER_CELL__ = []
FIRING_PROBABILITY = []
COMBINED_HISTS = []
original_HISTOGRAM_PER_CELL__ = []
CV2_PER_CELL__ = []
ID_LIST__ = []

bin_size = 5
stimulation_delay = 500
cell_id_ = 0


for k in range(len(DIR_CONTENT)):
    A = np.array(pd.read_csv(r'D:/DATA\GC-MLI SeqPatch\MLI FIRING/INDIVIDUALS/'+str(DIR_CONTENT[k]), sep=' '), dtype=np.float64)
    ISI_LIST = []
    BASELINE_ISI_ = []
    CV2_LIST__ = []
    time_index_of_first_spike_latency = []
    HISTOGRAM___ = []
    SPIKING_OCCURENCE = np.linspace(0, 0, 100)
    for i in range(len(A)):
        LvR = []
        ISIn = []
        ISI_ = []
        
        HISTOGRAM__ = np.linspace(0, 0, 2000/bin_size)
        ALL.append(A[i])
        for j in range(len(A[i])-1):
            if np.isnan(A[i][j+1]) == False:
                ISI_.append(np.float(A[i][j+1]-A[i][j]))
                if 100 < A[i][j] < stimulation_delay:
                    BASELINE_ISI_.append(np.float(A[i][j+1]-A[i][j]))
                for m in range(len(HISTOGRAM__)):
                    if m*bin_size < A[i][j]< (m+1)*bin_size:
                        HISTOGRAM__[m] = HISTOGRAM__[m] + 1000/np.float(A[i][j+1]-A[i][j])
        for l in range(len(HISTOGRAM__)):
            if HISTOGRAM__[l] == 0 and l > 0:
                if l*bin_size:
                    HISTOGRAM__[l] = HISTOGRAM__[l-1]
        HISTOGRAM___ .append(HISTOGRAM__)

        #THIS IS TO ASSESS REGULARITY
        if len(ISI_) > 0:
            for l in range(len(ISI_)-1):
                ISIn.append(2*(abs(ISI_[l+1] - ISI_[l]))/(ISI_[l+1] + ISI_[l]))
                LvR.append((1-(4*ISI_[l]*ISI_[l+1])/ np.square((ISI_[l]+ISI_[l+1])))*(1+(4*_R_)/(ISI_[l]+ISI_[l+1])))
            ISI_LIST.append(ISI_)
            try:
                LvR = np.sum(LvR)*(3/(len(ISI_)-1))
            except:
                LvR = 0
            if k == 4:
                DISPLAYED_ISI.append(A[i])
            CV2_LIST__.append(ISIn)

        cumsum_ISI = np.cumsum(ISI_)
        for l in range(len(cumsum_ISI)):
            if cumsum_ISI[l] > stimulation_delay:
                try:
                    previous_spike_lag = stimulation_delay-cumsum_ISI[l-1]
                except:
                    previous_spike_lag = 0
                break

        for l in range(len(ISI_)):
            for m in range(100):
                if stimulation_delay+m*5-previous_spike_lag < cumsum_ISI[l] < stimulation_delay+m*5+5-previous_spike_lag:
                    SPIKING_OCCURENCE[m] = SPIKING_OCCURENCE[m] + 1

    '''           
    if cell_id_ != str(DIR_CONTENT[k]):
        plt.figure(figsize=(4,4))
    plt.plot(np.linspace(0,2000,len(HISTOGRAM___[0])), np.mean(HISTOGRAM___, axis=0),c='black', alpha=0.9)
    '''
    
    mean = np.nanmean(HISTOGRAM___, axis=0)
    if True:
        HISTOGRAM_PER_CELL__.append(np.nan_to_num(mean))
        COMBINED_HISTS.append(HISTOGRAM___)
    else:
        HISTOGRAM_PER_CELL__.append(mean)

    original_HISTOGRAM_PER_CELL__ .append(np.nan_to_num(mean))
    CV2_PER_CELL__.append(np.nanmean(np.array(np.concatenate(CV2_LIST__), dtype=np.float64)))
    FIRING_PROBABILITY.append(SPIKING_OCCURENCE/len( HISTOGRAM___))

    cell_id_ = str(DIR_CONTENT[k])
    ID_LIST__.append(cell_id_)

    try:
        temp__ = [ISI_LIST[i][0:10] for i in range(len(ISI_LIST))]
        TOTAL_SPIKES.append(np.nanmedian([len(ISI_LIST[i]) for i in range(len(ISI_LIST))]))
        DISPLAYED_FREQ.append(np.mean(ISI_LIST[i][0:10], axis=0) for i in range(len(ISI_LIST)))

        if k == 4:
            DISPLAYED___.append(temp__)
        else:
            temp__ = np.concatenate(temp__)
        MED_ISI__.append(np.nanmedian(temp__))
        CV__.append(np.nanstd(temp__)/np.nanmean(temp__))
        CV2__.append(np.nanmean(ISIn))
        LVR__.append(LvR)
    except:
        pass

#This figure compares WT and Syn2KO instantaneous frequency after unitary GC photostimulation
plt.figure(figsize=(3, 3))
mean = np.mean([HISTOGRAM_PER_CELL__[i] for i in range(len(ID_LIST__)) if ('S2KO' in ID_LIST__[i]) == False], axis=0)
sem = sp.stats.sem([HISTOGRAM_PER_CELL__[i] for i in range(len(ID_LIST__)) if ('S2KO' in ID_LIST__[i]) == False], axis=0)
plt.plot(np.linspace(-500, 1500, len(mean)), mean)
plt.fill_between(np.linspace(-500, 1500, len(mean)), mean+sem, mean-sem, alpha=0.3)
LIST = [HISTOGRAM_PER_CELL__[i] for i in range(len(ID_LIST__)) if ('S2KO' in ID_LIST__[i])==True]
mean = np.mean(LIST, axis=0)
sem = sp.stats.sem([HISTOGRAM_PER_CELL__[i] for i in range(len(ID_LIST__)) if ('S2KO' in ID_LIST__[i])==True], axis=0)
plt.plot(np.linspace(-500, 1500, len(mean)), mean, color='red')
plt.fill_between(np.linspace(-500, 1500, len(mean)), mean+sem, mean-sem, alpha=0.3, color='red')
plt.xlim(0, 300)
plt.ylim(0,)
plt.xlabel('Time post-stimulus (ms)')
plt.ylabel('Frequency(Hz)')
plt.tight_layout()


plt.figure(figsize=(3, 3))
mean = np.nanmean([FIRING_PROBABILITY[i] for i in range(len(ID_LIST__)) if ('S2KO' in ID_LIST__[i]) == False and ('2017' in ID_LIST__[i]) == True], axis=0)
sem = sp.stats.sem(np.nan_to_num(np.array([FIRING_PROBABILITY[i] for i in range(len(ID_LIST__)) if ('S2KO' in ID_LIST__[i]) == False and ('2017' in ID_LIST__[i]) == True])), axis=0)
plt.plot(np.linspace(0, len(mean)*5, len(mean)), mean)
plt.fill_between(np.linspace(0, len(mean)*5, len(mean)), mean+sem, mean-sem, alpha=0.3)


LIST = [FIRING_PROBABILITY[i] for i in range(len(ID_LIST__)) if ('S2KO' in ID_LIST__[i]) == True]
mean = np.nanmean(LIST, axis=0)
sem = sp.stats.sem(np.nan_to_num(np.array([FIRING_PROBABILITY[i] for i in range(len(ID_LIST__)) if ('S2KO' in ID_LIST__[i]) == True])), axis=0)
plt.plot(np.linspace(0, len(mean)*5, len(mean)), mean, color='red')
plt.fill_between(np.linspace(0, len(mean)*5, len(mean)), mean+sem, mean-sem, alpha=0.3, color='red')

plt.xlim(0, 150)
plt.ylim(0,)
plt.xlabel("Time post-stimulus (ms)")
plt.ylabel("Firing Probability")
plt.tight_layout()


'''
plt.figure(figsize=(3, 3))
a, a_ = np.histogram([time_index_of_max_frequency[i] for i in range(len(ID_LIST__)) if ('S2KO' in ID_LIST__[i]) == False], bins=20)
b, b_ = np.histogram([time_index_of_max_frequency[i] for i in range(len(ID_LIST__)) if ('S2KO' in ID_LIST__[i]) == True], bins=20)
plt.plot(a_[1:len(a_)]*5,np.cumsum(a)/np.sum(a))
plt.plot(b_[1:len(b_)]*5,np.cumsum(b)/np.sum(b), color='red')
plt.xlabel("Frequency Peak Delay (ms)")
plt.ylabel("Cumulative Fraction")
plt.xlim(0, 250)
'''