# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 16:12:11 2017
@author: Kevin.Dorgans


/!\ REQUIRES RUNNING:
    1 > GC_photostimulation__MLI_FIRING_ANALYSIS.py 
    2 > GC_photostimulation_GC-MLI_CHARGE_ANALYSIS.py 
    
>> computes vector normalization > PCA > KMeans onto FULL GC-MLI eEPSC CHARGE DATASET
>> classifies GC-MLI eEPSC CHARGE DATA in class depending on STP profile
>> sorts WT versus Syn2KO (red) observations
>> retrieves MLI FIRING from experiments WHERE MLI firing recording is SUCESSFULLY COUPLED WITH MLI WCP Vclamp
>> COMPARES MLI PEAK FIRING DELAY to eEPSC SHORT TERM PLASTICITY PROFILE
"""

import sklearn
import numpy as np
import scipy as sp
from sklearn import preprocessing as pre
from statsmodels import robust
from matplotlib import cm
from matplotlib import pyplot as plt

N_Clust = 3
clr = cm.viridis(np.linspace(0.05, 0.85, N_Clust)).tolist()
PCA = sklearn.decomposition.PCA(n_components=5)
KMeans = sklearn.cluster.KMeans(n_clusters=N_Clust)

Shift_List = []
Shifted_HISTOGRAM = []

for i in range(len(eEPSC_HISTOGRAM)):
    j = 100
    NOIZE = np.std(eEPSC_HISTOGRAM[i][0:100])
    OFFSET = np.median(eEPSC_HISTOGRAM[i][0:100])
    while eEPSC_HISTOGRAM[i][j]-OFFSET > -3*NOIZE and i < 110:
        j += 1
    Shift_List.append(j)
    Shift_index = j-2
    Shifted_HISTOGRAM.append(eEPSC_HISTOGRAM[i][Shift_index:Shift_index+60])

NewSIGNAL_NAMEs = [eEPSC_HISTOGRAM_ID_LIST__[i] for i in range(len(Shifted_HISTOGRAM)) if ('S2KO' in eEPSC_HISTOGRAM_ID_LIST__[i]) == False]
NewSIGNAL_NAMEs_KO = [eEPSC_HISTOGRAM_ID_LIST__[i] for i in range(len(eEPSC_HISTOGRAM_ID_LIST__)) if ('S2KO' in eEPSC_HISTOGRAM_ID_LIST__[i]) == True]


PhotostimEPSCs = [Shifted_HISTOGRAM[i] for i in range(len(Shifted_HISTOGRAM)) if ('S2KO' in eEPSC_HISTOGRAM_ID_LIST__[i]) == False]

ToFIT = pre.normalize([PhotostimEPSCs[i][2:20] for i in range(len(PhotostimEPSCs))])
PCA_fit = PCA.fit(ToFIT)
X = PCA_fit.transform(ToFIT)
Y = PCA_fit.transform(pre.normalize([Shifted_HISTOGRAM[i][0:18] for i in range(len(Shifted_HISTOGRAM)) if ('S2KO' in eEPSC_HISTOGRAM_ID_LIST__[i]) == True]))



KMeansPredict = KMeans.fit_predict(ToFIT)

#Displays PC1 versus PC2 for WT DATA
plt.figure(figsize=(3, 3))
for i in range(len(X)):
    plt.scatter(X[:, 0][i], X[:, 1][i], c=clr[KMeansPredict[i]])
plt.xlabel('PC1 '+str(PCA.explained_variance_ratio_[0]*100))
plt.ylabel('PC2 '+str(PCA.explained_variance_ratio_[1]*100))
plt.tight_layout()

if True:
    #Displays PCA PARAMETERS FOR WT Data
    fig = plt.figure(figsize=(8, 2.5))
    ax = fig.add_subplot(131)
    for i in range(len(PCA.explained_variance_)):
        ax.bar(i, PCA.explained_variance_ratio_[i], color='black', alpha=0.7-i/10)
    plt.ylabel('% EXPLAINED VARIANCE')
    plt.xlabel('PRINCIPAL COMPONENTS')

    ax = fig.add_subplot(132)
    for i in range(len(PCA.components_[0])):
        ax.bar(i*5, PCA.components_[0][i], color='black')
    plt.ylabel('CONTRIBUTION IN PC1')
    plt.xlabel('Time post-stimulus (ms)')

    ax = fig.add_subplot(133)
    for i in range(len(PCA.components_[1])):
        plt.bar(i*5, PCA.components_[1][i], color='black')
    plt.ylabel('CONTRIBUTION IN PC2')
    plt.xlabel('Time post-stimulus (ms)')
    plt.tight_layout()
    plt.show()


#DISPLAYS GC-MLI eEPSC CHARGE FOLLOWING PHOTOSTIMULATION
plt.figure(figsize=(3, 3))
for i in range(N_Clust):
    THE_LIST = [PhotostimEPSCs[j]/5 for j in range(len(PhotostimEPSCs)) if KMeansPredict[j] == i]
    plt.plot(np.linspace(0, 300, len(PhotostimEPSCs[i])), np.mean(THE_LIST, axis=0), c=clr[i])
    THE_LIST_sem = sp.stats.sem(THE_LIST, axis=0)
    plt.fill_between(np.linspace(0, 300,len(PhotostimEPSCs[i])), np.mean(THE_LIST, axis=0), np.mean(THE_LIST, axis=0)-THE_LIST_sem, color=clr[i], alpha=0.3)

    DELAY__ = [PhotostimEPSCs[j].tolist().index(np.min(PhotostimEPSCs[j]))*300/len(PhotostimEPSCs[i]) for j in range(len(PhotostimEPSCs)) if KMeansPredict[j] == i]

    print("C'"+str(i+1)+ ": (n="+str(len(DELAY__))+") "+str(np.mean([np.min(THE_LIST[k]) for k in range(len(THE_LIST))])) + " +/- "+ str(sp.stats.sem([np.min(THE_LIST[k]) for k in range(len(THE_LIST))])))
    print("Delay: "+str(np.mean(DELAY__))+" ms +/-  "+str(sp.stats.sem(DELAY__)))

plt.xlim(0, 300)
plt.xlabel('Time post-stimulus (ms)')
plt.ylabel('eEPSC Charge (fC)')
plt.tight_layout()

DDG = scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(ToFIT, method='ward'), no_plot='True')
DDG = np.array(DDG['leaves'], dtype=np.int)
Normed_PhotoS = pre.normalize(PhotostimEPSCs)
MatS = [Normed_PhotoS[DDG[i]] for i in range(len(PhotostimEPSCs))]
plt.set_cmap('viridis_r')
#plt.matshow(MatS)

#========= NECESSITE ANALYSE DE FIRING + ANALYSE DE PATCH =====

To_Plot = []
Firing_IDs__ = []
Peak_Index_ = []
Rise_rime_10_90_ = []
eEPSC_DELAY_ = []

plt.figure(figsize=(5, 3))
for i in range(N_Clust):
    To_Plot.append([])
    Firing_IDs__.append([])
    Peak_Index_.append([])
    Rise_rime_10_90_.append([])
    eEPSC_DELAY_.append([])
    
for k in range(N_Clust):
    for i in range(len(NewSIGNAL_NAMEs)):
        if NewSIGNAL_NAMEs[i] in ID_LIST__:
            if KMeansPredict[i] == k:
                To_Append = sp.signal.medfilt(HISTOGRAM_PER_CELL__[ID_LIST__.index(NewSIGNAL_NAMEs[i])], 5)
                #plt.plot(np.linspace(-500,1500,len(To_Append)), To_Append,c=clr[k], alpha=0.3)
                if True:             
                    To_Plot[k].append((To_Append/np.max(To_Append[100::])))
                    max_index_ = To_Append[100::].tolist().index(np.max(To_Append[100::]))*5
                    Peak_Index_[k].append(max_index_)
                l = 101
                SD_NOIZE = np.std(HISTOGRAM_PER_CELL__[ID_LIST__.index(NewSIGNAL_NAMEs[i])][5:100])
                BASELINE_FIRING = np.median(HISTOGRAM_PER_CELL__[ID_LIST__.index(NewSIGNAL_NAMEs[i])][50:100])

                while To_Append[l] < BASELINE_FIRING:
                    l += 1
                    if l > 250:
                        break
                
                while To_Append[l] < np.max(To_Append[100:150])-(SD_NOIZE):
                    l += 1

                if (l-100)*5 > 0:
                    Rise_rime_10_90_[k].append((l-100)*5)
                    Firing_IDs__[k].append(NewSIGNAL_NAMEs[i])
                    l = 0
                    while PhotostimEPSCs[i][l] > np.min(PhotostimEPSCs[i])+abs(np.nanmean(PhotostimEPSCs[i][0:2])):
                        l += 1
                    temp__ = (PhotostimEPSCs[i].tolist().index(PhotostimEPSCs[i][l]))*300/len(PhotostimEPSCs[i])
                    eEPSC_DELAY_[k].append(temp__)
    SEM = sp.stats.sem(To_Plot[k], axis=0)
    MEAN = np.mean(To_Plot[k], axis=0)
    #for j in range(len(To_Plot[k])):
        #plt.plot(np.linspace(-500,1500,len(To_Plot[k][j])), To_Plot[k][j], c=clr[k])
    plt.plot(np.linspace(-500, 1500, len(MEAN)), MEAN, c=clr[k])
    plt.fill_between(np.linspace(-500, 1500, len(np.mean(To_Plot[k], axis=0))), MEAN+SEM, MEAN, alpha=0.3, color=clr[k])
plt.tight_layout()

#PeakIndexFor KO
Peak_Index_KO_ = []


Rise_rime_10_90_KO = []
Rise_rime_10_90_KO_GROUPED = []
PC1_KO_RETRIEVED = []
for i in [ID_LIST__[i] for i in range(len(Shifted_HISTOGRAM)) if ('S2KO' in ID_LIST__[i]) == True]:
    if i in [eEPSC_HISTOGRAM_ID_LIST__[i] for i in range(len(Shifted_HISTOGRAM)) if ('S2KO' in SIGNAL_NAMES[i]) == True]:
        j = ID_LIST__.index(i)
        Peak_Index_KO_.append(HISTOGRAM_PER_CELL__[j][100:170].tolist().index(np.max(HISTOGRAM_PER_CELL__[j][100:170]))*5)
        PC1_KO_RETRIEVED.append(Y[[ID_LIST__[i] for i in range(len(Shifted_HISTOGRAM)) if ('S2KO' in ID_LIST__[i]) == True].index(i)][0])
        To_Append = sp.signal.medfilt(HISTOGRAM_PER_CELL__[j], 3)

        SD_NOIZE = np.std(HISTOGRAM_PER_CELL__[j][0:100])
        l = 99
        while To_Append[l] < 3*SD_NOIZE:
            l += 1
            if l > 200:
                break
        m = l
        while To_Append[m] >= np.max(To_Append[100:150])-3*SD_NOIZE:
            m += 1
            if m > 170:
                break
        Rise_rime_10_90_KO.append((l+(m-l)/2-100)*5)

        _TEMP_ = []
        for k in range(len(COMBINED_HISTS[j])):
            l = 100
            SD_NOIZE = np.nanstd(COMBINED_HISTS[j][k][50:100])
            while COMBINED_HISTS[j][k][l] < 3*SD_NOIZE:
                l += 1
                if l > 170:
                    break
            _TEMP_.append(COMBINED_HISTS[j][k].tolist().index(np.max(COMBINED_HISTS[j][k][100:170])))
        Rise_rime_10_90_KO_GROUPED.append(_TEMP_)

#plt.ylim(0,1.1)
plt.xlim(0, 300)
plt.xlabel('Time post-stimulus (ms)')
plt.ylabel('Firing frequency acceleration (%)')


plt.figure(figsize=(1, 3))
for j in range(3):
    SAMPLE_a_ = Rise_rime_10_90_[j]

    for k in range(len(SAMPLE_a_)):
        plt.scatter(j,SAMPLE_a_[k],  color=clr[j], s=30)
    plt.bar(j,np.mean(SAMPLE_a_), alpha=0.3,  color=clr[j])
plt.xlim(-1, 3)
plt.ylim(0,)
plt.ylabel('Delay to Frequency Peak(ms)')
plt.tight_layout()


a_fit = []
b_fit = []
plt.figure(figsize=(2, 3))
for j in range(len(Firing_IDs__)):
    for i in range(len(Firing_IDs__[j])):
        a = X[np.array(NewSIGNAL_NAMEs).tolist().index(Firing_IDs__[j][i])][0]

        a_fit.append(a)
        b_fit.append(Rise_rime_10_90_[j][i])
        plt.scatter(a, Rise_rime_10_90_[j][i], c=clr[j], s=15)

fit = np.polyfit(a_fit, b_fit, 1)
fit_fn = np.poly1d(fit) 


plt.plot(np.linspace(-0.5, 0.5, 10), fit_fn(np.linspace(-0.5, 0.5, 10)), c='black', ls='--', alpha=0.3)
#plt.scatter(PC1_KO_RETRIEVED, Peak_Index_KO_ , color='red')
plt.xlabel('eEPSCs PC1')
plt.ylabel('Delay to Firing Peak (ms)')
plt.text(0, 0,str(sp.stats.pearsonr(a_fit, b_fit)))
plt.xlim(-0.5, 0.5)
plt.ylim(0, 200)
plt.tight_layout()


a_fit = []
b_fit = []
plt.figure(figsize=(2, 3))
for j in range(len(Firing_IDs__)):
    for i in range(len(Firing_IDs__[j])):
        a_fit.append(eEPSC_DELAY_[j][i])
        b_fit.append(Rise_rime_10_90_[j][i])
        plt.scatter(eEPSC_DELAY_[j][i],Rise_rime_10_90_[j][i], c=clr[j], s=15)

fit = np.polyfit(a_fit,b_fit,1)
fit_fn = np.poly1d(fit)
plt.xlabel('eEPSC delay (ms)')
plt.tight_layout()

plt.plot(np.linspace(-0.5, 0.5, 10), fit_fn(np.linspace(-0.5, 0.5, 10)), c='black', ls='--', alpha=0.3)
#plt.scatter(PC1_KO_RETRIEVED, Peak_Index_KO_ , color='red')
plt.xlabel('Delay to eEPSC Peak (ms)')
plt.ylabel('Delay to Firing Peak (ms)')
plt.text(0, 0, str(sp.stats.pearsonr(a_fit, b_fit)))
plt.ylim(0, 200)
plt.tight_layout()



#PCA map for SYN2KO
plt.figure(figsize=(3, 3))
for i in range(len(X)):
    plt.scatter(X[:,0][i], X[:,1][i], c=clr[KMeansPredict[i]], alpha=0.3)

for i in range(len(Y)):
    plt.scatter(Y[:,0][i], Y[:,1][i], c='red')
plt.tight_layout()

plt.figure(figsize=(3, 3))
for i in range(N_Clust):
    THE_LIST = pre.normalize([PhotostimEPSCs[j]/5 for j in range(len(PhotostimEPSCs)) if KMeansPredict[j] == i])
    plt.plot(np.linspace(0, 300, len(PhotostimEPSCs[i])), np.mean(THE_LIST, axis=0), c=clr[i])
    THE_LIST_sem = sp.stats.sem(THE_LIST, axis=0)
    plt.fill_between(np.linspace(0, 300, len(PhotostimEPSCs[i])), np.mean(THE_LIST, axis=0), np.mean(THE_LIST, axis=0)-THE_LIST_sem, color=clr[i], alpha=0.3)

HistKO = pre.normalize([Shifted_HISTOGRAM[i][0:20]/5 for i in range(len(Shifted_HISTOGRAM)) if ('S2KO' in eEPSC_HISTOGRAM_ID_LIST__[i]) == True ])
plt.plot(np.linspace(0, 150, len(HistKO[0])), np.mean(HistKO, axis=0), color='red')
THE_LIST_sem = sp.stats.sem(HistKO, axis=0)
plt.fill_between(np.linspace(0, 150, len(HistKO[0])), np.mean(HistKO, axis=0), np.mean(HistKO, axis=0)-THE_LIST_sem, color='red', alpha=0.3)
plt.xlim(0, 140)
plt.xlabel('Time post-stimulus (ms)')
plt.ylabel('Normalized Release Rate (fC/ms)')
plt.tight_layout()
