# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 18:47:27 2018
@author: Kevin.Dorgans

>> loads individual GC-MLI STP observations from *.csv files
>> sorts observations between WT and S2KO observations
>> labels observations with specific tags to retrieve original *.csv and *.wcp file
    ex. '2017-02-20 60S3C1Z1 Charge_temp_.csv'
        'YYYY-MM-DD CELL_ID, SLICE_ID, CELL_ID, ZONE_ID'
        
>> normalization using Vector Space Model     sklearn.preprocessing.normalize
>> computes PCA transformation                sklearn.decomposition.PCA
>> KMean clustering                           sklearn.cluster.vq.kmeans
        
LIST OF OUTPUT FIGURES
    WT PCA plot 
    S2KO PCA plot
    STP plot for eEPSC1-10 by KMean category
    Bar plot for eEPSC1
    Bar plot for PPR
    STP plot for normalized WT values
    Elbow plot for KMeans
    Pie chart for WT vs KO cluster identity
    STP plot for normalized S2KO values

IMPORTANT VARIABLES
    Experiments__                 full WT experimental data per synapse
    Experiment__Index__           index tag for WT data
    Experiments__S2KO             full Synapsin2 Knockout data
    Experiment__Index__S2KO       index for S2KO data
    center_normed_matrix          normalized matrix of WT data
    eigenvalue_WT_synapses        PC transformation of WT data
    center_normed_matrix_S2KO     normalized matrix of S2KO data
    eigenvalue_WT_synapses_S2KO   PC transformation of S2KO data
    KMC_clusters                  cluster identity for WT data
    KMC_clusters_S2KO             cluster identity for S2KO data
"""



import pandas as pd
import numpy as np
import scipy as sp
import os
import sklearn.decomposition
import sklearn.preprocessing
import sklearn.cluster
from scipy.cluster.vq import kmeans
from scipy.spatial.distance import cdist, pdist
from matplotlib import pyplot as plt
from matplotlib import cm



N_CLUST = 4 #Cluster number estimated from KMean elbows
clr = cm.viridis(np.linspace(0.05, 0.85, N_CLUST)).tolist()

#This function loads and sorts individual synapses from *.csv files between WT (Experiments__) ans S2KO (Exper)iments_S2KO
def _sort_experiments():
    path = os.listdir(r'D:\DATA\GC-MLI MINIMAL\GC-MLI MINIMAL 2019')
    Experiments__ = []
    Experiment__Index__ = []
    Experiments__S2KO = []
    Experiment__Index__S2KO = []

    for k in range(len(path)):
        full_data_matrix = np.array(pd.read_csv(r'D:\DATA\GC-MLI MINIMAL\GC-MLI MINIMAL 2019' +"/"+str(path[k]), sep=' ', header=None))       

        if ('S2KO' in str(path[k])) == False:
            Experiments__.append(full_data_matrix[0])
            Experiment__Index__.append(path[k])
            print(str(path[k]))
        elif  ('S2KO' in str(path[k])) == True:
            Experiments__S2KO.append(full_data_matrix[0])
            Experiment__Index__S2KO.append(path[k])
    return(Experiments__, Experiment__Index__, Experiments__S2KO, Experiment__Index__S2KO)


Experiments__, Experiment__Index__, Experiments__S2KO, Experiment__Index__S2KO = _sort_experiments()

Experiments__ = np.array(np.transpose(Experiments__)[1:11].transpose(), dtype=np.float64)
Experiments__S2KO = np.array(np.transpose(Experiments__S2KO)[1:11].transpose(), dtype=np.float64)

PCA = sklearn.decomposition.pca.PCA(n_components=5)
KMC = sklearn.cluster.KMeans(n_clusters=N_CLUST)

IMPUTER = sklearn.preprocessing.Imputer(strategy="median")
center_normed_matrix = sklearn.preprocessing.normalize(IMPUTER.fit_transform(Experiments__))
center_normed_matrix_S2KO = sklearn.preprocessing.normalize(IMPUTER.fit_transform(Experiments__S2KO))


PCA_fit = PCA.fit(center_normed_matrix)
eigenvalue_WT_synapses = PCA_fit.transform(center_normed_matrix)
eigenvalue_S2KO_synapses = PCA_fit.transform(center_normed_matrix_S2KO)
KMC_clusters = KMC.fit_predict(eigenvalue_WT_synapses)
KMC_clusters_S2KO = KMC.fit(eigenvalue_WT_synapses).predict(eigenvalue_S2KO_synapses)


'''
#Figure testing data import
MEAN = np.mean(Experiments__, axis=0)
SEM = sp.stats.sem(Experiments__, axis=0)
plt.figure(figsize=(3, 3))
plt.plot(np.linspace(1, 10, len(MEAN), len(MEAN)), MEAN)
plt.fill_between(np.linspace(1, 10, len(MEAN), len(MEAN)), MEAN+SEM, MEAN-SEM, alpha=0.2)
'''

#WT PCA transformed data plot components 1, 2
plt.figure(figsize=(3, 3))
for i in range(len(eigenvalue_WT_synapses)):
    plt.scatter(eigenvalue_WT_synapses[i][0], eigenvalue_WT_synapses[i][1], color=clr[KMC_clusters[i]])
    plt.xlabel('PC1')
    plt.ylabel('PC2')
plt.tight_layout()

#S2KO PCA transformed data plot components 1, 2
plt.figure(figsize=(3, 3))
for i in range(len(eigenvalue_WT_synapses)):
    plt.scatter(eigenvalue_WT_synapses[i][0], eigenvalue_WT_synapses[i][1], color=clr[KMC_clusters[i]], alpha=0.2)
for i in range(len(eigenvalue_S2KO_synapses)):
    plt.scatter(eigenvalue_S2KO_synapses[i][0], eigenvalue_S2KO_synapses[i][1], color='red')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
plt.tight_layout()

#Short-term-plasticity plot for eEPSC1-10 by KMean category
plt.figure(figsize=(4, 3))
for i in range(N_CLUST):
    plt.plot(np.linspace(1, 10, 10), np.mean([center_normed_matrix[j]/center_normed_matrix[j][0] for j in range(len(eigenvalue_WT_synapses)) if KMC_clusters[j] == i], axis=0), color=clr[i])
    plt.scatter(np.linspace(1, 10, 10), np.mean([center_normed_matrix[j]/center_normed_matrix[j][0] for j in range(len(eigenvalue_WT_synapses)) if KMC_clusters[j] == i], axis=0), color=clr[i])
    MEAN = np.mean([center_normed_matrix[j]/center_normed_matrix[j][0] for j in range(len(eigenvalue_WT_synapses)) if KMC_clusters[j] == i], axis=0)
    SEM = sp.stats.sem([center_normed_matrix[j]/center_normed_matrix[j][0] for j in range(len(eigenvalue_WT_synapses)) if KMC_clusters[j] == i], axis=0)
    plt.fill_between(np.linspace(1, 10, 10), MEAN+SEM, MEAN-SEM, alpha=0.3, color=clr[i])
plt.xlabel('eEPSC#')
plt.ylabel('eEPSCn/eEPSC1')
plt.tight_layout()

#Bar plot for eEPSC amplitude
plt.figure(figsize=(2, 3))
for i in range(N_CLUST):
    temp__ = [Experiments__[j][0] for j in range(len(eigenvalue_WT_synapses)) if KMC_clusters[j] == i]
    plt.scatter(np.linspace(i, i, len(temp__)), temp__, color=clr[i], alpha=0.2)
    plt.bar(i, np.mean(temp__), color=clr[i], alpha=0.3)
plt.xlim(-1, 4)
plt.ylabel("eEPSC1 Amplitude (pA)")
plt.xlabel("Category")
plt.tight_layout()

#Bar plot for PPR
plt.figure(figsize=(2, 3))
for i in range(N_CLUST):
    temp__ = [Experiments__[j][1]/Experiments__[j][0] for j in range(len(eigenvalue_WT_synapses)) if KMC_clusters[j] == i]
    plt.scatter(np.linspace(i, i, len(temp__)), temp__, color=clr[i], alpha=0.2)
    plt.bar(i, np.mean(temp__), color=clr[i], alpha=0.3)
    print(np.mean(temp__))
plt.xlim(-1, 4)
plt.ylim(0, 4)
plt.ylabel("eEPSC2/eEPSC1")
plt.xlabel("Category")
plt.tight_layout()

#STP plot for normalized WT values
plt.figure(figsize=(N_CLUST*2.5, 2.5))
for i in range(N_CLUST):
    plt.subplot(100+10*N_CLUST+i+1)
    temp__ = np.mean([center_normed_matrix[j] for j in range(len(eigenvalue_WT_synapses)) if KMC_clusters[j] == i], axis=0)
    temp__2 = sp.stats.sem([center_normed_matrix[j] for j in range(len(eigenvalue_WT_synapses)) if KMC_clusters[j] == i], axis=0)

    plt.scatter(np.linspace(1, 10, 10), temp__, color=clr[i])
    plt.plot(np.linspace(1, 10, 10), temp__, color=clr[i])
    plt.fill_between(np.linspace(1, 10, 10), temp__-temp__2, temp__+temp__2, color=clr[i], alpha=0.4)
    plt.xlabel('eEPSC#')
    plt.ylabel('Normalized Release')
    plt.ylim(0, 1)
    plt.xlim(0, 10)
plt.tight_layout()

if True:
    K_MAX = int(8)
    KK = range(1, K_MAX+1)
    
    kIdx = 1
    KM = [kmeans(eigenvalue_WT_synapses, k) for k in KK]
    centroids = [cent for (cent, var) in KM]
    D_k = [cdist(eigenvalue_WT_synapses, cent, 'euclidean') for cent in centroids]
    cIdx = [np.argmin(eigenvalue_WT_synapses, axis=1) for eigenvalue_WT_synapses in D_k]
    dist = [np.min(eigenvalue_WT_synapses, axis=1) for eigenvalue_WT_synapses in D_k]

    tot_withinss = [sum(d**2) for d in dist]  # Total within-cluster sum of squares
    totss = sum(pdist(eigenvalue_WT_synapses)**2)/eigenvalue_WT_synapses.shape[0]       # The total sum of squares
    betweenss = totss - tot_withinss          # The between-cluster sum of squares

    #Elbow plot for KMeans
    fig = plt.figure(figsize=(3, 3))
    ax = fig.add_subplot(111)
    ax.plot(KK, betweenss/totss*100)
    ax.scatter(KK, betweenss/totss*100, s=8)
    #ax.plot(KK[kIdx]+0.4, betweenss[kIdx]/totss*100, marker='o', markersize=12,
    #markeredgewidth=2, markeredgecolor='r', markerfacecolor='None')
    ax.set_ylim((0, 100))
    ax.scatter(KK[3], np.array(betweenss/totss*100)[3], color='red')
    ax.plot((4, 4), (0, 100), lw=1, c='red')
    plt.grid(True)
    plt.xlabel('Number of clusters')
    plt.ylabel('Percentage of variance explained (%)')
    plt.title('Elbow for KMeans clustering')
    plt.tight_layout()
    plt.show()

#Pie chart for WT vs KO cluster identity
fig = plt.figure(figsize=(6, 3))
ax = fig.add_subplot(121)
ax.pie([KMC_clusters.tolist().count(i) for i in range(N_CLUST)], labels=(range(N_CLUST)), colors=clr, autopct='%1.1f%%')
ax.set_title('WT PROFILES')
ax = fig.add_subplot(122)
ax.pie([KMC_clusters_S2KO.tolist().count(i) for i in range(N_CLUST)], labels=(range(N_CLUST)), colors=clr, autopct='%1.1f%%')
ax.set_title('Syn2-KO PROFILES')
plt.tight_layout()

#STP plot for normalized S2KO values
plt.figure(figsize=(N_CLUST*2.5, 2.5))
for j in range(N_CLUST):
    plt.subplot(100+10*N_CLUST+j+1)
    temp__ = [center_normed_matrix_S2KO[i] for i in range(len(center_normed_matrix_S2KO)) if KMC_clusters_S2KO[i] == j]
    temp__2 = sp.stats.sem([center_normed_matrix_S2KO[i] for i in range(len(center_normed_matrix_S2KO)) if KMC_clusters_S2KO[i] == j])

    plt.plot(np.linspace(1, 10, 10), np.mean(temp__, axis=0), color='red')
    plt.scatter(np.linspace(1, 10, 10), np.mean(temp__, axis=0), color='red')
    plt.fill_between(np.linspace(1, 10, 10), np.mean(temp__, axis=0)+temp__2, np.mean(temp__, axis=0)-temp__2, color='red', alpha=0.2)
    plt.xlabel('eEPSC#')
    plt.ylabel('Normalized Release')
    plt.ylim(0,1)
plt.tight_layout()
