# -*- coding: utf-8 -*-
"""
@author: KEVIN-DORGANS

>> generate downsampled and averaged *.csv files based on WCP voltagtage clamp data from unitary GC photostimulation

/!\ RUNS ONLY IN SynaptiQs with SQL database

http://synaptiqs.wixsite.com/synaptiqs/setup

"""

#DETECTION AND FILE CREATION

import numpy as np
import scipy as sp
import pandas as pd
from matplotlib import pyplot as plt


tagged_sweeps = Requete.tag['Selection']

bin_size = 5

HISTOGRAM__PER__CELL = []
plt.figure(figsize=(4, 4))
for i in range(len(tagged_sweeps)):
    if tagged_sweeps[int(i)] > 0:
		id = str(Requete.Analogsignal_ids[i])
		Navigate.Load_This_Trace(id)
		signal_name = str(Requete.Block_Info[i])
		sampling = Requete.timescale[1]
		RAW_SIGNAL_ =  np.array(Navigate.si)

		HISTOGRAM__ = np.array(np.linspace(0,0,2000/bin_size))
		for j in np.linspace(0,2000/bin_size-1,2000/bin_size):
			HISTOGRAM__[j]=np.trapz(RAW_SIGNAL_[j*bin_size/sampling:j*bin_size/sampling+(bin_size/sampling)])
	try:
		HISTOGRAM__PER__CELL.append(HISTOGRAM__)
		plt.plot(HISTOGRAM__, c='black', alpha=0.6)
	except:
		pass
	
HISTOGRAM__PER__CELL = np.array(HISTOGRAM__PER__CELL)
DATA_f = pd.DataFrame(data=HISTOGRAM__PER__CELL)
FILE_NAME =  'F:\DATA\GC-MLI SeqPatch\MLI PATCH\INDIVIDUALS_2\ '+str(signal_name) +' MLI_firing.csv'
DATA_f.to_csv(FILE_NAME, header=None, sep=' ', index=None)

plt.plot(np.std(HISTOGRAM__PER__CELL, axis=0), linewidth=3)
plt.tight_layout()
