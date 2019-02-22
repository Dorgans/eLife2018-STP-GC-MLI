# -*- coding: utf-8 -*-

import easygui
from neo import io

def load_WCPWave__(): 
    RAW = []
    directory = easygui.fileopenbox(default=r'D:\DATA\RAW')
    r = io.WinWcpIO(filename=directory)
    block = r.read_block()
    for i in range(len(block.segments)):
        RAW.append(block.segments[0].analogsignals[0])
    sampling__ = float(block.segments[0].analogsignals[0].sampling_rate)
    return RAW, 1000/sampling__


def Simple_eEPSC_Train_Analysis():
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib import pyplot as plt
    import numpy as np
    from scipy.optimize import curve_fit
    import pandas as pd
    import seaborn as sns
    import time
    import os

    def linear_func(x, a, b):
        y = a*x +b
        return y

    with PdfPages(str(os.getcwd())+r'\output.pdf') as pdf:
        __TIMESTART__ = time.time()

        AMP_VALUES = []
        CORRECTED_AMP_VALUES = []
        CHARGE_VALUES = []
        CORRECTED_CHARGE_VALUES = []
        SYN_name = []
        OFFSET = float
        sdNOISE = float
        FULL_EPSC_TRAIN_ = []
        FULL_PEAK_OCCURENCE_ = []
        FULL_RAW_TRACES = []
        FULL_CHARGE_NOISE_LIST = []
        FULL_AMP_NOISE_LIST = []
        popt = []
        tOFFSET = []
        occOFFSET = []
        _EPSC_NO_OFFSET = []
        EXTRACTED_CURRENT = []
        peak_occurence = []
        EXTRACTED_CURRENT = []
        event_latency__ = []
        
        ___anal_freq = 10
        ___delay = 500
        ___window = 8.3
        ___current_delay = 1.1
        ___number_of_pulses = 10

        RAW_WCP, sampling = load_WCPWave__()
        
        signal_name = __TIMESTART__
        print(str("ANALYZING "+str(signal_name)))
        print("----------------------------")
        print(str("SAMPLING IS "+str(sampling)))
        RAW_SIGNAL = np.concatenate(np.mean(RAW_WCP, axis=0))

        #This will estimate noise sd and charge related to time window
        EXTRACTED_CURRENT.append(RAW_SIGNAL[0:int(___window/sampling)])
        OFFSET = np.median(EXTRACTED_CURRENT)
        RAW_SIGNAL = RAW_SIGNAL - OFFSET        
        sdNOISE = np.std(np.concatenate(EXTRACTED_CURRENT, axis=0))
        chNOISE = abs(np.trapz(np.subtract(EXTRACTED_CURRENT, OFFSET))*sampling)

        #This estimates linear fit if current offset @eEPSC(n+1) is huge 
        for j in range(___number_of_pulses):
            EXTRACTED_CURRENT = []
            __RAW___ = RAW_SIGNAL[int((___delay + ___current_delay +(___anal_freq*j))/sampling):int((___delay + ___current_delay + ___window +(___anal_freq*j))/sampling)]
            event_latency__ .append(int(__RAW___.tolist().index(np.min(__RAW___[0:int(___window/sampling)])))*sampling)

            for k in range(int(___window/sampling)):
                EXTRACTED_CURRENT.append(np.mean(RAW_SIGNAL[int((___delay + ___current_delay + 0.5 +___anal_freq*j)/sampling)+k:int((___delay + ___current_delay + 0.5 +___anal_freq*j)/sampling)+k+5]))
            _EPSC_NO_OFFSET.append(EXTRACTED_CURRENT)

            temp__ = __RAW___[int(event_latency__[len(event_latency__)-1]/sampling)::]
            tOFFSET.append(temp__)
            occOFFSET.append(np.linspace(0, len(temp__), len(temp__)))
        for j in range(___number_of_pulses):
            popt_temp, pcov_temp = curve_fit(linear_func, occOFFSET[j], tOFFSET[j])
            popt.append(popt_temp)

        #This corrects values depending on fit
        for j in range(___number_of_pulses):
            EXTRACTED_CURRENT = []
            EXTRACTED_CURRENT = RAW_SIGNAL[int((___delay + ___current_delay  +___anal_freq*j)/sampling):int((___delay + ___window +___anal_freq*j)/sampling)]
            if j > 0:
                DECAY_FIT = []
                preciseOFFSET = 0
                __start__ = (___anal_freq)/sampling
                __stop__ = (___anal_freq +___window) /sampling

                for k in np.linspace(__start__, __stop__, 30):
                    if preciseOFFSET == 0:
                        preciseOFFSET = linear_func((___anal_freq- event_latency__[j-1])/sampling, popt[j-1][0], popt[j-1][1])

                    if linear_func(k+___anal_freq, popt[j-1][0], popt[j-1][1]) < 0:
                        DECAY_FIT.append(linear_func(k+___anal_freq, popt[j-1][0], popt[j-1][1]))

                    if len(DECAY_FIT) > 0:
                        OFFSET_CHARGE = (np.trapz(DECAY_FIT)*((___window/sampling)/len(DECAY_FIT))*sampling)

                    if preciseOFFSET > 0:
                        preciseOFFSET = 0
            else:
                OFFSET_CHARGE = 0
                preciseOFFSET = 0
                
            current_value = abs(np.min(EXTRACTED_CURRENT))-abs(3*sdNOISE)
            charge_value = np.trapz(EXTRACTED_CURRENT)*sampling
            
            __RAW___ = RAW_SIGNAL[int((___delay + ___current_delay +(___anal_freq*j))/sampling):int((___delay +___window+(___anal_freq*j))/sampling)]

            peak_occurence.append(int(__RAW___.tolist().index(np.min(__RAW___)))*sampling)
            current_value = abs(current_value)
            AMP_VALUES.append(current_value)
            CHARGE_VALUES.append(charge_value)
            CORRECTED_AMP_VALUES.append(current_value-(abs(preciseOFFSET)))
            CORRECTED_CHARGE_VALUES.append(abs(charge_value)-abs(OFFSET_CHARGE))
            FULL_EPSC_TRAIN_.append(EXTRACTED_CURRENT-preciseOFFSET)

        SYN_name.append('TEMP')

        FULL_PEAK_OCCURENCE_.append(peak_occurence)
        FULL_RAW_TRACES .append(_EPSC_NO_OFFSET)
        FULL_CHARGE_NOISE_LIST.append(chNOISE)
        FULL_AMP_NOISE_LIST.append(3*sdNOISE)

        fig = plt.figure(figsize=(8, 2), num=signal_name)
        ax = fig.add_subplot(132)
        ax.scatter(np.linspace(1, ___number_of_pulses, ___number_of_pulses), AMP_VALUES, c="black", alpha=0.1, linewidth=3, s=10)
        ax.scatter(np.linspace(1, ___number_of_pulses, ___number_of_pulses), CORRECTED_AMP_VALUES, c="red", alpha=0.5, s=15)
        ax.plot(np.linspace(1, ___number_of_pulses, ___number_of_pulses), CORRECTED_AMP_VALUES, c="red", alpha=0.5)
        ax.plot(np.linspace(1, ___number_of_pulses, ___number_of_pulses), CORRECTED_AMP_VALUES, c="red", alpha=0.75, linestyle='--')
        sns.despine(left=False, right=True, top=True, bottom=False)
        plt.xlabel("EPSC AMP (pA)")
        plt.ylim(0, np.max(CORRECTED_AMP_VALUES)*1.5)

        ax = fig.add_subplot(133)
        ax.scatter(np.linspace(1, ___number_of_pulses, ___number_of_pulses), CHARGE_VALUES, c="black", alpha=0.1, linewidth=3, s=10)
        ax.scatter(np.linspace(1, ___number_of_pulses, ___number_of_pulses), CORRECTED_CHARGE_VALUES, c="red", alpha=0.5, s=15)
        ax.plot(np.linspace(1, ___number_of_pulses, ___number_of_pulses), CORRECTED_CHARGE_VALUES, c="red", alpha=0.5)
        ax.plot(np.linspace(1, ___number_of_pulses, ___number_of_pulses), CORRECTED_CHARGE_VALUES, c="red", alpha=0.75, linestyle='--')
        sns.despine(left=False, right=True, top=True, bottom=False)
        plt.ylim(0, np.max(CORRECTED_CHARGE_VALUES)*1.5)
        plt.xlabel("EPSC CHARGE (fC)")


        ax = fig.add_subplot(131)
        ax.scatter(np.linspace(___delay, ___delay+___anal_freq*9, ___anal_freq), np.linspace(0, 0, ___anal_freq), color='black', s=100, alpha=0.9, marker='|')
        ax.plot(-np.max(AMP_VALUES)*1.5, np.max(AMP_VALUES)*1.5)
        ax.plot(np.linspace(___delay-3, ___delay+___anal_freq*9-3, 11), np.linspace(0, 0, 11), color="red", alpha=0.7, linewidth=2)
        for i in range(___number_of_pulses):
            if CORRECTED_AMP_VALUES[i] != 0:
                ax.scatter(___delay+___current_delay+peak_occurence[i]+i*___anal_freq, (0-CORRECTED_AMP_VALUES[i]), marker='o', color='orange', alpha=0.7)
                ax.plot(np.linspace((___delay+ ___current_delay+___anal_freq*i),(___delay+ ___current_delay+ ___anal_freq*i+len(FULL_EPSC_TRAIN_[i])*sampling), len(FULL_EPSC_TRAIN_[i])), FULL_EPSC_TRAIN_[i], c='black', alpha=0.8)


        ax.scatter(np.linspace(___delay, ___delay+___anal_freq*9, ___number_of_pulses), np.linspace(0, 0, ___number_of_pulses), color='black', s=100, alpha=0.9, marker='|')
        ax.set_xlim(___delay-20, ___delay+___anal_freq*9+20)
        ax.set_ylim(10,-np.max(CORRECTED_AMP_VALUES)*1.5)
        sns.despine(left=False, right=True, top=True, bottom=False)
        plt.xlabel("RAW TRACE")
        plt.tight_layout()
        pdf.savefig(fig)

    DATA_f = pd.DataFrame(data=CORRECTED_AMP_VALUES)
    DATA_f.to_csv(str(signal_name)+' Amplitude_temp_.csv', header=None, sep=' ')
    DATA_f2 = pd.DataFrame(data=CORRECTED_CHARGE_VALUES)
    DATA_f2.to_csv(str(signal_name)+' Charge_temp_.csv', header=None, sep=' ')
    
    print("__ANALYSIS SUCCEDED__")
    print("TIME : "+str(time.time()-__TIMESTART__)+" Sec.")

    return