import os
from os.path import exists
import subprocess
import time
from model import op2_reading
import numpy as np
from pyNastran.op2.op2 import read_op2, OP2
from scipy.signal import resample, find_peaks, savgol_filter
import scipy
from scipy.interpolate import interp1d
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt

op2_file = read_op2("Y:/Documents/NX Simulations/Nastran_test/r1-50_msh-6_tw-0.5_ts-2.0_tt-25.0_x-400_y-600_nsm-None/model_pynastran.op2", debug=False)
s_v = list()
for i in range(op2_file.get_result("rms.cquad8_stress")[2].ntotal):
    s11 = op2_file.get_result("rms.cquad8_stress")[2].data[0, i, 0]/1000
    s22 = op2_file.get_result("rms.cquad8_stress")[2].data[0, i, 1]/1000
    s12 = op2_file.get_result("rms.cquad8_stress")[2].data[0, i, 2]/1000
    s_v.append(np.sqrt(s11**2 - s11*s22 + s22**2 + 3*s12**2))
#count = 0
#for ele_11 in s_v:
#    if np.round(ele_11, 4) == 0.2864:
#        print(count)
#        print(ele_11)
#    count += 1


np.max(s_v)

#plt.plot(s_v)
#plt.show()
### Ancien bloc ###

#freqs = np.array(op2_file.get_result("modal_contribution.cquad8_stress")[1].freqs)
#max_von_skin = np.zeros(len(freqs))
#max_von_core = np.zeros(len(freqs))
#omax_core = np.zeros(len(freqs))
#for i, f in enumerate(freqs):
#    max_von_core[i] = np.max(op2_file.get_result("modal_contribution.cquad8_stress")[1].data[i, :, 7]) / 1000
#    max_von_skin[i] = np.max(op2_file.get_result("modal_contribution.ctria6_stress")[1].data[i, :, 7]) / 1000
#
## freq_max = freqs[np.max(max_von_skin)==max_von_skin]
## spl_size = 1000
##max_von_int = interp1d(freqs, max_von_skin)
##max_von_resample = max_von_int(np.linspace(20, 2000, spl_size))
#eig_f = np.round(op2_file.eigenvalues[''].radians / 2 / np.pi, 4)
#eigen_i = list()
#for f in freqs:
#    eigen_i.append(np.any(np.round(f, 4)==eig_f))
#Mean2_Von = np.sqrt(np.sum(max_von_skin[eigen_i]))  #
#Mean_Von = trapezoid(max_von_skin, freqs)/(2000-20)
#
#TSS_Von = np.sqrt(np.sum(max_von_skin**2))  #
#
#print("Sum sqrt eig: " + str(Mean2_Von))
#print("Moyenne : " + str(Mean_Von))
#print("RSS Von-Mises? : " + str(TSS_Von))
#
#plt.figure()
#plt.plot(freqs, max_von_core)
#plt.plot(freqs, freqs**0*Mean2_Von)
#plt.plot(freqs, freqs**0*Mean_Von)
#plt.legend(['Von Mises', "sqrt(Int)", "Mean"])
#plt.show()


### Plus ancien ###

#peaks1, _ = find_peaks(max_von_skin, distance=1)
#peaks_filt = list()
#for p in peaks1:
#    # VOIR SI ÇA FUCK AVEC MODE à 20HZ ou 2000Hz
#    if max_von_skin[p] > max_von_skin[p-2] and max_von_skin[p] > max_von_skin[p+2]:
#        peaks_filt.append(p)
#print(freqs[peaks_filt])
#print(max_von_skin[peaks_filt])


#filtered = savgol_filter(max_von_skin, 51, 3)





#op2_filename = "Y:/Documents/NX Simulations/Nastran_Mesh/r1-50_msh-4_tw-0.5_ts-2.0_tt-25.0_x-400_y-600_nsm-None/model_pynastran.op2"
#model_op2 = read_op2(op2_filename)
#type(model_op2.get_result("modal_contribution.cquad8_stress")[1].data)
#print(model_op2.eigenvalues[''].radians/2/np.pi)




