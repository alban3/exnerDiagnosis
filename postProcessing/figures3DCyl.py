#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filename: figures3D.py
Date: 2025-11-13
Version: 1.0
License: GPL-3.0-or-later
Copyright 2025 Alban Gilletta, Cyrille Bonamy, Julien Chauchat
This file is part of ExnerAnalysis.

ExnerAnalysis is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

ExnerAnalysis is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with ExnerAnalysis. If not, see <https://www.gnu.org/licenses/>.

Contact: alban.gilletta-de-saint-joseph@univ-grenoble-alpes.fr
"""

#
# Import section
#

import subprocess
import sys
import numpy as np
import os, glob
import fluidfoam
from pylab import figure, subplot, axis, xlabel, ylabel, show, savefig, plot
from pylab import title, matplotlib
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, LogNorm
import matplotlib.cm
from matplotlib.lines import Line2D
from netCDF4 import Dataset
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from massBalance3DCyl import Threedimsimu
from func1 import create_point_2Dcyl, save_point, read_point
from func1 import uns2cylvec, uns2cyl
import h5py

def mobileAvg(Tp, array, time):
    import numpy as np
    n = 0
    dt = time[1] - time[0]
    Nt = len(array)
    arrayAvg = np.zeros(Nt)
    while n < Nt-1:
        if time[n] < Tp:
            arrayAvg[n] = np.sum(array[0:n+1],axis=0)/(n+1)
            n += 1
        else:
            arrayAvg[n] = np.sum(array[int(n+1-Tp/dt):n+1],axis=0)/(Tp/dt)
            n += 1
            if n == Nt-1 :
                break
    arrayAvg[Nt-1] = np.sum(array[int(n+1-Tp/dt):],axis=0)/(Tp/dt)
    return arrayAvg

Tp = 2

caseM0 = "M0"
caseM1 = "M1"
caseM2 = "M2"
caseM3 = "M3"
caseM4 = "M4"
caseM5 = "M5"

caseList = [caseM5]

read_point_nc = True        # True if the file pointerpostproc.nc has been created
save_point_nc = False       # False if the file pointerpostproc.nc has been saved
meshCreated = True          # True if the file read_mesh.h5 has been created
resultCreated = True        # True if the file read_massBalance'+str(self.asint).replace('.','p')+'.h5 has been created

turbModel = 'LESdyn'
mesh = "M5"
case = turbModel+mesh
basepath = "../Configuration/3DCyl/"+ turbModel + "/"

### USEFUL IF RESULTS HAVE BEEN DIVIDED IN TIME -> JOIN RESULTS ###

mergeResult = False
if mergeResult:
    asintList = [0.60]
    for asint in asintList:
        with h5py.File(basepath+case+'/constant'+'/read_massBalanceZ1'+str(asint).replace('.','p')+'.h5',mode='w') as h5fw:
            ldset1 = 26
            for i in range(ldset1):
                k = -1
                for h5name in sorted(glob.glob(os.path.join(basepath, case, 'constant', 'read_massBalanceZ1'+str(asint).replace('.','p')+'_t*.h5'))):
                    k = k + 1
                    h5fr = h5py.File(h5name, 'r')
                    if k == 0:
                        dset1 = list(h5fr.keys())
                        arr_data = h5fr[dset1[i]][:-1]
                        dslen = arr_data.shape[0]
                    if np.logical_and(k > 0, k < 3): #if more than 2 . nc uncomment
                        arr_data = np.concatenate((arr_data, h5fr[dset1[i]][:-1]))
                        dslen = arr_data.shape[0]
                    if k == 3: # if 2 .nc -> k != 0: // if more than 2 .nc -> k == 3:
                        arr_data = np.concatenate((arr_data, h5fr[dset1[i]][:]))
                        dslen = arr_data.shape[0]
                h5fw.create_dataset(dset1[i]+'_final',data=arr_data[:])  

simuList = []
k=-1
for case in caseList:
    k = k + 1
    if case == caseM0 :
        gC = 2.7416556078516723
        gZw = 41.432951063305104
        gZ1 = 0.5830419509076277
        gZ1b = 0.265286536915622
        nzw = 11
        nzbw = 42
        nz1 = 176
        nz1b = 52
    elif case == caseM1 :
        gC = 2.4448365874972424
        gZw = 62.33183979817063
        gZ1 = 0.21697282450595096
        gZ1b = 0.5419339335730116
        nzw = 11
        nzbw = 46
        nz1 = 148
        nz1b = 44
    elif case == caseM2 :
        gC = 2.5854939547329243
        gZw = 62.33183979817063
        gZ1 = 0.21697282450595096
        gZ1b = 0.5419339335730116
        nzw = 11
        nzbw = 46
        nz1 = 148
        nz1b = 44
    elif case == caseM3 :
        gC = 2.3292823634703046
        gZw = 112.1832662895043
        gZ1 = 0.0775727693174146
        gZ1b = 0.9351357050411135
        nzw = 12
        nzbw = 58
        nz1 = 150
        nz1b = 42
    elif case == caseM4 :
        gC = 2.43820920811315
        gZw = 112.1832662895043
        gZ1 = 0.0775727693174146
        gZ1b = 0.9351357050411135
        nzw = 12
        nzbw = 58
        nz1 = 150
        nz1b = 42
    elif case == caseM5 :
        gC = 2.522311151573334
        gZw = 112.1832662895043
        gZ1 = 0.0775727693174146
        gZ1b = 0.9351357050411135
        nzw = 12
        nzbw = 58
        nz1 = 150
        nz1b = 42
    simuList.append(Threedimsimu(basepath, case, gC, gZw, gZ1, gZ1b, nzw, nzbw, nz1, nz1b, read_point_nc, save_point_nc, meshCreated, resultCreated, 0.60))

simu5 = simuList[0]

matplotlib.rc('text', usetex = True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{amssymb} \usepackage{physics} \usepackage[dvipsnames]{xcolor}')

# plot creation and parameters
matplotlib.rcParams.update({'font.size': 12})
mpl.rcParams["lines.linewidth"] = 1.5
mpl.rcParams["lines.markersize"] = 10

plt.figure(figsize=[12, 15])
gs = plt.GridSpec(1, 1)
gs.update(left=0.08, right=0.975, top=0.95, bottom=0.15, wspace=0.1,
          hspace=0.1)
font_legend = matplotlib.font_manager.FontProperties(family='Comic Sans MS',
                                                     weight='bold',
                                                     style='normal', size=14)

fignum1 = plt.figure(num=1,figsize=(12,15),dpi=300)
ax1 = fignum1.add_subplot(111)
plt.xlim(-3e-5,3e-5) 
plt.ylim(-3e-5,3e-5)
plt.xlabel(r'\textcircled{1} (m$^{3}$ s$^{-1}$)',fontsize = '20')
plt.ylabel(r'-(\textcircled{2} + \textcircled{3} + \textcircled{4} + \textcircled{5}) (m$^{3}$ s$^{-1}$)', fontsize = '20') 
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.grid()

fignum2 = plt.figure(num=2,figsize=(12,15),dpi=300)
ax2 = fignum2.add_subplot(111)
plt.xlim(0,120)
plt.ylim(-1.8e-5,1.8e-5)
plt.xlabel(r'Time (s)',fontsize = '20')
plt.ylabel(r'Volume evolution (m$^{3}$ s$^{-1}$)', fontsize = '20')
plt.grid()
plt.tick_params(labelsize=20)

fignum3 = plt.figure(num=3,figsize=(12,15),dpi=300)
ax3 = fignum3.add_subplot(111)
plt.xlim(0,120)
plt.ylim(-1.8e-5,1.8e-5)
plt.xlabel(r'Time (s)',fontsize = '20')
plt.ylabel(r'Residual (m$^{3}$ s$^{-1}$)', fontsize = '20')
plt.grid()
ax3_0 = ax3.inset_axes([0.25,0.05,0.4,0.25])
ax3_0.axis([0, 120, -300e0, 300e0])
ax3_0.grid(linestyle=':')
ax3_0.set_ylabel(r'Residual ($\%$)', fontsize = '20')

fignum4 = plt.figure(num=4,figsize=(12,15),dpi=300)
ax4 = fignum4.add_subplot(111)
plt.xlim(0,120)
plt.ylim(-1.8e-5,1.8e-5)
plt.xlabel(r'Time (s)',fontsize = '20')
plt.ylabel(r'Volume evolution (m$^{3}$ s$^{-1}$)', fontsize = '20')
plt.grid()

#
# FIGURE 1
#

color = simu5.time[0:-1]
sc1 = ax1.scatter(simu5.phi_dzbed_sur_dt, -(simu5.intdVolPhi_sur_dt + simu5.intSlat_PindS + simu5.intS_PindS + simu5.intS2_PindS + simu5.intS3_PindS + simu5.intStop_PindS), marker ='o', s=50, c=color, cmap='seismic', edgecolors='none', vmin= simu5.time[0], vmax = simu5.time[-1])
ax1.plot(np.linspace(-10,10,100),np.linspace(-10,10,100),'--k',label ='y=x')
ax1.plot(np.linspace(-20,20,10),np.linspace(0,0,10),'-k', lw = 1)
ax1.plot(np.linspace(0,0,10),np.linspace(-20,20,10),'-k', lw = 1)
cbar1 = plt.colorbar(sc1,ax=ax1)
cbar1.set_label(r'Time (s)',fontsize= '20')
ax1.legend(prop={"size": 20})

#
# FIGURE 2
#

ax2.plot(simu5.time[:-1], simu5.phi_dzbed_sur_dt , alpha=0.4, color = 'orange', ls = '-', lw =1.25)
ax2.plot(simu5.time[:-1], simu5.intdVolPhi_sur_dt , alpha=0.4, color = 'b', ls = '-', lw =1.25)
ax2.plot(simu5.time[:-1], simu5.intSlat_PindS, alpha=0.4, color = 'g', ls = '-', lw =1.25)
ax2.plot(simu5.time[:-1], simu5.intS_PindS + simu5.intS2_PindS + simu5.intS3_PindS, alpha=0.4, color = 'm', ls = '-', lw =1.25)
ax2.plot(simu5.time[:-1], simu5.intStop_PindS, alpha=0.4, color = 'c', ls = '-', lw =1.25)
ax2.plot(simu5.time[:-1], mobileAvg(Tp, simu5.phi_dzbed_sur_dt, simu5.time) , color = 'orange', ls = '-', label = r'\textcircled{1} Bed evolution')
ax2.plot(simu5.time[:-1], mobileAvg(Tp, simu5.intdVolPhi_sur_dt, simu5.time) , color = 'b', ls = '-', label = r'\textcircled{2} Storage evolution')
ax2.plot(simu5.time[:-1], mobileAvg(Tp, simu5.intSlat_PindS, simu5.time) , color = 'g', ls = '-', label = r' \textcircled{3} Lateral flux')
ax2.plot(simu5.time[:-1], mobileAvg(Tp, simu5.intS_PindS + simu5.intS2_PindS + simu5.intS3_PindS, simu5.time), color = 'm', ls = '-', label = r' \textcircled{4} Bed flux')
ax2.plot(simu5.time[:-1], mobileAvg(Tp, simu5.intStop_PindS, simu5.time), color = 'c', ls = '-', label = r' \textcircled{5} Top flux')
ax2.plot(np.linspace(0,120,10),np.linspace(0,0,10),'-k', lw = 1)
ax2.legend(prop={"size": 24}, ncol = 2)

#
# FIGURE 3
#

ax3.plot(simu5.time[:-1], simu5.phi_dzbed_sur_dt + simu5.intdVolPhi_sur_dt + simu5.intSlat_PindS + simu5.intS_PindS + simu5.intS2_PindS + simu5.intS3_PindS, color = 'c', label = r'$\phi_b=0.60$' )
ax3.legend(prop={"size": 20})
ax3_0.plot(simu5.time[:-1], (simu5.phi_dzbed_sur_dt + simu5.intdVolPhi_sur_dt + simu5.intSlat_PindS + simu5.intS_PindS + simu5.intS2_PindS + simu5.intS3_PindS)/simu5.phi_dzbed_sur_dt*100, color = 'c' )

#
# FIGURE 4
#

ax4.plot(simu5.time[:-1], simu5.intSlat1stbed_PindS, ls = '-', color = 'b', label = r' \textcircled{3abed}')
ax4.plot(simu5.time[:-1], simu5.intSlat1stsusp_PindS, ls = ':', color = 'b', label = r' \textcircled{3asusp}')
ax4.plot(simu5.time[:-1], simu5.intSlat1endbed_PindS, ls = '-', color = 'r', label = r' \textcircled{3bbed}')
ax4.plot(simu5.time[:-1], simu5.intSlat1endsusp_PindS, ls = ':', color = 'r', label = r' \textcircled{3bsusp}')
ax4.plot(simu5.time[:-1], simu5.intSoutbed_PindS, ls = '-', color = 'g', label = r' \textcircled{3cbed}')
ax4.plot(simu5.time[:-1], simu5.intSoutsusp_PindS, ls = ':', color = 'g', label = r' \textcircled{3csusp}')
ax4.legend(prop={"size": 20})

for i in plt.get_fignums():
    if i == 1:
        plt.figure(i).set_size_inches(10,8)
        plt.figure(i).savefig('Figures/3DCyl/'+turbModel+'/'+mesh+'/'+'massBalance_'+case+'60'+'.eps', 
                facecolor='w', edgecolor='w', dpi=300, format='eps')
    if i == 2:
        ax2.set_rasterized(True)
        ax2.set_rasterization_zorder(0)
        plt.figure(i).set_size_inches(12,10)
        plt.figure(i).savefig('Figures/3DCyl/'+turbModel+'/'+mesh+'/'+'contrib_'+case+'60'+'.eps', 
                facecolor='w', edgecolor='w', dpi=300, format='eps')
    if i == 3:
        plt.figure(i).set_size_inches(10,8)
        plt.figure(i).savefig('Figures/3DCyl/'+turbModel+'/'+mesh+'/'+'residual_'+case+'60'+'.eps', 
                facecolor='w', edgecolor='w', dpi=300, format='eps')
    if i == 4:
        plt.figure(i).set_size_inches(10,8)
        plt.figure(i).savefig('Figures/3DCyl/'+turbModel+'/'+mesh+'/'+'LatFlux_'+case+'60'+'.eps', 
                facecolor='w', edgecolor='w', dpi=300, format='eps')     
