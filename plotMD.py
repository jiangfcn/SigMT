# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 14:14:08 2021

@author: AJITHABH K. S.
Last modified: 21-07-2022

This script can be used to see the Mahalanobis distance (MD) for
all time windows/events for all target frequencies.

Any impedance
value can be selected for plotting. The color of the data points
indicate the coherency value.

Change following three lines to see different impedance values
and MD for Ex and Ey components
Z2_all = bandavg.get('Zxy_single') 'Change xy to yx'
Z2_mcd = bandavg.get('Zxy_mcd') 'Change xy to yx'
maha = bandavg.get('mahal_robustEx') 'Change Ex to Ey'

"""

### Colorbar 
import matplotlib as mpl
def reverse_colourmap(cmap, name = 'my_cmap_r'):
    reverse = []
    k = []   

    for key in cmap._segmentdata:    
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []

        for t in channel:                    
            data.append((1-t[0],t[2],t[1]))            
        reverse.append(sorted(data))    

    LinearL = dict(zip(k,reverse))
    my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL) 
    return my_cmap_r


### Ploting of Mahalanobis distance
import matplotlib
import numpy as np
#matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
cdict = {'red': ((0.0, 0.0, 0.0),
                 (0.1, 0.5, 0.5),
                 (0.2, 0.0, 0.0),
                 (0.4, 0.2, 0.2),
                 (0.6, 0.0, 0.0),
                 (0.8, 1.0, 1.0),
                 (1.0, 1.0, 1.0)),
        'green':((0.0, 0.0, 0.0),
                 (0.1, 0.0, 0.0),
                 (0.2, 0.0, 0.0),
                 (0.4, 1.0, 1.0),
                 (0.6, 1.0, 1.0),
                 (0.8, 1.0, 1.0),
                 (1.0, 0.0, 0.0)),
        'blue': ((0.0, 0.0, 0.0),
                 (0.1, 0.5, 0.5),
                 (0.2, 1.0, 1.0),
                 (0.4, 1.0, 1.0),
                 (0.6, 0.0, 0.0),
                 (0.8, 0.0, 0.0),
                 (1.0, 0.0, 0.0))}
my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
my_cmap = reverse_colourmap(my_cmap)

#coh_selected_all = bandavg.get('coh_selectedEx')
coh_selected_all = np.ones(np.shape(bandavg.get('ExExc')),dtype=float)
Z2_all = bandavg.get('Zxy_single')
Z2_mcd = bandavg.get('Zxy_mcd')
maha = bandavg.get('mahal_robustEx')
for fnum in range(np.size(ftlist)):
    coh_selected = coh_selected_all[fnum,:].reshape(-1,1)
    ind_coh = np.where(coh_selected==0)[0].reshape(-1,1)
    Z2 = Z2_all[fnum,:]
    Z2 = np.delete(Z2,ind_coh).reshape(-1,1)
    #cc = AllcohEx[fnum,:]
    cc = maha[fnum,:]
    cc = cc.reshape(-1,1)
    cc = np.delete(cc,ind_coh).reshape(-1,1)
    fig, ax = plt.subplots()
    sc = plt.scatter(Z2.real, Z2.imag,c=cc,cmap=my_cmap)
    plt.colorbar(sc)
    plt.clim(0,10) 
    plt.scatter(Z2_mcd.real[fnum],Z2_mcd.imag[fnum],c='black',marker=(5, 1))
    plt.title(procinfo.get('selectedsite') + ' - ' + procinfo.get('meas')+
             ' ('+str(procinfo.get('fs'))+' Hz) f='+ str(round(ftlist[fnum][0],2)) +' Hz')
    #plt.scatter(np.mean(Zxy_singleR),np.mean(Zxy_singleI))
    # plt.savefig('C:/Users/Ajithabh/Desktop/myImagePDF.eps', format='eps', dpi=1200)
    

