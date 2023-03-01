# -*- coding: utf-8 -*-
"""
Created on Mon May  4 16:49:35 2020
@author: AJITHABH K. S.
Last modified: 06-01-2023

##########
Program: SigMT
Program for the processing of MT time series.
Authors: Ajithabh K. S. and Prasanta K. Patro
CSIR - National Geophysical Research Institute, Hyderabad, India.
##########

This code controls all operations in this package.

Give the path of the folder where your sites are kept in the
'project_path' variable.

Then inputs such as site name and measurements will be asked in the console.

In case decimation is required, make 'dflag = 1' and give decimation sequence 
in 'decimate' variable.

Make 'ctflag =1' to enable coherency threshold.

Make 'pdflag =1' to enable polarization direction based data selection.

Read user manual for more details about data selection
 
"""
# Importing necessary modules
import config, mtproc, var, data_sel, tipper, mahaDist, plotting
from matplotlib import pyplot as plt
from scipy import signal
import numpy as np
import time
import math
#
config = config.configuration()
## Provide project path where sites are kept
project_path = 'A:/ProgramCode/SigMT/data/'
# #
# #========= Selection of site and setting a path =========
sites, selectedsite, measid, all_meas, select_meas, \
    proc_path = mtproc.makeprocpath(project_path)
#
#========= Site is selected and path is created =========
#========= Time series reading starts =========
procinfo = {}
ts, procinfo['fs'], procinfo['sensor_no'], timeline, \
    procinfo['ChoppStat'], loc = mtproc.ts(proc_path)
#========= Decimation section ================= 
# Keep dflag = 0 if decimation is not required
dflag = 0
if dflag == 1:
    decimate = [8,8,4]
    for d in decimate:
        ts['tsEx'] = signal.decimate(ts.get('tsEx'), d, n=None, ftype='iir')
        ts['tsEy'] = signal.decimate(ts.get('tsEy'), d, n=None, ftype='iir')
        ts['tsHx'] = signal.decimate(ts.get('tsHx'), d, n=None, ftype='iir')
        ts['tsHy'] = signal.decimate(ts.get('tsHy'), d, n=None, ftype='iir')
        ts['tsHz'] = signal.decimate(ts.get('tsHz'), d, n=None, ftype='iir')
        procinfo['fs'] = procinfo.get('fs')/d
#========= Decimation section end =============
#
# Some calculations and printing some information
# No need to edit
procinfo['nofs'] = len(ts['tsEx'])
procinfo['notch'] = 0 # Notch flag 1 - On, 0 - Off
print('--------------------')
print('MT site: ' + selectedsite)
print('Measurement directory: ' + all_meas[select_meas])
print('fs= '+ str(procinfo.get('fs')) +' Hz')
print('Sensor numbers:')
print(procinfo.get('sensor_no'))
print('Length of time series = ' + str(procinfo.get('nofs')))
print('--------------------')
procinfo['meas'] = all_meas[select_meas]
procinfo['proc_path'] = proc_path
procinfo['selectedsite'] = selectedsite
del all_meas, select_meas, selectedsite
del proc_path, project_path, sites
print('Unused variables deleted.')
print('\n\n--------------------')
# Find out window length
if config.get('FFT_Length') == 0:
    procinfo['WindowLength'] = mtproc.FFTLength(procinfo.get('nofs'))
else:
    procinfo['WindowLength'] = config.get('FFT_Length')
procinfo['overlap'] = 50 #Input in percentage %
print('\nWindow Length selected: '+ str(procinfo.get('WindowLength')))
procinfo['nstacks'] = math.floor(procinfo.get('nofs')/procinfo.get('WindowLength'))
procinfo['nstacks'] = (procinfo.get('nstacks') * 2) - 1
print('Time series overlap: ' + str(procinfo.get('overlap'))+'%')
print('No. of windows: '+ str(procinfo.get('nstacks')))
print('--------------------')
#
#
#==================== Start band averaging ====================
# No need to edit
# Band average value after calibration and averaging using parzen window
ftlist,bandavg = mtproc.bandavg(ts,procinfo,config)
#
#==================== Band averaging finished =================
#
#
#====Data selection tools section. Coherency threshold & Polarization direction
cohMatrixEx = np.ones(np.shape(bandavg.get('ExExc')),dtype=float)
cohMatrixEy = np.ones(np.shape(bandavg.get('ExExc')),dtype=float)
pdmat = np.ones(np.shape(bandavg.get('ExExc')),dtype=float)
# Calculation of coherency values for all time windows
AllcohEx = data_sel.cohEx(bandavg)
AllcohEy = data_sel.cohEy(bandavg)
# Calculation of polarization directions for all time windows
alpha_degH,alpha_degE = data_sel.pdvalues(bandavg)
#
#====== Coherency threshold ======
#JF,20230117
#precentage of small coherencies to be deleted: 
nprecent = 0.8
nleast = 16  #at least number of segments
sortcohEx = np.sort(AllcohEx,axis=1)
sortcohEy = np.sort(AllcohEy,axis=1)

ncolumn = np.shape(sortcohEx)[1]
nperiod = np.shape(sortcohEx)[0]
if (1-nprecent)*ncolumn < nleast:
    nprecent = 1-(nleast/ncolumn)
if ncolumn < nleast: 
    nprecent = 0

CohThre = np.ones(nperiod)
for icoh in range(nperiod):
    thre_Ex = sortcohEx[icoh,int(nprecent*ncolumn)]
    thre_Ey = sortcohEy[icoh,int(nprecent*ncolumn)]
    CohThre[icoh] = np.amin([thre_Ex,thre_Ey])
    if CohThre[icoh] > 0.9995: #JF
        CohThre[icoh] = 0.9995

ctflag = 0 # Give '1' to perform coherency threshold based selection
if ctflag == 1:
    #CohThre = [0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9]
    for i in range(np.shape(AllcohEx)[0]):
        for j in range(np.shape(AllcohEx)[1]):
            if AllcohEx[i,j] < CohThre[i]:
                cohMatrixEx[i,j] = 0
            else:
                cohMatrixEx[i,j] = 1
            if AllcohEy[i,j] < CohThre[i]:
                cohMatrixEy[i,j] = 0
            else:
                cohMatrixEy[i,j] = 1
#
#====== Polarization direction ======
# Degrees
Hpdflag = 10 # Give '0' to discard polarization angles
Epdflag = 10
if Hpdflag == 0: # H-based
    pdlim = [-72,-68] # discard
    alpha = alpha_degH # Use either alpha_degE or alpha_degH
    for i in range(np.shape(pdmat)[0]):
        for j in range(np.shape(pdmat)[1]):
            if alpha[i,j] > pdlim[0] and alpha[i,j] < pdlim[1]:
                pdmat[i,j] = 0
            else:
                pdmat[i,j] = 1
if Epdflag == 0: # E-based
    pdlim = [-65,-55] # discard
    alpha = alpha_degE # Use either alpha_degE or alpha_degH
    for i in range(np.shape(pdmat)[0]):
        for j in range(np.shape(pdmat)[1]):
            if alpha[i,j] > pdlim[0] and alpha[i,j] < pdlim[1]:
                pdmat[i,j] = 0
            else:
                pdmat[i,j] = 1
# time-windows            
pdflag = 10
timewindow_limits = [150,350] # time-window discard 
# Give '1' -Save time-windows
# Give '0' -Discard time-windows 
if pdflag == 0:    
    for i in range(np.shape(pdmat)[0]):
        for j in range(np.shape(pdmat)[1]):
            if j > timewindow_limits[0] and j < timewindow_limits[1]:
                pdmat[i,j] = 0
            else:
                pdmat[i,j] = 1
#JF
if pdflag == 1:    
    for i in range(np.shape(pdmat)[0]):
        for j in range(np.shape(pdmat)[1]):
            if j > timewindow_limits[0] and j < timewindow_limits[1]:
                pdmat[i,j] = 1
            else:
                pdmat[i,j] = 0
#====== END Polarization direction ==========================================

#
bandavg['cohMatrixEx'] = cohMatrixEx
bandavg['cohMatrixEy'] = cohMatrixEy
bandavg['pre_sel_matEx'] = cohMatrixEx * pdmat
bandavg['pre_sel_matEy'] = cohMatrixEy * pdmat
bandavg['mdmatrixEx'],bandavg['Zxx_mcd'],bandavg['Zxy_mcd'],bandavg['mahal_robustEx'] = mahaDist.mcd(bandavg,'Ex',config)
bandavg['mdmatrixEy'],bandavg['Zyx_mcd'],bandavg['Zyy_mcd'],bandavg['mahal_robustEy'] = mahaDist.mcd(bandavg,'Ey',config)
bandavg['selectedEx'] = bandavg.get('mdmatrixEx')
bandavg['selectedEy'] = bandavg.get('mdmatrixEy')
bandavg['avgt'] = np.sum((bandavg.get('selectedEx'))!=0,axis=1)
del cohMatrixEx,cohMatrixEy
#
#
#=========== Tipper estimation ===================================
[TxAll, TyAll] = tipper.tippall(bandavg)
mahaWtTx, Tx_mcd_mean = tipper.mcd(TxAll,config)
mahaWtTy, Ty_mcd_mean = tipper.mcd(TyAll,config)
bandavg['tipp_selected'] = mahaWtTx * mahaWtTy
[Tx, Ty] = tipper.tipper(bandavg)
[TxVar, TyVar] = tipper.tipperVar(bandavg)
#=========== Tipper estimation DONE===============================
#
#==================== Robust estimation begins ===================
#
Z_huber = mtproc.perform_robust(ftlist,bandavg)
#
#==================== Calculation of variance ===================
Zvar = {}
Zvar['xx'],Zvar['xy'],cohEx = var.ZExvar(Z_huber,bandavg)
Zvar['yx'],Zvar['yy'],cohEy = var.ZEyvar(Z_huber,bandavg)
#
#
### Plotting figures ###
plotting.plotfigs(procinfo, ftlist, Z_huber, Zvar, Tx, Ty, cohEx, cohEy)
print('Finished.')