# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 12:27:42 2021
@author: AJITHABH K. S.
Last modified: 06-01-2023

##########
Program: SigMT
Program for the processing of MT time series.
Authors: Ajithabh K.S. and Prasanta K. Patro
CSIR - National Geophysical Research Institute, Hyderabad, India.
##########

This code controls remote reference method in MT data processing.

Give the path of the folder where your sites are kept in the
'project_path' variable.

Then inputs such as site name and measurements will be asked in the console.
First, the details of the local site will be asked. Then, the details of the 
remote site need to be entered.

In case decimation is required, make 'dflag = 1' and give decimation sequence 
in 'decimate' variable.

Make 'ctflag =1' to enable coherency threshold.

Make 'pdflag =1' to enable polarization direction based data selection.

Read user manual for more details about the data selection
"""
# Importing necessary modules
import mtprocRR, var, data_sel, tipperR, mahaDist, mtproc, config, plotting
from matplotlib import pyplot as plt
from scipy import signal
import numpy as np
import os
import time
import math
config = config.configuration()
# define project path where sites are kept
project_path = 'A:/ProgramCode/SigMT/data/'
# #
# #========= Selection of site and setting a path =========
sites, selectedsite, measid, all_meas, select_meas, proc_path = mtproc.makeprocpath(project_path)
# #========= Site is selected and path is created =========
# #========= Selection of remote site and setting a path ==
sites, selectedsiteR, measidR, all_measR, select_measR, proc_pathR = mtproc.makeprocpath(project_path)
# #========= Remote site is selected and path is created ==
#-------- Time series reading starts ------------------
#
#
procinfo = {}
procinfoR = {}
ts, procinfo['fs'], procinfo['sensor_no'], timeline, procinfo['ChoppStat'], loc = mtprocRR.ts(proc_path)
tsR, procinfoR['fs'], procinfoR['sensor_no'], timelineR, procinfoR['ChoppStat'], locR = mtprocRR.ts(proc_pathR)
#==============
#==============
# Prepare ts data for RR 
timeline = np.asarray(timeline)
timelineR = np.asarray(timelineR)
del tsR['tsEx']
del tsR['tsEy']
del tsR['tsHz']
tsRx = tsR.get('tsHx')
tsRy = tsR.get('tsHy')
#Rstart = timeline[0]
#Rend = timeline[-1]
#Rstart_ind = np.where(np.equal(timelineR,Rstart))[0][0]
#Rend_ind = np.where(np.equal(timelineR,Rend))[0][0]
#JF:
tsEx = ts.get('tsEx')
tsEy = ts.get('tsEy')
tsHx = ts.get('tsHx')
tsHy = ts.get('tsHy')
tsHz = ts.get('tsHz')
# Rstart = timeline[0]
# Rend = timeline[-1]
# Rstart_ind = np.where(np.equal(timelineR,Rstart))[0][0]
# Rend_ind = np.where(np.equal(timelineR,Rend))[0][0]
# tsRx = tsRx[Rstart_ind:Rend_ind+1]
# tsRy = tsRy[Rstart_ind:Rend_ind+1]

#JF:2023-02-07
Lstart = timeline[0]
Lend = timeline[-1]
Rstart = timelineR[0]
Rend = timelineR[-1]
if Rend <= Lstart:
    raise Exception('no overlap between Reference and local sites!')
elif Rend > Lstart and Rend <= Lend:
    Rclipend = Rend
    Rclipstart = max([Lstart,Rstart])
elif Rend > Lend:
    Rclipend = Lend
    Rclipstart = max([Lstart,Rstart])
Rstart_ind = np.where(np.equal(timelineR,Rclipstart))[0][0]
Rend_ind = np.where(np.equal(timelineR,Rclipend))[0][0]
Lstart_ind = np.where(np.equal(timeline,Rclipstart))[0][0]
Lend_ind = np.where(np.equal(timeline,Rclipend))[0][0]

tsRx = tsRx[Rstart_ind:Rend_ind+1]
tsRy = tsRy[Rstart_ind:Rend_ind+1]

tsEx = tsEx[Lstart_ind:Lend_ind+1]
tsEy = tsEy[Lstart_ind:Lend_ind+1]
tsHx = tsHx[Lstart_ind:Lend_ind+1]
tsHy = tsHy[Lstart_ind:Lend_ind+1]
tsHz = tsHz[Lstart_ind:Lend_ind+1]
#end of JF:2023-02-07


del tsR['tsHx']
del tsR['tsHy']
tsR['tsRx'] = tsRx
tsR['tsRy'] = tsRy

#JF, clipping ts
ts['tsEx'] = tsEx
ts['tsEy'] = tsEy
ts['tsHx'] = tsHx
ts['tsHy'] = tsHy
ts['tsHz'] = tsHz
# clipstart = time.strftime('%Y/%m/%d %H:%M:%S',Rclipstart[0]*24*3600)
# clipend = time.strftime('%Y/%m/%d %H:%M:%S',Rclipend[0]*24*3600)
# print('RR time intervel:',+ str(clipstart),'-',+str(clipend))
#end of JF,clipping ts
del timelineR, Rstart, Rend, Rstart_ind, Rend_ind, Lstart, Lend, Lstart_ind, Lend_ind
del tsRx, tsRy, tsEx, tsEy, tsHx, tsHy, tsHz
#del timelineR, Rstart, Rend
#del Rstart_ind, Rend_ind, tsRx, tsRy



#========= Decimation section ================= 
# Keep dflag = 0 if decimation is not required
#RR site
print('RR: input sampling freq = ' + str(procinfoR.get('fs')) +' Hz') 
Rdflag = 1
if Rdflag == 1:
    Rdecimate = [8,8,4]
    for Rd in Rdecimate:        
        tsR['tsRx'] = signal.decimate(tsR.get('tsRx'), Rd, n=None, ftype='iir')
        tsR['tsRy'] = signal.decimate(tsR.get('tsRy'), Rd, n=None, ftype='iir')
        procinfoR['fs'] = procinfoR.get('fs')/Rd
print('RR: output sampling freq = ' + str(procinfoR.get('fs')) +' Hz') 
#local site
print('Local: input sampling freq = ' + str(procinfo.get('fs')) +' Hz')

dflag = 1
if dflag == 1:
    decimate = [8,8,4]
    for d in decimate:
        ts['tsEx'] = signal.decimate(ts.get('tsEx'), d, n=None, ftype='iir')
        ts['tsEy'] = signal.decimate(ts.get('tsEy'), d, n=None, ftype='iir')
        ts['tsHx'] = signal.decimate(ts.get('tsHx'), d, n=None, ftype='iir')
        ts['tsHy'] = signal.decimate(ts.get('tsHy'), d, n=None, ftype='iir')
        ts['tsHz'] = signal.decimate(ts.get('tsHz'), d, n=None, ftype='iir')
        procinfo['fs'] = procinfo.get('fs')/d
print('Local: output sampling freq = ' + str(procinfo.get('fs')) +' Hz') 
print('decimation have done!') 
#end of JF's edit
#========= Decimation section end =============
#
# Some calculations and printing some information
# No need to edit
procinfo['nofs'] = len(ts['tsEx']) #number of samples
procinfo['notch'] = 0 # Notch flag 1 - On, 0 - Off
print('--------------------')
#print('\nSaved time series after trend & bias removal.')
print('MT site: ' + selectedsite)
print('Measurement directory: ' + all_meas[select_meas])
print('fs= '+ str(procinfo.get('fs')) +' Hz')
print('Sensor numbers:')
print(procinfo.get('sensor_no'))
print('Length of time series = ' + str(procinfo.get('nofs')))
print('--------------------')
procinfo['meas'] = all_meas[select_meas]
procinfo['proc_path'] = proc_path
procinfoR['proc_path'] = proc_pathR
procinfo['selectedsite'] = selectedsite
del all_meas, select_meas, selectedsite
del all_measR, select_measR, selectedsiteR, proc_pathR
del proc_path, project_path, sites
print('Unused variables deleted.')
print('\n\n--------------------')
# Find out window length
if config.get('FFT_Length') == 0:
    procinfo['WindowLength'] = mtproc.FFTLength(procinfo.get('nofs'))
else:
    procinfo['WindowLength'] = config.get('FFT_Length')
procinfo['overlap'] = 50 #Input in percentage 
print('\nWindow Length selected: '+ str(procinfo.get('WindowLength')))
procinfo['nstacks'] = math.floor(procinfo.get('nofs')/procinfo.get('WindowLength'))
procinfo['nstacks'] = (procinfo.get('nstacks') * 2) - 1
print('Time series overlap: ' + str(procinfo.get('overlap'))+'%')
print('No. of windows: '+ str(procinfo.get('nstacks')))
print('--------------------')
#==================== Start band averaging ====================
# No need to edit
# Band average value after calibration and averaging using parzen window
ftlist,bandavg = mtprocRR.bandavg(ts,procinfo,tsR,procinfoR,config)
#
#==================== Band averaging finished ====================
#########
spmat = mtprocRR.cleanSpec(bandavg)
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
nprecent = 0.7
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
        
        
ctflag = 1
if ctflag == 1:
    #CohThre = [0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9]
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
pdflag = 0
if pdflag == 1:
    pdlim = [-65,-55]
    alpha = alpha_degE #Use either alpha_degE or alpha_degH
    for i in range(np.shape(pdmat)[0]):
        for j in range(np.shape(pdmat)[1]):
            if alpha[i,j] > pdlim[0] and alpha[i,j] < pdlim[1]:
                pdmat[i,j] = 0
            else:
                pdmat[i,j] = 1
#            
pdflag = 0
if pdflag == 1:
    timewindow_limits = [140,350]
    for i in range(np.shape(pdmat)[0]):
        for j in range(np.shape(pdmat)[1]):
            if j > timewindow_limits[0] and j < timewindow_limits[1]:
                pdmat[i,j] = 0
            else:
                pdmat[i,j] = 1
#
bandavg['cohMatrixEx'] = cohMatrixEx
bandavg['cohMatrixEy'] = cohMatrixEy
bandavg['pre_sel_matEx'] = cohMatrixEx * pdmat
bandavg['pre_sel_matEy'] = cohMatrixEy * pdmat
bandavg['mdmatrixEx'],bandavg['Zxx_mcd'],bandavg['Zxy_mcd'],bandavg['mahal_robustEx'] = mahaDist.mcd(bandavg,'Ex',config)
bandavg['mdmatrixEy'],bandavg['Zyx_mcd'],bandavg['Zyy_mcd'],bandavg['mahal_robustEy'] = mahaDist.mcd(bandavg,'Ey',config)
bandavg['selectedEx'] = bandavg.get('mdmatrixEx') * spmat
bandavg['selectedEy'] = bandavg.get('mdmatrixEy') * spmat
#bandavg['tipp_selected'] = cohMatrixEx * cohMatrixEy * selmatTx * selmatTy

# bandavg['coh_selected'] = selmatEx * selmatEy
bandavg['avgt'] = np.sum((bandavg.get('selectedEx'))!=0,axis=1)
del cohMatrixEx,cohMatrixEy
#
#=========== Tipper estimation ===================================
[TxAll, TyAll] = tipperR.tippall(bandavg)
mahaWtTx, Tx_mcd_mean = tipperR.mcd(TxAll,config)
mahaWtTy, Ty_mcd_mean = tipperR.mcd(TyAll,config)
bandavg['tipp_selected'] = mahaWtTx * mahaWtTy
[Tx, Ty] = tipperR.tipper(bandavg)
[TxVar, TyVar] = tipperR.tipperVar(bandavg)
#=========== Tipper estimation DONE===============================
#==================== Robust estimation begins ===================
#
Z_huber = mtprocRR.perform_robust(ftlist,bandavg)
#
#==================== Calculation of variance ===================
Zvar = {}
Zvar['xx'],Zvar['xy'],cohEx = var.ZExvar(Z_huber,bandavg)
Zvar['yx'],Zvar['yy'],cohEy = var.ZEyvar(Z_huber,bandavg)
# ----------------------------------------------------------------------
#
#
### Plotting figures ###
plotting.plotfigs(procinfo, ftlist, Z_huber, Zvar, Tx, Ty, cohEx, cohEy)
print('Finished.')