# -*- coding: utf-8 -*-
"""
Created on Sun Jul 24 13:36:21 2022

@author: AJITHABH K. S.
Last modified: 21-07-2022

This module plots four figures after the processing is completed.

1. Apparent resistivity and phase.
2. Coherency values for Ex and Ey components.
3. Tipper amplitudes and phase
4. Tipper real and imaginary components.
"""
  
def plotfigs(procinfo, ftlist, Z_huber, Zvar, Tx, Ty, cohEx, cohEy):
    import numpy as np
    from matplotlib import pyplot as plt
    # ----------------------------------------------------------------------
    # Apparant resistivities and phase
    rho_xy = (0.2/ftlist) * ((abs(Z_huber.get('Zxy'))) ** 2)
    rho_yx = (0.2/ftlist) * ((abs(Z_huber.get('Zyx'))) ** 2)
    phase_xy = np.degrees(np.arctan(Z_huber.get('Zxy').imag/Z_huber.get('Zxy').real))
    phase_yx = np.degrees(np.arctan(Z_huber.get('Zyx').imag/Z_huber.get('Zyx').real))
    #
    #
    # Errors for app. resistivity and phase
    err_rxy = (0.4/ftlist) * ((abs(Z_huber.get('Zxy')))) * Zvar.get('xy')
    err_ryx = (0.4/ftlist) * ((abs(Z_huber.get('Zyx')))) * Zvar.get('yx')
    err_pxy = np.degrees((Zvar.get('xy')) / abs(Z_huber.get('Zxy')))
    err_pyx = np.degrees((Zvar.get('yx')) / abs(Z_huber.get('Zyx')))
    err_rxy = err_rxy.reshape((-1,))
    err_ryx = err_ryx.reshape((-1,))
    err_pxy = err_pxy.reshape((-1,))
    err_pyx = err_pyx.reshape((-1,))
    #
    # Plotting section
    #
    plt.figure(num=1)
    plt.subplot(211)
    plt.scatter(ftlist,rho_xy,c='r',s=10)
    plt.scatter(ftlist,rho_yx,c='b',s=10)
    plt.errorbar(ftlist,rho_xy,err_rxy,ecolor='r',fmt="none")
    plt.errorbar(ftlist,rho_yx,err_ryx,ecolor='b',fmt="none")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim((10000, 0.001)) 
    plt.ylim(0.1, 100000)
    plt.yticks([0.1, 1, 10, 100, 1000, 10000])
    plt.grid(which='both',linestyle='-.', linewidth=0.4)
    plt.title(procinfo.get('selectedsite') + ' - ' + procinfo.get('meas')+' ('+str(procinfo.get('fs'))+' Hz)')
    plt.subplot(212)
    plt.scatter(ftlist,phase_xy,c='r',s=10)
    plt.scatter(ftlist,phase_yx,c='b',s=10)
    plt.errorbar(ftlist,phase_xy,err_pxy,ecolor='r',fmt="none")
    plt.errorbar(ftlist,phase_yx,err_pyx,ecolor='b',fmt="none")
    plt.xscale('log')
    plt.xlim((10000, 0.001))
    plt.ylim((0, 90))
    plt.yticks([0,15,30,45,60,75,90])
    plt.grid(which='both',linestyle='-.', linewidth=0.4)
    # plt.savefig('C:/Users/ajithabh/Desktop/MTProc/'+ procinfo.get('meas')+'.png', format='png', dpi=600)
    plt.figure(num=2)
    plt.scatter(ftlist,cohEx,c='r')
    plt.scatter(ftlist,cohEy,c='b')
    plt.xscale('log')
    plt.xlim((10000, 0.001)) 
    plt.ylim(0, 1)
    plt.yticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
    plt.grid(which='both',linestyle='-.', linewidth=0.4)
    plt.title(procinfo.get('selectedsite') + ' - ' + procinfo.get('meas')+' ('+str(procinfo.get('fs'))+' Hz)')
    # plt.savefig('C:/Users/ajithabh/Desktop/MTProc/'+ 'coh' + procinfo.get('meas')+'.png', format='png', dpi=600)
    # Plot tipper =====
    TxA = np.sqrt((np.real(Tx) ** 2) + (np.imag(Tx) ** 2))
    TyA = np.sqrt((np.real(Ty) ** 2) + (np.imag(Ty) ** 2))
    TxP = np.degrees(np.arctan2(Tx.imag,Tx.real))
    TyP = np.degrees(np.arctan2(Ty.imag,Ty.real))
    plt.figure(num=3)
    plt.subplot(211)
    plt.scatter(ftlist,TxA,c='r')
    plt.scatter(ftlist,TyA,c='b')
    plt.ylim(0, 1)
    plt.xscale('log')
    plt.xlim((10000, 0.001))
    plt.grid(which='both',linestyle='-.', linewidth=0.4)
    plt.title(procinfo.get('selectedsite') + ' - ' + procinfo.get('meas')+' ('+str(procinfo.get('fs'))+' Hz)')
    plt.subplot(212)
    plt.scatter(ftlist,TxP,c='r')
    plt.scatter(ftlist,TyP,c='b')
    plt.xscale('log')
    plt.xlim((10000, 0.001))
    plt.grid(which='both',linestyle='-.', linewidth=0.4)
    #----------
    plt.figure(num=4)
    plt.subplot(211)
    plt.scatter(ftlist,Tx.real,c='r')
    plt.scatter(ftlist,Ty.real,c='b')
    #plt.ylim(0, 1)
    plt.xscale('log')
    plt.xlim((10000, 0.001))
    plt.grid(which='both',linestyle='-.', linewidth=0.4)
    plt.title(procinfo.get('selectedsite') + ' - ' + procinfo.get('meas')+' ('+str(procinfo.get('fs'))+' Hz)')
    plt.subplot(212)
    plt.scatter(ftlist,Tx.imag,c='r')
    plt.scatter(ftlist,Ty.imag,c='b')
    plt.xscale('log')
    plt.xlim((10000, 0.001))
    plt.grid(which='both',linestyle='-.', linewidth=0.4)