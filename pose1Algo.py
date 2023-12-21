# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 21:40:22 2023

@author: premr
"""

def getMeanRotFromSum(rSum):
    #   rSum is 3x3 mat = simple algebraic sum / N
    #   Ravg is rotation matrix
    #   the averaged calculated as in Eq.(3.7) Moakher 2002 "Means & Averaging
    #   in the Group of Rotations', SIAM J.Matrix Analysis and Applications,
    #   vol.24, pp.1-16
    Rc = rSum.T@rSum;
    [U,D,Vh] = np.linalg.svd(Rc);
    dd = 1/np.sqrt(D)
    Lambda = np.diag(dd)
    Rbar = rSum@U@Lambda@U.T;
    return Rbar



import glob
import numpy as np
import matplotlib.pyplot as plt
folderAll  = ['.\Artifact1'];
sensors = ['Sensor1','Sensor2']

for ii in range(len(folderAll)):
    folder1 = folderAll[ii]
    filesAll = [];
    
    #First get all the filenames
    for fname in glob.glob(folder1+'\*.txt'):
        filesAll.append(fname)
        print(fname)
    if len(filesAll) == 0:
        print('No files were found - check path')
        exit()
    rAll = []

    #Now process all the files
    print('Processing all files ...')
    for fname1 in filesAll:
        dataT = np.genfromtxt(fname1,delimiter=",")
        data1 = dataT[:,:3] #Only get the xyz coordinates
        covMatrix = np.cov(data1.T)
        [r1, s1, vh1] = np.linalg.svd(covMatrix); 
        #v1 = vh1.T
        rAll.append(r1)
    Rsum = np.mean(rAll,axis=0);
    Rbar = getMeanRotFromSum(Rsum);
    
    # Peforms checks
    print('Check if the rotation matrix is orthogonal')
    checkOk = Rbar@Rbar.T - np.eye(3)
    err1 = np.median(np.abs(np.ndarray.flatten(checkOk)));
    if err1 < 1E-6:
        print('OK - Rotation matrix is orthogonal')
    else:
        print('Rotation matrix is NOT orthogonal')
        exit()
    
    
    #Calculate the angle alpha
    alpha1 = []; test1 = []
    for jj in range(len(filesAll)):
        rm = rAll[jj];
        deltaRM = Rbar@rm.T;
        test1.append(0.5*np.trace(deltaRM)-0.5)
        alphaTemp = np.arccos(0.5*np.trace(deltaRM)-0.5);
        alpha1.append(alphaTemp)
    
    U = 95.0;
    P1 = np.percentile(alpha1,95)
    MULT = 1000;
    print('Range of angles = %2.3f - %2.3f milliradian, %2.1f percentile = %2.3f milliradian\n' % (min(alpha1)*MULT,max(alpha1)*MULT, U, P1*MULT) )
    
    plt.figure(ii)
    hh = plt.hist(alpha1,20,edgecolor='black', linewidth=1.2)
    y1 = np.arange(np.max(hh[0]))
    x1 = y1*0+P1
    plt.plot(x1,y1)
    plt.legend(['95th percentile'])
    plt.xlabel('Angle in radians')
    plt.ylabel('Frequency')
    plt.title(sensors[ii]+': Histogram and '+str(U)+' percentile')
    pltName = 'PYY_'+sensors[ii]+'.png'
    plt.savefig(pltName)
    
    #break;