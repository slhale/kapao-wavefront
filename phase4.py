# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 22:50:49 2015
Edited on Thur Mar 3 2016

@author: Emily Yang
@editor: Sarah Hale
"""

# -*- coding: utf-8 -*-


import numpy as np
#import pyximport
#pyximport.install()
from kapaolibplus import Slopes
from FTR.utils import circle_aperture, remove_piston, remove_tiptilt
from FTR import FourierTransformReconstructor as FTRecon
import matplotlib.pyplot as plt

SLOPE_X_FILE = './runs/slope_x_' + '20140710_224918' + '.tel'
SLOPE_Y_FILE = './runs/slope_y_' + '20140710_224918' + '.tel'

#SLOPE_X_FILE = './run_421/slope_x_' + '20131217_023427' + '.tel'
#SLOPE_Y_FILE = './run_421/slope_y_' + '20131217_023427' + '.tel'

#SLOPE_X_FILE = './TelemetryOpenLoop/slope_x_' + filenum + '.tel'
#SLOPE_Y_FILE = './TelemetryOpenLoop/slope_y_' + filenum + '.tel'

slope_x, slope_y = Slopes(SLOPE_X_FILE), Slopes(SLOPE_Y_FILE)



def change_files():
    global SLOPE_X_FILE
    global SLOPE_Y_FILE
    global slope_x
    global slope_y
    
    # get user input
    print "Please input the date associated with the file"
    print "   for example, '20140710_224918'"
    date = raw_input(" > ")
    
    # if the user didn't actually put anything, then default to this
    if date == '':
        date = '20140710_224918'
   
    # redo the slope file names with this new date 
    SLOPE_X_FILE = './runs/slope_x_' + date + '.tel'
    SLOPE_Y_FILE = './runs/slope_y_' + date + '.tel'
    
    # redo the slopes with the new data  
    slope_x, slope_y = Slopes(SLOPE_X_FILE), Slopes(SLOPE_Y_FILE)


def slope_to_recon(slope_frame):
    """
    Turn a 1D list of slope values to an 11x11 grid with zeros in empty spots
    
    """
    grid = np.zeros((11,11))
    grid[0][3:8] = slope_frame[0:5]
    grid[1][2:9] = slope_frame[5:12]
    grid[2][1:10] = slope_frame[12:21]
    grid[3] = slope_frame[21:32]
    grid[4] = slope_frame[32:43]
    grid[5] = slope_frame[43:54]
    grid[6] = slope_frame[54:65]
    grid[7] = slope_frame[65:76]
    grid[8][1:10] = slope_frame[76:85]
    grid[9][2:9] = slope_frame[85:92]
    grid[10][3:8] = slope_frame[92:97]
    grid = grid.T # we're filling in row-by-row from the top, but numbering starts
                  # in the bottom left with zero and proceeds column-by-column
    return grid

'''
xs = slope_to_recon(slope_x.data[0])
ys = slope_to_recon(slope_y.data[0])
'''



shape = (11, 11)
r = 5.5
ap = circle_aperture(shape, r)


def remove_ttp(data):
    pr, p = remove_piston(ap, data)
    tr, tx, ty = remove_tiptilt(ap, pr)
    return tr


def recon2(timestep):
    xs = slope_to_recon(slope_x.data[timestep])
    ys = slope_to_recon(slope_y.data[timestep])
    recon = FTRecon(ap, filter='mod_hud', suppress_tt=True)
    phi = recon(xs, ys)
    
    reconflat = remove_ttp(phi).flatten()
    
    return reconflat


def rmsplot():
    rmsarray = np.zeros(0)    
    for i in xrange(len(slope_x.data)):
        rec = recon2(i)
        rms = np.sqrt(np.mean(np.square(rec)))
        rmsarray = np.append(rmsarray,rms)
    
    plt.plot(range(len(rmsarray)),rmsarray)
    plt.xlabel('Timestep')
    plt.ylabel('RMS')
    plt.show()
