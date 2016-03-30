###  To run this code on the AO Lab Linuxes:
##  1) Move to the directory with TV3.py in a terminal window.
##  2) type "ur_setup"
##  3) type "python"
##  4) type "from TV3 import telviz" at the python prompt
##  5a) then telviz('subdirectory','filenum') at the python prompt.
##  5b) then telviz('subdirectory','filenum',True) at the python prompt to send
##      plots to .pngs.
##
##  Both the subdirectory and the filenum should have single quotes around them
##  to make them a string. So as it is currently in the Dropbox, you'd type 
##  telviz('TelemetryOpenLoop', '20140621_000545') as an example.
##
##
##  To run on the lab Mac, skip step 2.

## Modified by Sarah Hale
## can email at shale@hmc.edu


def telviz(subdirectory,filenum,save=False):
    '''
    Visualizes slopes, intensities, and newpos data in a heatmap grid, and 
    plots tip/tilt and pinned actuators as a function of time.
    Inputs: subdirectory & filenum = strings, save = optional boolean
    Outputs: 1 figure with first and last time step visualizations of each tel
            file type, and 2 plots as a function of time
    '''
    
    
    from matplotlib import pyplot as plt
    import numpy as np
    import os.path
    
    from kapaolibplus import (subaps_to_grid, overlay_indices, Slopes, IntensityMap, DMPositions, newpos_to_grid, overlay_indices_newpos, slope_to_grid, overlay_indices_slope)
    
    
    
    filenum1, filenum2 = filenum.split('_')
    


    SLOPE_X_FILE = subdirectory + '/' + 'slope_x_' + filenum + '.tel'
    SLOPE_Y_FILE = subdirectory + '/' + 'slope_y_' + filenum + '.tel'
    INTENSITY_MAP = subdirectory + '/' + 'intensity_map_' + filenum + '.tel'
    NEWPOS_FILE = subdirectory + '/' + 'new_pos_' + filenum + '.tel'
    
    if os.path.isfile(SLOPE_X_FILE)==False:
        SLOPE_X_FILE = subdirectory + '/' + 'slope_x_' + filenum1 + '_' + str(int(filenum2)-1) + '.tel'
        if os.path.isfile(SLOPE_X_FILE) == False:
            SLOPE_X_FILE = subdirectory + '/' + 'slope_x_' + filenum1 + '_' + str(int(filenum2)+1) + '.tel'
    if os.path.isfile(SLOPE_Y_FILE)==False:
        SLOPE_Y_FILE = subdirectory + '/' + 'slope_y_' + filenum1 + '_' + str(int(filenum2)-1) + '.tel'
        if os.path.isfile(SLOPE_Y_FILE) == False:
            SLOPE_Y_FILE = subdirectory + '/' + 'slope_y_' + filenum1 + '_' + str(int(filenum2)+1) + '.tel'
    if os.path.isfile(INTENSITY_MAP)==False:
        INTENSITY_MAP = subdirectory + '/' + 'intensity_map_' + filenum1 + '_' + str(int(filenum2)-1) + '.tel'
        if os.path.isfile(INTENSITY_MAP) == False:
            INTENSITY_MAP = subdirectory + '/' + 'intensity_map_' + filenum1 + '_' + str(int(filenum2)+1) + '.tel'
    if os.path.isfile(NEWPOS_FILE)==False:
        NEWPOS_FILE = subdirectory + '/' + 'new_pos_' + filenum1 + '_' + str(int(filenum2)-1) + '.tel'
        if os.path.isfile(NEWPOS_FILE) == False:
            NEWPOS_FILE = subdirectory + '/' + 'new_pos_' + filenum1 + '_' + str(int(filenum2)+1) + '.tel'



#    SLOPE_X_FILE = './' + subdirectory + '/' + 'slope_x_' + filenum + '.tel'
#    SLOPE_Y_FILE = './' + subdirectory + '/' + 'slope_y_' + filenum + '.tel'
#    INTENSITY_MAP = './' + subdirectory + '/' + 'intensity_map_' + filenum + '.tel'
#    NEWPOS_FILE = './' + subdirectory + '/' + 'new_pos_' + filenum + '.tel'
#    
#    if os.path.isfile(SLOPE_X_FILE)==False:
#        SLOPE_X_FILE = './' + subdirectory + '/' + 'slope_x_' + filenum1 + '_' + str(int(filenum2)-1) + '.tel'
#        if os.path.isfile(SLOPE_X_FILE) == False:
#            SLOPE_X_FILE = './' + subdirectory + '/' + 'slope_x_' + filenum1 + '_' + str(int(filenum2)+1) + '.tel'
#    if os.path.isfile(SLOPE_Y_FILE)==False:
#        SLOPE_Y_FILE = './' + subdirectory + '/' + 'slope_y_' + filenum1 + '_' + str(int(filenum2)-1) + '.tel'
#        if os.path.isfile(SLOPE_Y_FILE) == False:
#            SLOPE_Y_FILE = './' + subdirectory + '/' + 'slope_y_' + filenum1 + '_' + str(int(filenum2)+1) + '.tel'
#    if os.path.isfile(INTENSITY_MAP)==False:
#        INTENSITY_MAP = './' + subdirectory + '/' + 'intensity_map_' + filenum1 + '_' + str(int(filenum2)-1) + '.tel'
#        if os.path.isfile(INTENSITY_MAP) == False:
#            INTENSITY_MAP = './' + subdirectory + '/' + 'intensity_map_' + filenum1 + '_' + str(int(filenum2)+1) + '.tel'
#    if os.path.isfile(NEWPOS_FILE)==False:
#        NEWPOS_FILE = './' + subdirectory + '/' + 'new_pos_' + filenum1 + '_' + str(int(filenum2)-1) + '.tel'
#        if os.path.isfile(NEWPOS_FILE) == False:
#            NEWPOS_FILE = './' + subdirectory + '/' + 'new_pos_' + filenum1 + '_' + str(int(filenum2)+1) + '.tel'
        
        


    slope_x, slope_y, intensity_map = Slopes(SLOPE_X_FILE), Slopes(SLOPE_Y_FILE), IntensityMap(INTENSITY_MAP)
    new_pos = DMPositions(NEWPOS_FILE)
    
 
    lenx = len(slope_x.data)
    leny = len(slope_y.data)
    leni = len(intensity_map.data)
    lenn = len(new_pos.data)
    

    
    # Plots

     # Summary Panel
        
    figall = plt.figure(figsize=(12,8))#(10,7.5))#(20,15))

    # Set font size to be smaller 
    plt.rcParams.update({'font.size': 10})
 
    plt.subplot2grid((3,2),(0,0))#, colspan=2)
    plt.imshow(slope_to_grid(slope_x.data[0]), origin='lower')
    #overlay_indices_slope()
    plt.colorbar()
    plt.title('X Slope - first time step')
   
    ''' 
    plt.subplot2grid((3,4),(0,1))
    plt.imshow(slope_to_grid(slope_x.data[lenx-1]), origin='lower')
    overlay_indices_slope()
    plt.colorbar()
    plt.title('X Slope - last time step')
    '''

    plt.subplot2grid((3,2),(0,1))#2), colspan=2)
    plt.imshow(slope_to_grid(slope_y.data[0]), origin='lower')
    #overlay_indices_slope()
    plt.colorbar()
    plt.title('Y Slope - first time step')

    '''
    plt.subplot2grid((3,4),(0,3))
    plt.imshow(slope_to_grid(slope_y.data[leny-1]), origin='lower')
    overlay_indices_slope()
    plt.colorbar()
    plt.title('Y Slope - last time step')
    '''

    plt.subplot2grid((3,2),(1,0))#, colspan=2)
    plt.imshow(newpos_to_grid(new_pos.data[0]), origin='lower')
    #overlay_indices_newpos()
    plt.colorbar()
    plt.title('DM Pos - first time step')

    '''    
    plt.subplot2grid((3,4),(1,1))
    plt.imshow(newpos_to_grid(new_pos.data[lenn-1]), origin='lower')
    overlay_indices_newpos()
    plt.colorbar()
    plt.title('DM Pos - last time step')
    '''

    plt.subplot2grid((3,2),(1,1))#2), colspan=2)
    plt.imshow(subaps_to_grid(intensity_map.data[0]), origin='lower')
    #overlay_indices()
    plt.colorbar()
    plt.title('Intensity - first time step')

    '''    
    plt.subplot2grid((3,4),(1,3))
    plt.imshow(subaps_to_grid(intensity_map.data[leni-1]), origin='lower')
    overlay_indices()
    plt.colorbar()
    plt.title('Intensity - last time step')
    '''

    
    #Tip/tilt as a function of time
    
    #pulling out the 120th and 122nd index for each time step
    tt_1 = np.zeros(lenn)
    tt_2 = np.zeros(lenn)
    for i in range(0,lenn):
        tt_1[i] = new_pos.data[i][120] # Channel 2 (Left/Right on PDVShow & Andor)
        tt_2[i] = new_pos.data[i][122] # Channel 1 (Up/Down on PDVShow & Andor)
        


    plt.subplot2grid((3,2),(2,0))#,colspan=2)
    plt.plot(new_pos.timestamps - new_pos.timestamps[0],tt_1,'.', new_pos.timestamps - new_pos.timestamps[0],tt_2,'.')
    plt.ylabel("Tip/tilt")
    plt.title('Tip/Tilt as  Function of Time')
    plt.xlabel("Time (ms)")
    plt.legend(['120 (L/R? PDVShow)', '122 (U/D? PDVShow)'],'upper left')    
        
    #pinned actuators as a function of time
    
    pinned = np.zeros(lenn)
    for i in range(0,lenn):
        for j in range(0,120):
            if new_pos.data[i][j] <= 100 or new_pos.data[i][j] >= 64900:
                pinned[i] = pinned[i] + 1
                
    plt.subplot2grid((3,2),(2,1))#2),colspan=2)
    plt.plot(new_pos.timestamps - new_pos.timestamps[0],pinned,'.')
    plt.ylabel("Number of pinned actuators")
    plt.xlabel("Time (ms)")
    plt.title("Pinned Actuators as a Function of Time")


    plt.suptitle('All, Unfixed ('+filenum+')',fontsize=16)        

    if save == True:
        figall.savefig('./' + subdirectory + '/' + 'fig_all_' + filenum + '.png', dpi=300)

    
    plt.show()
