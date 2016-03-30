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


# default values for subdirectory nd filenum match Sarah's setup 
def telviz(subdirectory='runs',filenum='20140701_214748',save=False):
    '''
    Visualizes slopes, intensities, and newpos data in a heatmap grid, and 
    plots tip/tilt and pinned actuators as a function of time.
    Inputs: subdirectory & filenum = strings, save = optional boolean
    Outputs: 1 figure with first and last time step visualizations of each tel
            file type, and 2 plots as a function of time
    '''
    
    
    from matplotlib import pyplot as plt
    from matplotlib.widgets import Slider
    import numpy as np
    from FTR import FourierTransformReconstructor as FTRecon
    from FTR.utils import circle_aperture
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
        

    slope_x, slope_y, intensity_map = Slopes(SLOPE_X_FILE), Slopes(SLOPE_Y_FILE), IntensityMap(INTENSITY_MAP)
    new_pos = DMPositions(NEWPOS_FILE)
    
    
    ### Helper functions here 
    ## Copied from Sarah's modified phase4.py
    
    def slope_to_recon(slope_frame=slope_x.data[0]):
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
        # transpose the grid
        grid = grid.T # we're filling in row-by-row from the top, but numbering starts
                      # in the bottom left with zero and proceeds column-by-column
        return grid
    
    def recon(timestep):
        '''
            Reconstructs the wavefront at a given timestamp using x and y
            slope data. 
        '''
        xs = slope_to_recon(slope_x.data[timestep])
        ys = slope_to_recon(slope_y.data[timestep])
        shape = (11, 11)
        r = 5.5
        ap = circle_aperture(shape, r)
        recon = FTRecon(ap, filter='mod_hud', suppress_tt=True)
        phi = recon(xs, ys)
        
        return phi
    
    
    lenx = len(slope_x.data)
    leny = len(slope_y.data)
    leni = len(intensity_map.data)
    lenn = len(new_pos.data)
    
    # Initialize our data variables
    # the variables will be updated when the time slider changes
    slope_x_data = slope_x.data[0]
    slope_y_data = slope_y.data[0]
    recon_data = recon(0)
    new_pos_data = newpos_to_grid(new_pos.data[0])
    intensity_map_data = subaps_to_grid(intensity_map.data[0])
    
    # Plots

     # Summary Panel
        
    figall = plt.figure(figsize=(12,8))

    # Set font size to be smaller 
    plt.rcParams.update({'font.size': 10})
 
    plt.subplot2grid((3,2),(0,0))
    recon_im = plt.imshow(recon_data, origin='lower', interpolation='none')
    plt.colorbar()
    plt.title('Wave Reconstruction')
    
    '''
    plt.subplot2grid((3,2),(0,0))
    slope_x_im = plt.imshow(slope_to_grid(slope_x_data), origin='lower', interpolation='none')
    plt.colorbar()
    plt.title('X Slope - first time step')
   
    plt.subplot2grid((3,2),(0,1))
    slope_y_im = plt.imshow(slope_to_grid(slope_y_data), origin='lower', interpolation='none')
    plt.colorbar()
    plt.title('Y Slope - first time step')
    '''

    plt.subplot2grid((3,2),(1,0))
    new_pos_im = plt.imshow(new_pos_data, origin='lower', interpolation='none')
    plt.colorbar()
    plt.title('DM Pos - first time step')

    plt.subplot2grid((3,2),(1,1))
    intensity_map_im = plt.imshow(intensity_map_data, origin='lower', interpolation='none')
    plt.colorbar()
    plt.title('Intensity - first time step')

    # Add axes for time slider
    axes = figall.add_axes([0.25, 0.02, 0.5, 0.02])
    max_time = len(slope_x.data) - 1 # the maximum index that exists for the time 
    timeslider = Slider(axes, 'Time', 0, max_time, valinit=0, valfmt='%i')


    def update(val):
        # Update the data
        time_index = int(val)
        recon_data = recon(time_index)
        new_pos_data = newpos_to_grid(new_pos.data[time_index])
        intensity_map_data = subaps_to_grid(intensity_map.data[time_index])
        
        # Set the image array to this
        recon_im.set_array(recon_data)
        new_pos_im.set_array(new_pos_data)
        intensity_map_im.set_array(intensity_map_data)
        
        # Redraw the plot
        figall.canvas.draw()
        
    # Whe the slider is slid, update the plot
    timeslider.on_changed(update)
   
 
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
