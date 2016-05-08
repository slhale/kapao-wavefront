###  To run this code on the AO Lab Linuxes:
##  1) Move to the directory with TV5.py in a terminal window.
##  2) type "ur_setup"
##  3) type "python"
##  4) type "from TV5 import telviz" at the python prompt
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


# default values for subdirectory and filenum match Sarah's setup 
def telviz(subdirectory='runs',filenum='20140701_214748',normalize_recon=True,normalize_dm=True,normalize_intensity=True,run_summary=False,save=False):
    '''
    Visualizes slopes, intensities, and newpos data in a heatmap grid, and 
    plots tip/tilt and pinned actuators as a function of time.
    Inputs: subdirectory & filenum = strings, save = optional boolean
    Outputs: 1 figure with first and last time step visualizations of each tel
            file type, and 2 plots as a function of time
    '''
    
    
    from matplotlib import pyplot as plt
    from matplotlib.widgets import Slider
    from matplotlib.font_manager import FontProperties
    import numpy as np
    from FTR import FourierTransformReconstructor as FTRecon
    from FTR.utils import circle_aperture, remove_piston, remove_tiptilt
    import os.path
    import imp
    from kapaolibplus import (subaps_to_grid, overlay_indices, Slopes, IntensityMap, DMPositions, newpos_to_grid, overlay_indices_newpos, slope_to_grid, overlay_indices_slope)
    
    # If we are told to run the ConfigSummary, do so
    if run_summary:
        # Non-sketch way does not work :(
        #summary = imp.load_source('ConfigSummary', '/home/student/Pomona_v3/Config/ConfigSummary.py')
        #summary.main() 
        # Use sketch way
        execfile('/home/student/Pomona_v3/Config/ConfigSummary.py')
    
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
    ## Copied from Sarah's modified phase5.py
    
    shape = (11, 11)
    r = 5.5
    ap = circle_aperture(shape, r)
    
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
        recon = FTRecon(ap, filter='mod_hud', suppress_tt=True)
        phi = recon(xs, ys)
        
        return phi
    
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
    
    
    lenx = len(slope_x.data)
    leny = len(slope_y.data)
    leni = len(intensity_map.data)
    lenn = len(new_pos.data)
    
    # Go through the slopes data and recon for each timestep (and store the recon)
    wave_recon = []
    intensity = []
    new_dm = []
    max_time = len(slope_x.data) - 1 # the maximum index that exists for the time 
    for t in range(max_time):
        wave_recon.append(recon(t))
        # Go through the intensity data and new pos data and store the processed version too
        intensity.append(subaps_to_grid(intensity_map.data[t]))
        new_dm.append(newpos_to_grid(new_pos.data[t]))
 
    # Find the min and max values for each heatmap plot 
    # do this in order to get the correct normalization for the colorbars
    #print wave_recon[4][1] # prints an array
    #print intensity[10][1] # prints an array
    #print new_dm[10][1] # prints an array
    # the wave reconstruction is an array of 2d arrays
    # the intensity is an array of arrays
    min_recon_val = 1000000  # big number
    max_recon_val = -1000000 # small number
    max_slope_x_val = -1000000
    min_dm_val = 1000000
    max_dm_val = -1000000
    min_intensity_val = 1000000
    max_intensity_val = -1000000

    print 'max time', max_time
    #print 'len of intensity', len(intensity)

    # loop over all the times
    for t in range(max_time):
        if normalize_recon: 
            # loop over all the individual data values in the wave reconstruction
            wave_max_i = len(wave_recon[0])
            for i in range(wave_max_i):
                wave_max_j = len(wave_recon[0][0])
                for j in range(wave_max_j):
                    current_recon_val = wave_recon[t][i][j]
                    if min_recon_val > current_recon_val:
                        min_recon_val = current_recon_val
                    if max_recon_val < current_recon_val:
                        max_recon_val = current_recon_val
                #
            #
        if normalize_intensity:
            # loop over all the individual data values in the intensity
            inten_max_i = len(intensity[0])
            for i in range(inten_max_i):
                inten_max_j = len(intensity[0][0])
                for j in range(inten_max_j):
                    current_inten_val = intensity[t][i][j]
                    if min_intensity_val > current_inten_val:
                        min_intensity_val = current_inten_val
                    if max_intensity_val < current_inten_val:
                        max_intensity_val = current_inten_val
                #
            #
        if normalize_dm:
            # loop over all the individual data values in the new dm position
            dm_max_i = len(new_dm[0])
            for i in range(dm_max_i):
                dm_max_j = len(new_dm[0][0])
                for j in range(dm_max_j):
                    current_dm_val = new_dm[t][i][j]
                    if min_dm_val > current_dm_val:
                        min_dm_val = current_dm_val
                    if max_dm_val < current_dm_val:
                        max_dm_val = current_dm_val
                #
            #
    
    #print "wave recon slice", wave_recon[0]
    print "max recon val", max_recon_val
    print "min recon val", min_recon_val
    #print "intensity slice", intensity[0]
    print "max intensity val", max_intensity_val
    print "min intensity val", min_intensity_val
    #print "new dm slice", new_dm[0]
    print "max dm val", max_dm_val
    print "min dm val", min_dm_val
    
    # Make initial data so that the colorbars are ok
    recon_init = wave_recon[0]
    if normalize_recon:
        recon_init[0][0] = max_recon_val
        recon_init[0][1] = min_recon_val
    dm_init = new_dm[0]
    if normalize_dm:
        dm_init[0][0] = max_dm_val
        dm_init[0][1] = min_dm_val
    inten_init = intensity[0]
    if normalize_intensity:
        inten_init[0][0] = max_intensity_val
        inten_init[0][1] = min_intensity_val

    # Initialize our data variables
    # the variables will be updated when the time slider changes
    #slope_x_data = slope_x.data[0]
    #slope_y_data = slope_y.data[0]
    recon_data = recon_init#wave_recon[0]
    new_pos_data = dm_init#newpos_to_grid(new_pos.data[0])
    intensity_map_data = inten_init#subaps_to_grid(intensity_map.data[0])
    
    # Plots

     # Summary Panel
        
    figall = plt.figure(figsize=(10,8))
    
    # Set font size to be smaller 
    plt.rcParams.update({'font.size': 10})

    # Set the subplots to have spacing between them 
    plt.subplots_adjust(hspace=0.4, wspace=0.4)

    # Plot the three 'square' plots vertically on the left side 
 
    plt.subplot2grid((3,2),(0,0))
    intensity_map_im = plt.imshow(intensity_map_data, origin='lower', interpolation='none')
    plt.colorbar()
    plt.title('Intensity')
    
    plt.subplot2grid((3,2),(1,0))
    new_pos_im = plt.imshow(new_pos_data, origin='lower', interpolation='none')
    plt.colorbar()
    plt.title('DM Position')

    plt.subplot2grid((3,2),(2,0))
    recon_im = plt.imshow(recon_data, origin='lower', interpolation='none')
    plt.colorbar()
    plt.title('Wave Reconstruction')

    # Add axes for time slider
    axes = figall.add_axes([0.25, 0.02, 0.5, 0.02])
    timeslider = Slider(axes, 'Time', 0, max_time, valinit=0, valfmt='%i')


    def update(val):
        # Update the data
        time_index = int(val)
        recon_data = recon(time_index)#wave_recon[time_index]
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
    
    # Plot the three 'rectangular' plots vertically on the right side 
    
    ## RMS as a funtion of time plot
    
    rms_x_array = np.zeros(0)    
    for i in xrange(len(slope_x.data)):
        rec = recon2(i)
        rms = np.sqrt(np.mean(np.square(rec)))
        rms_x_array = np.append(rms_x_array,rms)
    
    rms_y_array = np.zeros(0)    
    for i in xrange(len(slope_y.data)):
        rec = recon2(i)
        rms = np.sqrt(np.mean(np.square(rec)))
        rms_y_array = np.append(rms_y_array,rms)
    
    plt.subplot2grid((3,2),(0,1))
    #plt.plot(range(len(rms_x_array)),rms_x_array, 'b')
    #plt.plot(range(len(rms_y_array)),rms_y_array, 'g')
    plt.plot(slope_x.timestamps - slope_x.timestamps[0], rms_x_array, 'b')
    plt.plot(slope_y.timestamps - slope_y.timestamps[0], rms_y_array, 'g')
    #plt.xlabel('Timestep')
    plt.xlabel('Time (ms)')
    plt.ylabel('RMS')
    plt.title('RMS a Function of Time')
    
    
    ## Tip/tilt as a function of time plot
    
    # pulling out the 120th and 122nd index for each time step
    tt_1 = np.zeros(lenn)
    tt_2 = np.zeros(lenn)
    for i in range(0,lenn):
        tt_1[i] = new_pos.data[i][120] # Channel 2 (Left/Right on PDVShow & Andor)
        tt_2[i] = new_pos.data[i][122] # Channel 1 (Up/Down on PDVShow & Andor)
    
    plt.subplot2grid((3,2),(1,1))
    plt.plot(new_pos.timestamps - new_pos.timestamps[0],tt_1,'.', new_pos.timestamps - new_pos.timestamps[0],tt_2,'.')
    plt.ylabel("Tip/tilt")
    plt.title('Tip/Tilt as  Function of Time')
    plt.xlabel("Time (ms)")
    fontP = FontProperties()
    fontP.set_size('x-small')
    #plt.legend(['120 (L/R? PDVShow)', '122 (U/D? PDVShow)'],'best', prop=fontP)    
    lgn = plt.legend(['120 (L/R?)', '122 (U/D?)'],'best', prop=fontP)
    lgn.draggable()

    ## Pinned actuators as a function of time plot
    
    pinned = np.zeros(lenn)
    for i in range(0,lenn):
        for j in range(0,120):
            if new_pos.data[i][j] <= 100 or new_pos.data[i][j] >= 64900:
                pinned[i] = pinned[i] + 1
                
    plt.subplot2grid((3,2),(2,1))
    plt.plot(new_pos.timestamps - new_pos.timestamps[0],pinned,'.')
    plt.ylabel("Number of pinned actuators")
    plt.xlabel("Time (ms)")
    plt.title("Pinned Actuators as a Function of Time")


    plt.suptitle('All, Unfixed ('+filenum+')',fontsize=16)        

    if save == True:
        figall.savefig('./' + subdirectory + '/' + 'fig_all_' + filenum + '.png', dpi=300)
    
    update(0)
    
    plt.show()
