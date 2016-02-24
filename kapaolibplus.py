import re
import glob
import datetime
import dateutil
import dateutil.parser
#import pyfits
import astropy.io.fits as pyfits

from matplotlib import pyplot
from matplotlib import pyplot as plt
import numpy
import numpy as np

UTC_TZ = dateutil.tz.gettz('UTC')

TOTAL_SUBAPS = 97
TOTAL_XSLOPES = 97
TOTAL_YSLOPES = 97
TOP_N_SUBAPS = 50
TOTAL_NEWPOS = 120

def subaps_to_grid(subap_frame, plot=False):
    """
    Turn a 1D list of subaperture values to an 11x11 grid
    
    See figure from Rudy / Contreras for how we count our subaps.
    
    Or just do this:
    
    >>> subaps_to_grid(range(97), plot=True)
    >>> overlay_indices()
    
    (n.b. sensible image plotting requires the
    origin='bottom' keyword argument unless you set
    it as default in matplotlibrc)
    """
    grid = np.ndarray((11,11))
    grid[:,:] = np.nan
    grid[0][3:8] = subap_frame[0:5]
    grid[1][2:9] = subap_frame[5:12]
    grid[2][1:10] = subap_frame[12:21]
    grid[3] = subap_frame[21:32]
    grid[4] = subap_frame[32:43]
    grid[5] = subap_frame[43:54]
    grid[6] = subap_frame[54:65]
    grid[7] = subap_frame[65:76]
    grid[8][1:10] = subap_frame[76:85]
    grid[9][2:9] = subap_frame[85:92]
    grid[10][3:8] = subap_frame[92:97]
    grid = grid.T # we're filling in row-by-row from the top, but numbering starts
                  # in the bottom left with zero and proceeds column-by-column
    if plot:
        plt.imshow(grid, origin='bottom')
    return grid

from matplotlib import patheffects

def overlay_indices():
    '''
    Number subapertures 0..96 on an image by overlaying text
    '''
    labels = subaps_to_grid(range(97))
    for col in range(11):
        for row in range(11):
            if not np.isnan(labels[row][col]):
                plt.text(
                     col - 0.25, row - 0.25, int(labels[row][col]),
                     color='white',
                     fontsize=10,
                     path_effects=[patheffects.withStroke(linewidth=3, foreground="black")]
                )

def overlay_subap_values(values):
    '''
    Overlay 97 values on their respective subapertures
    by drawing text in the appropriate spots
    '''
    labels = subaps_to_grid(values)
    for col in range(11):
        for row in range(11):
            if not np.isnan(labels[row][col]):
                plt.text(
                     col - 0.25, row - 0.25, labels[row][col],
                     color='white',
                     fontsize=10,
                     path_effects=[patheffects.withStroke(linewidth=3, foreground="black")]
                )

def parse_telem_table(filename):
    dates = []
    data_frames = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.split('\t')
            
            dt = dateutil.parser.parse(line[0])
            # dt = datetime.datetime.strptime(line[0], '%Y-%m-%d %H:%M:%S.%f')
            dt = dt.replace(tzinfo=UTC_TZ) # be very explicit about wanting UTC on every platform
            dt = np.datetime64(dt) # numpy prints local-tz dates/times, but handles tz-aware datetimes correctly
            
            data = line[1:-1] # chop off date from front, \n from end
            data = map(float, data)
            dates.append(dt)
            data_frames.append(data)
    
    timestamps = numpy.array(dates, dtype='datetime64[ms]')
    # timestamps not guaranteed to be exactly equal, so we take the 99th percentile
    # then to get seconds, divide by a timedelta of 1 s
    dt = np.percentile(np.diff(timestamps), 99) / np.timedelta64(1, 's')
    data = numpy.array(data_frames)
    return data, timestamps, dt

class Slopes(object):
    '''
    Slopes - A collection of (x or y) slope values loaded
    from a telemetry file (named as the first argument to __init__).
    '''
    def __init__(self, filename):
        self.filename = filename
        self.data, self.timestamps, self.dt = parse_telem_table(filename)
        self.median_by_slope = np.median(self.data, axis=0)#EYang
        
        self.overall_slopes = np.average(self.data, axis=1)
        self.subtracted_data = self.data - self.overall_slopes[:,np.newaxis]
    
    def do_fft(self, data=None, subtracted_data=None, fromidx=0, toidx=None):
        if data is None:
            data = self.data
        if subtracted_data is None:
            subtracted_data = self.subtracted_data
        
        freqs = np.fft.fftfreq(len(self.timestamps[fromidx:toidx]), self.dt)
        power = np.abs(np.fft.fft(self.data[fromidx:toidx], axis=0))
        subtracted_power = np.abs(np.fft.fft(self.subtracted_data[fromidx:toidx], axis=0))
        return freqs, power, subtracted_power
    
    def show_fft(self, *args, **kwargs):
        freqs, power, subtracted_power = self.do_fft(*args, **kwargs)
        power_reference = (freqs)**(-2.0/3.0)
        fig = plt.figure(figsize=(16,8))
        plt.plot(freqs, power_reference, 'g--', label="$f^{-2/3}$")
        plt.plot(freqs, np.sum(power, axis=1), label="Power (summed over subaps)")
        plt.plot(freqs, np.sum(subtracted_power, axis=1), label="Tilt-subtracted Power")
        plt.legend()
        plt.yscale('log')
        plt.xscale('log', nonposx='mask')
        plt.xlim(0.01, np.max(freqs))
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Power')
        return fig
    
    @property
    def freqs(self):
        if not hasattr(self, '_freqs'):
            self._freqs, self._power, self._subtracted_power = self.do_fft()
        return self._freqs
    
    @property
    def power(self):
        if not hasattr(self, '_power'):
            self._freqs, self._power, self._subtracted_power = self.do_fft()
        return self._power
    
    @property
    def subtracted_power(self):
        if not hasattr(self, '_subtracted_power'):
            self._freqs, self._power, self._subtracted_power = self.do_fft()
        return self._subtracted_power
    
    def __unicode__(self):
        return u'<Slopes: {0}>'.format(self.filename)
    __repr__ = __str__ = __unicode__

class IntensityMap(object):
    '''
    IntensityMap - A collection of per-subap intensities loaded
    from a telemetry file (named as the first argument to __init__).
    '''
    def __init__(self, filename, min_light=None):
        self.filename = filename
        self.min_light = min_light
        
        self.data, self.timestamps, self.dt = parse_telem_table(filename)
        self.median_by_subap = np.median(self.data, axis=0)
        if not min_light:
            print 'Guessing min_light from median intensity values ( taking top', TOP_N_SUBAPS, ' subaps)'
            median_sort = self.median_by_subap.copy()
            median_sort.sort()
            self.min_light = median_sort[TOP_N_SUBAPS - 1] # get median intensity for Nth-most-intense subap
            print 'min_light =', self.min_light
        
        self.include = set()
        for idx, flag in enumerate(self.median_by_subap > self.min_light):
            if flag:
                self.include.add(idx)
    
    @property
    def include_mask(self):
        mask = np.zeros((TOTAL_SUBAPS,), dtype='bool')
        for idx in self.include:
            mask[idx] = True
        return mask
    
    def show_included_subaps(self):
        plt.imshow(subaps_to_grid(self.include_mask), cmap='Greens', origin='bottom')
        overlay_indices()
    
    def show_frame(self, frame):
        return plt.imshow(subaps_to_grid(self.data[frame]), cmap='gray', origin='bottom')

    def write_FITS(self, filename):
        dest = np.zeros((self.data.shape[0], 11, 11))
        for i in range(self.data.shape[0]):
            dest[i] = subaps_to_grid(self.data[i])
        pyfits.writeto(filename, dest)
        print 'wrote to', filename

    def __unicode__(self):
        return u'<IntensityMap: {0}>'.format(self.filename)
    __repr__ = __str__ = __unicode__

class DMPositions(object):
    '''
    DMPositions - Actuator positions computed in each loop iteration
    (n.b. this is for use with `new_pos_` files, not `dm_pos_` which
    includes additional unnecessary data)
    '''
    def __init__(self, filename):
        self.filename = filename
        
        self.data, self.timestamps, self.dt = parse_telem_table(filename)
        self.median_by_pos = np.median(self.data, axis=0) ##### EY

    def __unicode__(self):
        return u'<DMPositions: {0}>'.format(self.filename)
    __repr__ = __str__ = __unicode__

def load_reconstructor(filename):
    reconstructor = []
    with open(filename, 'r') as f:
        rows, cols = map(int, f.readline().split('\t'))
        for line in f:
            line = line.strip()
            split_line = map(float, map(str.strip, line.split('\t')))
            reconstructor.append(split_line)
    return np.array(reconstructor)




'''
Emily Yang's additions to Joseph Long's kapaolib file
'''

#Turn newpos data into a heatmap grid

def newpos_to_grid(newpos_frame, plot=False):
    """
    Turn a 1D list of newpos values to a 12x12 grid
    
    """
    grid = np.ndarray((12,12))
    grid[:,:] = np.nan
    grid[0][3:9] = newpos_frame[0:6]
    grid[1][2:10] = newpos_frame[6:14]
    grid[2][1:11] = newpos_frame[14:24]
    grid[3] = newpos_frame[24:36]
    grid[4] = newpos_frame[36:48]
    grid[5] = newpos_frame[48:60]
    grid[6] = newpos_frame[60:72]
    grid[7] = newpos_frame[72:84]
    grid[8] = newpos_frame[84:96]
    grid[9][1:11] = newpos_frame[96:106]
    grid[10][2:10] = newpos_frame[106:114]
    grid[11][3:9] = newpos_frame[114:120]
    grid = grid.T # we're filling in row-by-row from the top, but numbering starts
                  # in the bottom left with zero and proceeds column-by-column
    if plot:
        plt.imshow(grid, origin='bottom')
    return grid


def overlay_indices_newpos():
    '''
    Number actuator positions 0..119 on an image by overlaying text
    '''
    labels = newpos_to_grid(range(120))
    for col in range(12):
        for row in range(12):
            if not np.isnan(labels[row][col]):
                plt.text(
                     col - 0.2, row - 0.2, int(labels[row][col]),
                     color='white',
                     fontsize=8,
                     path_effects=[patheffects.withStroke(linewidth=3, foreground="black")]
                )



#Turn slope data into a heatmap grid (basically same as subaps_to_grid, renamed
#to reduce confusion)


def slope_to_grid(slope_frame, plot=False):
    """
    Turn a 1D list of slope values to an 11x11 grid
    
    """
    grid = np.ndarray((11,11))
    grid[:,:] = np.nan
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
    if plot:
        plt.imshow(grid, origin='bottom')
    return grid


def overlay_indices_slope():
    '''
    Number slopes 0..96 on an image by overlaying text
    '''
    labels = slope_to_grid(range(97))
    for col in range(11):
        for row in range(11):
            if not np.isnan(labels[row][col]):
                plt.text(
                     col - 0.25, row - 0.25, int(labels[row][col]),
                     color='white',
                     fontsize=10,
                     path_effects=[patheffects.withStroke(linewidth=3, foreground="black")]
                )

