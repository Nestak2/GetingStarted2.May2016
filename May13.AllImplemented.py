# What are the variables we want set as default?
# Taking out the user input leaves much less comparison options? what do we want to compare actually?
# Now bin_width=500, wanted_bin=0 occure as kew words in main() and in binning() function. Is this OK?


import GAS
from spectral_cube import SpectralCube
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as nd
import astropy.units as u


files =['NGC1333_NH3_11_DR1.fits', 'NGC1333_NH3_22_DR1.fits', 'NGC1333_NH3_33_DR1.fits', 
       'NGC1333_C2S_DR1.fits',
       'NGC1333_HC5N_DR1.fits',
        'NGC1333_HC7N_21_20_DR1.fits', 'NGC1333_HC7N_22_21_DR1.fits']


def main(bin_width=500, thisbin=0):   
    # This is the main routine.
    loop_count = 0
    fig = plt.figure(figsize=(10,7))
    # for plotting plots in a single figure put plt.figure before all the .plot commands
    # otherwise you will get the plots in a different figures
    
    ax1 = plt.subplot(len(files), 1, 1)
    
    f_nam = 'NGC1333_NH3_11_DR1.fits'
    y, x = binning(f_nam)
    # This takes in a 2D map and returns a 2D map where the image values are the bin to which a pixel belongs.
    
    for file_name in files:            
        print 'in loop', loop_count, 'the values of y[100:105] and x[100:105] are:', y[100:105], x[100:105]
        sp_av, cube = averaging(file_name, y, x)
        #Change to velocity axis and such here as well.  Then average all spectra with that bin label.
        ax = fig.add_subplot(len(files),1,loop_count+1)
        plt.plot(cube.spectral_axis,sp_av)
        plt.yticks(np.arange(min(sp_av), max(sp_av), (max(sp_av)-min(sp_av))/4))
        
        if file_name != files[-1]:
            plt.setp(ax.get_xticklabels(), visible=False)

        plt.legend([file_name], bbox_to_anchor=(1.05, 1.1), prop={'size':8})
        loop_count += 1        

#     plt.legend(leg, loc = 5, prop={'size':8})
#     ax.legend(bbox_to_anchor=(1.05, 0), loc='lower left', borderaxespad=0.)
#     fig.legend(curves_names, leg, bbox_to_anchor=(0.5, -0.15))
#    plt.savefig("May13.AllImplementedBin_width=%rThisbin=%r.png" %(bin_width, thisbin))    
#     figtext(.0,.0,'From top to botom:\n' + str(files) , fontsize=8)
    fig.suptitle("average spectrum of pixels in brightnes bin = %r for bin width = %r" %(thisbin, bin_width), fontsize=12)
    plt.show()

def binning(f_nam, bin_width=500, thisbin=0):
    """A function creating brightness bins of pixels, and eventualy a map, in the given spectral cube"""
    cube = SpectralCube.read(f_nam)
    cube = cube.with_spectral_unit(u.km/u.s,velocity_convention='radio')
    Tmax = cube.apply_numpy_function(np.nanmax,axis=0) # array of the maximum values in the spectra of each pixel
    baddata = nd.morphology.binary_dilation(np.isnan(Tmax),np.ones((25,25)))
    Tmax[baddata]=0.0
    Tmax[np.isfinite(Tmax)]
    bin_arr = np.sort(Tmax[np.isfinite(Tmax)])
    bin_arr2 = bin_arr[:: - bin_width] # this creates an array of the bin margins, in which every bin has a width of "bin_width"  
    np.digitize(Tmax,bin_arr2)
    bins = np.digitize(Tmax,bin_arr2)
    y, x = np.where(bins==thisbin)
    return y, x

def averaging(file_name, y, x):
    cube = SpectralCube.read(file_name)
    cube = cube.with_spectral_unit(u.km/u.s,velocity_convention='radio')
    sp_av = np.nanmean(cube.filled_data[:,y,x].value,axis=1)
    return sp_av, cube

main()
