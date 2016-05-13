# What are the variables we want set as default?
# Taking out the user input leaves much less comparison options? what do we want to compare actually?


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


def main():   
    loop_count = 0
    leg = []
    curves_names = []
    fig = plt.figure(figsize=(10,7))
    # for plotting plots in a single figure put plt.figure before all the .plot commands
    # otherwise you will get the plots in a different figures
    
    ax1 = plt.subplot(len(files), 1, 1)
    
    for loop_count in range(len(files)):
        if loop_count == 0:
            _, _, _, _, _, _, _, y, x = binning(loop_count)
            # this creates the initial brightnes-bins map. For that the first name in 'files' is used. So if you want to
            # plot a map regarding to a different file, put this file on position 0 in 'files'
            
            
        Tmax, _, bin_arr2, cube, file_name, bin_width, wanted_bin, y, x = binning(loop_count)

        sp_av = np.nanmean(cube.filled_data[:,y,x].value,axis=1)
        ax = fig.add_subplot(len(files),1,loop_count+1)
        plt.plot(cube.spectral_axis,sp_av)
        plt.yticks(np.arange(min(sp_av), max(sp_av), (max(sp_av)-min(sp_av))/4))
        
        if loop_count < (len(files) - 1):
            plt.setp(ax.get_xticklabels(), visible=False)

#         leg.append(str(file_name))
#         leg[-1] = leg[-1] + ', bw=' + str(bin_width) + ', bin#=' + str(wanted_bin)
        leg = [file_name]
        plt.legend(leg, bbox_to_anchor=(1.05, 1.1), prop={'size':8})
        
        loop_count += 1
        
    print 'leg =', leg
    print 'curves_names =', curves_names
    

#     plt.legend(leg, loc = 5, prop={'size':8})
#     ax.legend(bbox_to_anchor=(1.05, 0), loc='lower left', borderaxespad=0.)
#     fig.legend(curves_names, leg, bbox_to_anchor=(0.5, -0.15))
#     plt.savefig("bin_width=500wanted_bin=0.png")    
#     figtext(.0,.0,'From top to botom:\n' + str(files) , fontsize=8)

    fig.suptitle("average spectrum of pixels in brightnes bin = %r for bin width = %r" %(wanted_bin, bin_width), fontsize=12)
    plt.show()

   

    
def binning(loop_count, bin_width=500, wanted_bin=0):
    """A function creating brightness bins of pixels, and eventualy a map, in the given spectral cube"""
    # 1. create the array of bins - bin_arr
    # maybe use for that the hint below!
#     global bin_width, bin_arr, bin_arr2, cube, baddata, Tmax, legend_content
    
#     print "What file do you want to work with? (e.g. NGC1333_NH3_11_DR1.fits)"
#     file_name = raw_input()
    file_name = files[loop_count]
    cube = SpectralCube.read(file_name)
    cube = cube.with_spectral_unit(u.km/u.s,velocity_convention='radio')
    Tmax = cube.apply_numpy_function(np.nanmax,axis=0) # array of the maximum values in the spectra of each pixel
    baddata = nd.morphology.binary_dilation(np.isnan(Tmax),np.ones((25,25)))
    Tmax[baddata]=0.0
    Tmax[np.isfinite(Tmax)]
    
#     print "How many pixels would you like as a bin width in brightness?"
#     bin_width = input()
    bin_arr = np.sort(Tmax[np.isfinite(Tmax)])
    bin_arr2 = bin_arr[:: - bin_width] # this creates an array of the bin margins, in which every bin has a width of "bin_width"
#     print "The margins of the bins are at:", bin_arr2

#     2. use the array of bins for labeling all the pixels in which of the bins they belong
#     and also doing a plot of the pixel labels
    
    np.digitize(Tmax,bin_arr2)
#     plt.title('plot1')
#     plt.figure(figsize=(5,10))
#     plt.imshow(np.digitize(Tmax,bin_arr2))
#     plt.clf()


#     print "What bin value would you like to average the spectrum from?"
#     wanted_bin = input()
    bins = np.digitize(Tmax,bin_arr2)
#     bins = np.digitize(Tmax,binning()[1])
    y, x = np.where(bins==wanted_bin)

    return Tmax, np.digitize(Tmax,bin_arr2), bin_arr2, cube, file_name, bin_width, wanted_bin, y, x

# 3. creating a function that shows you the label of the pixel you are looking for, when you give the coordinates
def show_pix_label():
    global bin_arr, bin_arr2, bin_width, cube, Tmax, baddata
    print "What pixel are you interested in, for gatting the bin-index of its brigthness?"
    print "Provide its coordinates - first x then y on the image: "
    x = input()
    y = input()

    print "The bin-index of the brightness of the pixel is:", np.digitize(Tmax,bin_arr2)[y][x]
    

main()
