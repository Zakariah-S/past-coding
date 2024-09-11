'''The below code pertains to the output of an X-ray Fluorescence Spectrometer (XRF) that I worked on in the Summer of 2023.
The code contains functions that load relevant data from an input file, calibrate data based on the analysis of a known material,
and plot the spectra of an observation accordingly. 
There is also code that writes out the spectrum histograms as .hyperc files (a custom file format for our XRF), 
and that compares two spectra to locate similarities and differences.'''

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

#Get plottable lists of bins, counts from .txt files. These bins and counts can then be plotted with matplotlib.

def load_data(file):
    #Put in full file path as argument, get out two lists of bins and counts. The entry at bins[i] corresponds to the entry at counts[i]
    bins = []
    counts = []
    with open(file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.split('\t')[1]:
                [bin, count] = line.split('\t')
                bins.append(int(bin))
                counts.append(int(count))
    return bins, counts

''' Uses a chosen file to calibrate data by matching known peaks to the MCA bins at which they appear.
    At least two peaks must be inserted manually. 
    The function uses these to perform a linear fit and give an equation for energy in terms of the bin number.
    The parameters of this fit are returned and used in plot_histo.'''

def calibration(cal_file = '', gauss = False, do_write = False, write_txt='xrfpeaks.txt'):

    ''' gauss:      When set to false, use manually-inserted bin number as the putative peak location.
                    When set to true, run Gaussian fit using the manually-inserted bin number to get a computed estimate of the peak location.
        do_write:   When set to true, write out a .txt file with peak information. 
                    When set to false, do nothing.'''

    ''' Initialise dictionary of peaks that we already know.
        To start an entry: 'Name of peak' : [Bin no. of peak, Counts at bin no., Peak Energy (keV)]
        If gauss = True, program later fills this in with a Gaussian estimation as follows:
        'Name of peak' : [Inputted peak bin, Inputted counts, Peak Energy (keV), Estimated peak bin, Estimated max counts, Bin error]'''
    
    maxima = {
        'CuKa' : [2139, 225, 8.040],
        'ZnKa' : [2292, 123, 8.637]
    }

    #Line function which we'll use to derive a conversion from bin number to energy
    def line(x, m, b):
        return m*x + b

    if gauss == True: 
    
        #Insert calibration file, load bins vs. counts data from it
        cbins, ccounts = load_data(cal_file)

        #Define gauss function to be fit to each peak
        def gauss(x, M, A, S):
            y = A * np.exp(-((x-M)**2)/(2*S**2))
            return y

        #Fit Gaussians around each maximum, save params to dictionary
        for max in maxima:
            xdata = list(range(maxima[max][0]-40, maxima[max][0]+40))
            ydata = ccounts[maxima[max][0]-40:maxima[max][0]+40]
            par, cov = curve_fit(gauss, xdata, ydata, [maxima[max][0], maxima[max][1], 1.])
            print(par)
            maxima[max].append(par[0])
            maxima[max].append(par[1])
            maxima[max].append(par[2])
    
        #Write to file, if needed
        if do_write == True:
            with open(write_txt, 'w') as f:
                f.write('Inserted peaks:\n\n')
                for max in maxima:
                    f.write('{}: Energy {} keV at bin {}, Peak counts = {}\n'.format(max, maxima[max][2], maxima[max][0], maxima[max][1]))
                f.write('\nUsing Gaussians:\n\n')
                f.write('')
                for max in maxima:
                    f.write('{}: Energy {} keV at bin {} +/- {}, Peak counts = {}'.format(max, maxima[max][2], maxima[max][3], maxima[max][5], maxima[max][4]))

        bin_nos = [maxima[max][3] for max in maxima]
        bin_error = [maxima[max][5] for max in maxima]
        energies = [maxima[max][2] for max in maxima]

        par, cov = curve_fit(line, energies, bin_nos, bin_error)

        m, b = 1/par[0], par[1]/par[0]

        print(m, b)
        return m, b

    if gauss == False:
        #Then we just leave the working dictionary untouched lol
        if do_write == True:
            with open(write_txt, 'w') as f:
                f.write('Inserted peaks:\n\n')
                for max in maxima:
                    f.write('{}: Energy {} keV at bin {}, Peak counts = {}\n'.format(max, maxima[max][2], maxima[max][0], maxima[max][1]))

        bin_nos = [maxima[max][0] for max in maxima]
        energies = [maxima[max][2] for max in maxima]
        print(bin_nos)
        print(energies)
        par, cov = curve_fit(line, bin_nos, energies)
        m, b = par[0], par[1]

        print(m, b)
        return m, b

''' Plots spectra; can take either .txt or .hyperc files as input. 
    For .txt files, use the 'Full' mode. For .hyperc files, use the 'Separate' mode.'''

def plot_histo(file, mode = 'Full', ROI = [], cal_file = '/Users/zaksaeed/Downloads/xrfdata/20230615testcyl_detector1'):
    
    ''' file: File (.txt or .hyperc file with full path) containing histogram that you want to plot.
        mode: Selecting 'Full' plots the entire spectrum from the entered .txt file. 
                Only currently works for .txt files, so if you only have a .hyperc you'll need ot convert it!
              Selecting 'Separate' plots the spectrum over a selected rectangle of a raster image. 
                The rectangle corners must be defined in the arguments of separate_point_spectra().
                This mode only works with .hyperc files.
        ROI: Region of interest. Sets limits of x-axis.
        cal_file: File that you're using to calibrate (match MCA bin to energy). 
              At the time of my writing this, this should be a file with the spectrum of brass.
              If there is no calibration file, then the x-axis will simply be the MCA bin number.'''
    
    print(mode)
    #Call different functions, based on the mode.
    if mode == 'Full':
        bins, counts = load_data(file)
    elif mode == 'Separate':
        bins, counts = separate_spectra(file)
    else:
        raise Exception("Please enter a valid mode for plotting-- either 'Full' or 'Separate'.")

    #Convert from MCA bin to energy (keV) if a calibration file is provided. Plot accordingly.
    if cal_file:
        m, b = calibration(cal_file)
        energies = []
        for bin in bins:
            energies.append(m*bin + b)
        plt.plot(energies, counts, '-', linewidth = 0.7)
        plt.xlabel('Energy (keV)')
    else: 
        plt.plot(bins, counts, '-', linewidth = 0.7)
        plt.xlabel('Bin')
    plt.ylabel('Counts')
    plt.title("X-ray Spectrum")

    #If ROI is defined, set x limits accordingly.
    if ROI != []:
            plt.xlim((max(ROI[0], 0), ROI[1]))
    plt.show()

    print(sum(counts))

#Returns a bins and counts list for a rectangular selection over a raster image. Only compatible with .hyperc files for raster scans.

def separate_spectra(file, rectangle = []):

    ''' file: Must be a .hyperc file.
        rectangle: Enter the ordered-pair coordinates of two opposite corners of the rectangle whose spectrum you want. 
            Be sure that your chosen rectangle doesn't exceed the limits of the image. 
            An example of what might be entered here is: [(1, 1), (3, 3)]. This writes out a histogram for the rectangle defined by these corners.'''
    
    #Open file, extract lines with spectral data
    with open(file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if "<Analysis_Data>" in line:
                start = lines.index(line) + 1
            if "</Analysis_Data>" in line:
                end = lines.index(line)
            if "<Image_Info>" in line:
                index = lines.index(line)
                totwidth = int(lines[index+1].split("</Width>")[0].split("<Width>")[-1])
                totheight = int(lines[index+2].split("</Height>")[0].split("<Height>")[-1])
        pixels = lines[start:end]
        
    #Check that rectangle is well-defined and falls within image. If there isn't a rectangle, plot the whole spectrum.
    if rectangle != []:
        if rectangle[0][0] < 0 or rectangle[0][1] < 0 or rectangle[1][0] >= totwidth or rectangle[1][1] >= totheight:
            raise Exception("Please make sure your area of interest lies within the bounds of the image.")
        
        try:
            width = rectangle[1][0] - rectangle[0][0]
            height = rectangle[1][1] - rectangle[0][1]
        except:
            raise Exception("Make sure to enter numbers for the rectangle corners.")

    else:
        width = totwidth
        height = totheight

    #Setup bins and counts lists
    bins = list(range(16384))
    counts = [0] * 16384

    #Add data from each pixel, one by one
    for i in range(width):
        for j in range(max(height, 1)):
            for pixel in pixels:
                if pixel.split(",")[0:2] == [str(i), str(j)]:
                    data = pixel.split(",")[2:]
                    #print(data)
                    #print(len(data))
                    n = 0
                    for item in data:
                        #print("n = " + str(n))
                        #print(item)
                        if 'R' in item:
                            zeroes = int(item.split('R')[0].strip(), base=16)
                            n += zeroes
                        else: 
                            item = int(item.strip(), base=16)
                            counts[n] += item
                            n += 1
    
    '''with open('./xrfcyldet1.txt', 'w') as f:
        for i in range(len(bins)):
            f.write('{}\t{}\n'.format(bins[i], counts[i]))'''

    return bins, counts

#Input a .hyperc file from a raster scan and get a folder with one histogram .txt file for each pixel of the raster scan.

def write_pixel_spectra_files(file):

    '''file: Must be a .hyperc file.'''

    import os

    #Read file, get pixel spectra
    with open(file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if "<Analysis_Data>" in line:
                start = lines.index(line) + 1
            if "</Analysis_Data>" in line:
                end = lines.index(line)
            if "<Image_Info>" in line:
                index = lines.index(line)
                width = int(lines[index+1].split("</Width>")[0].split("<Width>")[-1])
                height = int(lines[index+2].split("</Height>")[0].split("<Height>")[-1])
        pixels = lines[start:end]

    #Define names for new directory and files, set up new directory (warn user if directory already exists)
    filename = file.split('/')[-1].split(".")[0]
    newpath = file.split('.')[0] + '/'
    print(file)
    print(filename)
    print(newpath)

    if width == 1 and height == 1:
        raise Exception("Given .hypercfile is for a point scan. Please use data from a raster scan.")
    #If file already processed by this program, give option to end the program here. 
    if os.path.isdir(newpath) == False:
        os.mkdir(newpath)
    else:
        answer = input("A folder with this file name has already been made. Proceeding WILL save over these data. Enter 'Y' to continue anyway, and enter anything else to abort: ")
        if answer != 'Y':
            raise Exception("Ended program.")

    #Get bins, counts list for each pixel, save to file with format {original file name}-{xcoord}{ycoord}.txt
    for pixel in pixels:
        temp = pixel.split(",")
        x = temp[0]
        y = temp[1]
        data = temp[2:]
        bins = list(range(16384))
        counts = [0] * 16384
        print(len(data))
        n = 0
        for item in data:
            #print("n = " + str(n))
            if 'R' in item:
                zeroes = int(item.split('R')[0].strip(), base=16)
                n += zeroes
            else: 
                item = int(item.strip(), base=16)
                counts[n] += item
                n += 1
        print(newpath + filename + '-' + x + y + '.txt')
        with open(newpath + filename + '-' + x + y + '.txt', 'w') as f:
            for i in range(len(bins)):
                f.write('{}\t{}\n'.format(bins[i], counts[i]))

#Compare two files and tell us whether there are any differences (as well as how many and where they are).

def are_spectra_same(file1, file2, verbose = True):

    ''' file1, file2: .txt files to be compared.
        verbose: If True, print feedback about differences. 
            (Default is true if this program is run alone, but is set to false in the below programme which can loop over this one a LOT.)'''

    #Open the two files, get the contents of each. 
    with open(file1, 'r') as f:
        with open(file2, 'r') as g:
            lines1 = f.readlines()
            lines2 = g.readlines()
            if len(lines1) != len(lines2):
                if verbose: print("Files have a different number of lines!")
                return 1
            #Count the number of lines that are different, record the line numbers that are different
            n = 0
            diffs = []
            for i in range(len(lines1)):
                if lines1[i] == lines2[i]:
                    continue
                else:
                    n += 1
                    diffs.append(i + 1)
    
    '''Print info about differences:
        If n = 0, tells us the files are the same.
        If n > 0, tells us that the files are different, as well as the lines at which they differ.'''
    if n > 1:
        if verbose:
            print(file1 + ' & ' + file2 + ':')
            statement = "Files are different! They have {} differences, at lines ".format(n)
            for j in range(len(diffs)-1):
                statement += '{}, '.format(diffs[j])
            statement += ' and {}.'.format(diffs[-1])
            print(statement)
    elif n == 1:
        if verbose: 
            print(file1 + ' & ' + file2 + ':')
            statement = "Files are different! They differ at line {}.".format(diffs[0])
            print(statement)
    elif n == 0:
        if verbose: 
            print(file1 + '\n' + file2)
            print("File contents are the same!")
    return n

#Check an entire directory full of files to see if any of them are the same.

def are_many_spectra_same(dir):

    #dir: Directory of files that you want to check.

    #Get list of files in chosen directory
    import os
    files = [f for f in os.listdir(dir) if not f.startswith('.')]
    print(files)
    print('\n')

    n = len([name for name in files if os.path.isfile(os.path.join(dir, name))])
    
    #Run are_spectra_same for all possible combinations of files. For each, record whether the files are the same.
    z = 0
    sames = []
    for i in range(n):
        for j in range(i+1, n):
            #print(os.path.join(dir, files[i]) + '\n' + os.path.join(dir, files[j]) + '\n')
            if are_spectra_same(os.path.join(dir, files[i]), os.path.join(dir, files[j]), verbose = False) == 0:
                z += 1
                sames.append((files[i], files[j]))

    #Print out any pairs of files that match, or let us know that all files are different. 
    if z == 0:
        print("Everything is different!")

    else:
        print("Pairs of matching files:\n")
        for pair in sames:
            print(pair[0])
            print(pair[1] + '\n')

#Recreate PyMCA plot from a .csv file generated by PyMCA.

def plot_from_pymcacsv(file, name, elims, save=False, show=True):
    ''' file: .csv file (produced by PyMCA) from which you wish to produce a plot.
        name: Title of plot (will appear at top).
        elims: (Minimum energy on plot, Maximum energy on plot)
        save: Boolean indicating whether you want to save the plot as a .png.
        show: Boolean indicating whether you want to show the plot.'''
    
    import pandas as pd
    
    data = pd.read_csv(file, index_col=1) #Reads in csv as dataframe

    #print(data)
    #data.info()

    #Rename stuff so we have better plot titles
    data.rename(columns={
        "counts": "Data",
        "fit": "Cumulative Fit",
    }, inplace=True)

    #Drop columns that we don't need, like the channel number and the (unused) pileup column
    data.drop(columns=["channel", "continuum", "pileup"], inplace=True)

    colors = [
        'black',
        'gray',
        'r',
        'b',
        'g',
        'm',
        'orange',
        'magenta',
        'blueviolet',
        'aqua',
        'gold',
    ]

    #Ensure that every column gets assigned a colour
    while len(colors) < len(data.columns):
        colors += colors

    colors = colors[0:len(data.columns)]

    #Plot data, give special formatting to the curves representing the data and total fit
    ax = data.plot(legend=True, title = name, lw=1.5, xlim=elims, figsize=(14, 5), 
              color=colors, grid=True, xlabel="Energy (keV)", ylabel="Counts")
    
    for line in ax.get_lines():
        if line.get_label() == "Data": 
            line.set_linewidth(0.5)
            line.set_zorder(-1)
        if line.get_label() == "Cumulative Fit":
            line.set_linewidth(2)
            line.set_zorder(-2)
            line.set_linestyle("dashed")

    if save==True: plt.savefig("/Users/zaksaeed/testfile.png", dpi=600)

    if show == True: plt.show()

if __name__ == '__main__':
    #calibration(gauss = False, do_write = False)
    #calibration(gauss = True, do_write = False)
    #plot_histo('/Users/zaksaeed/Downloads/20230615testcyl2_detector1', ROI = [6, 12])
    #plot_histo('/Users/zaksaeed/Downloads/xrf/20230615test_detector1', ROI = [0, 12])
    #separate_point_spectra('/Users/zaksaeed/Downloads/20230615testcyl_detector1.xhyperc')
    #plot_histo('/Users/zaksaeed/Downloads/xrfdata/20230615testcyl_detector1', mode = 'Full', ROI = [0, 12])
    #plot_histo('/Users/zaksaeed/Downloads/20230615testcyl_detector1.xhyperc', mode = 'Separate', ROI = [0, 12])
    #write_pixel_spectra_files('/Users/zaksaeed/Downloads/20230615testcyl_detector1.xhyperc')
    #are_spectra_same('/Users/zaksaeed/XRF/20230622_GeneralPointTest-20kV-20microA_detector0', '/Users/zaksaeed/XRF/20230622_GeneralPointTest-20kV-20microA_detector1')
    #are_many_spectra_same('/Users/zaksaeed/XRF/NewYellowGreen')
    plot_from_pymcacsv("/Users/zaksaeed/XRF/CometLike/CometLike_Yellow_10m_30kV_30uA_det0.csv", 
                       name = "Comet-like Sample, Yellow Side",
                       elims=(0, 9))
