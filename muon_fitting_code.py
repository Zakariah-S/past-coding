import numpy as np
import matplotlib.pyplot as plt
from scipy.odr import *
from scipy.optimize import curve_fit

def gausscalmaxima(do_write = False):

    global bins
    global counts

    cbins = []
    ccounts = []
    cdict = {}
    peak_times = [1.0036, 2.003, 4.004, 6.005, 8.005, 10.0037, 12.002, 14.020, 16.013, 18.014, 19.022]
        
    with open('cali2.txt', 'r') as f:
        lines = f.readlines()
        for line in lines:
            cbins.append(int(line.split(': ')[0]))
            ccounts.append(int(line.split(': ')[1]))
            cdict[int(line.split(': ')[0])] = int(line.split(': ')[1])
            if int(line.split(': ')[0]) == 8000:
                break

    def gauss(x, A, M, S):
        y = A * np.exp(-((x-M)**2)/(2*S**2))
        return y
    
    cmaxima = []
    with open('calpeaks2.txt', 'r') as f:
        lines = f.readlines()
        for line in lines:
            try: cmaxima.append(int(line.split(' ')[2]))
            except: continue
    #print(cmaxima)
    n = 1
    peakdict = {}

    for max in cmaxima:
        xdata = list(range(max-11, max+12))
        ydata = ccounts[max-12:max+11]
        par, cov = curve_fit(gauss, xdata, ydata, [cdict[max], max, 0.5])
        #print("Peak {}: \nMax Bin = {}   Error = {}".format(n, par[1], par[2]))
        #print("Alt. error: {}".format(np.sqrt(cov[1][1])))
        peakdict[n] = [par[1], par[2]]
        #print(np.sqrt(cov[2][2]))
        n += 1
    
    if do_write == True:
        with open('calpeaks2.txt', 'a') as f:
            f.write('\nUsing Gaussians:\n\n')
            for i in range(len(peak_times)):
                f.write('{} us: {} {} {}\n'.format(peak_times[i], cmaxima[i], peakdict[i+1][0], peakdict[i+1][1]))

    return peakdict, peak_times

def avgcalmaxima(do_write = False):
    counts = []
    #with open('calibration.txt', 'r') as f:
    with open('cali2.txt', 'r') as f:
        lines = f.readlines()
        for line in lines:
            counts.append(int(line.split(': ')[1]))
    maxima = {}
    avgmaxima = {}
    error = {}
    peakdict = {}
    peak_times = [1.0036, 2.003, 4.004, 6.005, 8.005, 10.0037, 12.002, 14.020, 16.013, 18.014, 19.022]

    n = 1
    for i in range(len(counts)):
        if counts[i] > 1000:
            if counts[i] > counts[i-1] and counts[i] > counts[i+1]:
                maxima[n] = i+1
                sum = 0
                sumsquares = 0
                totcounts = 0
                for k in range(i-11, i+12):
                    sum += counts[k]*(k+1)
                    totcounts += counts[k]
                    sumsquares += counts[k]*(k+1)**2
                avg = sum/totcounts
                avgmaxima[n] = avg
                stderr = (sumsquares/totcounts-avg**2)**(0.5)
                error[n] = stderr
                peakdict[n] = [avg, stderr]
                n += 1
        if i == 8000:
            break
    #print(maxima)

    #for i in range(len(maxima)):
        #print('Peak {} at {}. Error: {}'.format(i+1, avgmaxima[i+1], error[i+1]))

    if do_write == True:
        with open('calpeaks2.txt', 'w') as f:
            f.write('Using averages:\n\n')
            for i in range(len(maxima)):
                f.write("{} us: {} {} {}".format(peak_times[i], maxima[i+1], avgmaxima[i+1], error[i+1]))
                f.write('\n')

    return peakdict, peak_times

def t_to_bin(method):
    '''If method is set to 1, use maxima that were obtained through fitting Gaussians.
       If method is set to 0, use maxima obtained through averaging.'''
    if method == 1:
        peakdict, peak_times = gausscalmaxima(False)
    elif method == 0:
        peakdict, peak_times = avgcalmaxima(False)

    def line(x, A, B):
        return A*x + B
    
    tdata = peak_times
    bindata = []
    error = []

    for i in range(len(peak_times)):
        bindata.append(peakdict[i+1][0])
        error.append(peakdict[i+1][1])
    
    '''for i in range(len(peak_times)):
        print('Time {}: Bin {}, Error {}'.format(peak_times[i], bindata[i], error[i]))'''
    par, cov = curve_fit(line, tdata, bindata, None, error, True)

    binfit = []

    for i in range(len(peak_times)):
        binfit.append(line(tdata[i], par[0], par[1]))

    print(par)
    print(cov)

    chisq = 0
    for i in range(len(tdata)):
        chisq += ((binfit[i]-bindata[i])**2)/binfit[i]

    print('Chi-squared = {}'.format(chisq))
    print("Reduced Chi-squared = " + str(chisq/(len(tdata)-2)))

    t = np.array(tdata)

    plt.errorbar(tdata, bindata, error, None, '.', label = 'Data', ecolor = 'black', elinewidth = 0.4, capsize = 0.4)
    plt.plot(tdata, binfit, '-', label = 'Best Fit: \n({0:.3f} \u00B1 {1:.3f})t + ({2:.3f} \u00B1 {3:.3f})'
             .format(par[0], np.sqrt(cov[0][0]), par[1], np.sqrt(cov[1][1])))
    plt.xlabel('Time elapsed (\u03BCs)')
    plt.ylabel('Bin Number')
    plt.legend(loc = 'upper left')
    plt.title('Bin Number vs. Time Interval')
    plt.show()

    return par[0], par[1], cov

def load_data(file, A, B, err, write=False):

    bins = []
    counts = []

    dA = err[0][0]
    dB = err[1][1]

    with open(file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            bins.append(int(line.split(': ')[0]))
            counts.append(int(line.split(': ')[1]))
            if int(line.split(': ')[0]) == 8000:
                break

    '''for i in range(len(bins)):
        print(str(bins[i]) + ', ' + str(counts[i]))'''

    tdata = []
    terr = []

    for i in range(len(bins)):

        #Remember: Bin_number = A * (time) + B. Thus, time = (1/A)(Bin_number - B)
        t = (bins[i] - B)/A + 0.032
        err = np.sqrt((((B-bins[i])*dA)/A**2)**2 + (dB/A)**2)
        tdata.append(t)
        terr.append(err)

    #print(terr)

    '''for i in range(len(bins)):
        print('Bin {}/{} us: {}'.format(bins[i], tdata[i], counts[i]))'''

    if write == True:
        with open('timedatafinal2.txt', 'w') as f:
            for i in range(len(bins)):
                f.write('Bin {}/{} us: {}'.format(bins[i], tdata[i], counts[i]))
                f.write('\n')

    '''for i in range(len(tdata)):
        print("{}, {}".format(tdata[i], terr[i]))'''

    return bins, tdata, counts, terr

def rebin(binno, tdata, counts, terr, avg = False):

    #binno = number of original bins to be placed into each new one.

    newt = []
    newc = []
    newterr = []

    n = 1
    count = 0
    
    if avg == False:
        for i in range(len(tdata)):
            count += counts[i]
            if n == 1:
                t = tdata[i]
                err = terr[i]
            if n == binno:
                newc.append(count)
                newt.append(t)
                newterr.append(err)
                n = 0
                count = 0
            n += 1
    
    elif avg == True:
        sumtimes = 0
        for i in range(len(tdata)):
            count += counts[i]
            sumtimes += tdata[i]
            if n == binno:
                newc.append(count)
                newt.append(sumtimes/binno)
                n = 0
                count = 0
                sumtimes = 0
            n += 1
    
    binwidth = newt[1]-newt[0]
    print(binwidth)

    return newt, newc, newterr

def expo(x, A, B, T):
    y = A*(np.exp(-x/T)) + B
    return y

def expo2(x, A, B, C, D, T):
    return A*np.exp(-x/T) + B + C*np.exp(-x/D)

def fitmca(file, binno):

    A, B, err = t_to_bin(1)
    bins, tdata, counts, terr = load_data(file, A, B, err)
    #newbins, newcounts, newt = rebin(binno, bins, tdata, counts)
    totcounts = sum(counts)
    tdata, counts, terr = rebin(binno, tdata, counts, terr, False)
    counterr = [np.sqrt(count) for count in counts]

    np.asarray(tdata)
    np.asarray(counts)

    '''for i in range(len(tdata)):
        print("Bin {}/{} us: {}".format(i+1, tdata[i], counts[i]))'''

    #cut off data points from the beginning that might fall below the curve

    maxbinindex = 0
    for i in range(len(tdata)):
        if counts[i] == max(counts):
            maxbinindex = i
            break
    
    print(maxbinindex)
    plt.bar(tdata, counts, width = 20/(8000/binno), label = 'Histogram Data')
    #plt.plot(tdata, counts, '.')

    #par, cov = curve_fit(expo, tdata[maxbinindex:], counts[maxbinindex:])
    par, cov = curve_fit(expo, tdata[maxbinindex:], counts[maxbinindex:])
    '''parA = par[0]
    parB = par[1]
    parC = par[2]
    parD = par[3]
    parT = par[4]'''
    parA = par[0]
    parB = par[1]
    parT = par[2]

    print(par)
    print(cov)

    fit_y = [expo(tdata[i], parA, parB, parT) for i in range(len(tdata))]
    '''for i in range(len(tdata)):
        fit_y.append(expo(tdata[i], parA, parB, parT))'''

    opt_fit = [expo(tdata[i], parA, parB, 2.117) for i in range(len(tdata))]
    '''for i in range(len(tdata)):
        opt_fit.append(expo(tdata[i], parA, parB, 2.196))'''
    chisq = 0
    for i in range(len(tdata[maxbinindex:])):
        chisq += ((fit_y[i+maxbinindex]-counts[i+maxbinindex])**2)/fit_y[i+maxbinindex]

    res = [fit_y[i]-counts[i] for i in range(len(fit_y))]

    print('Total counts: {}'.format(totcounts))
    print("Number of bins: " + str(len(tdata[maxbinindex:])))
    print('Chi-squared = {}'.format(chisq))
    print("Reduced Chi-squared = " + str(chisq/(len(tdata[maxbinindex:])-3)))
    print("Binwidth = {}".format(tdata[1]-tdata[0]))

    plt.title('Muon Decay Time')
    plt.xlabel('Decay time (\u03BCs)')
    plt.ylabel('Counts')
    plt.errorbar(tdata, counts, counterr, terr, fmt = 'none', ecolor = 'darkslategray', elinewidth = 0.8, capsize = 1)
    plt.plot(tdata, opt_fit, '-' 'g', label = 'Theoretical Fit: T = 2.117 \u03BCs', lw = 1.5)
    plt.plot(tdata, fit_y, '-' 'r', label = 'Best Fit: T = {0:.3f} \u00B1 {1:0.3f} \u03BCs'.format(parT, np.sqrt(cov[2][2])), lw = 2)
    plt.legend(loc = 'upper right')
    plt.show()

    plt.title('Residuals')
    plt.xlabel('Decay time (\u03BCs)')
    plt.ylabel('Residual')
    plt.errorbar(tdata[maxbinindex:], res[maxbinindex:], counterr[maxbinindex:], terr[maxbinindex:], fmt = 'none', ecolor = 'darkslategrey')
    plt.show() 

def fit_osc(file, binno):
    tdata = []
    counts = []
    with open(file, 'r') as f:
        lines = f.readlines()
        if lines[4].strip() != 'Time,Ampl':
            raise Exception("This doesn't seem to be an unaltered file from the oscilloscope.")
        for line in lines[5:]:
            counts.append(int(line.split(',')[1]))
            #time = line.split(',')[0]
            #val = float(time.split('e')[0])
            #exp = float(time.split('e')[1])
            t = (float(line.split(',')[0]) + 32e-9)
            tdata.append(t)
            if float(line.split(',')[0]) > 2e-5:
                break
    print("Original number of bins: " + str(len(tdata)))
    allcounts = sum(counts)
    newt = []
    newc = []
    count = 0
    n = 1

    for i in range(len(tdata)):
        count += counts[i]
        if n == 1:
            t = tdata[i]
        if n == binno:
            newc.append(count)
            newt.append(t)
            n = 0
            count = 0
        n += 1

    tdata = newt
    counts = newc
    newt = []
    newc = []

    for i in range(len(tdata)):
        if counts[i] != 0:
            newt.append(tdata[i])
            newc.append(counts[i])
    
    tdata = newt
    counts = newc
    tdata = [t*10**6 for t in tdata]
    cerror = [np.sqrt(count) for count in counts]

    startindex = 1
    endindex = 15
    par, cov = curve_fit(expo, tdata[startindex:], counts[startindex:], [200., 5., 2.1], cerror[startindex:])
    A, B, T = par[0], par[1], par[2]

    print(par)
    print(cov)

    fitdata = [expo(t, A, B, T) for t in tdata]
    optfit = [expo(tdata[i], A, B, 2.117) for i in range(len(tdata))]

    chisq = 0
    for i in range(len(tdata[startindex:endindex])):
    #for i in range(len(tdata[startindex:])):
        '''print(str(((fitdata[i+startindex]-counts[i+startindex])**2)/fitdata[i+startindex]) + ' ' + 
              '{} {}'.format(fitdata[i+startindex], counts[i+startindex]))'''
        chisq += ((fitdata[i+startindex]-counts[i+startindex])**2)/fitdata[i+startindex]

    print("Total counts: {}".format(allcounts))
    print("Number of bins: " + str(len(tdata[startindex:])))
    print('Chi-squared = {}'.format(chisq))
    print("Reduced Chi-squared = " + str(chisq/(len(tdata[startindex:endindex])-3)))
    #print("Reduced Chi-squared = " + str(chisq/(len(tdata[startindex:])-3)))
    print("Binwidth = {}".format(tdata[1]-tdata[0]))

    res = [fitdata[i]-counts[i] for i in range(len(fitdata))]
    #lnres = [np.log(abs(res[i])) for i in range(len(res))]

    plt.bar(tdata, counts, width = 0.75, label = 'Histogram Data')
    plt.errorbar(tdata, counts, cerror, None, fmt = 'none', ecolor = 'darkslategrey')
    plt.plot(tdata, optfit, '-' 'g', label = 'Theoretical Fit: T = 2.117 \u03BCs')
    plt.plot(tdata, fitdata, '-' 'r', label = 'Best Fit: T = {0:.3f} \u00B1 {1:0.3f} \u03BCs'.format(T, np.sqrt(cov[2][2])))
    plt.legend(loc = 'upper right')
    plt.xlabel('Decay time (\u03BCs)')
    plt.ylabel('Counts')
    plt.title('Muon Decay Time')
    plt.show()

    plt.title('Residuals')
    plt.xlabel('Decay time (\u03BCs)')
    plt.ylabel('Residual')
    plt.errorbar(tdata[startindex:], res[startindex:], cerror[startindex:], None, fmt = '.', ecolor = 'darkslategrey')
    plt.show() 

'''    plt.title('Residuals')
    plt.xlabel('Decay time (\u03BCs)')
    plt.ylabel('Residual')
    plt.errorbar(tdata[startindex:], lnres[startindex:], cerror[startindex:], None, fmt = '.', ecolor = 'darkslategrey')
    plt.show()
'''
if __name__ == '__main__':
    fitmca(file = 'data2.txt', binno = 100)
    #fitmca(file = 'datafinal.txt', binno = 100)
    #fit_osc('F1--TracedecayTimeGoodfinalmorebins.csv', binno=10)
    #fit_osc('F4--Trace--00001.csv', binno=30)
    #t_to_bin(1)
