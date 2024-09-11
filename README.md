# past-coding-unorganised

Welcome! This is a fairly unorganised repository of code that I've made in the past for things (at least, it's all the extensive code files that I could find.)

The files muon_fitting_code.py and read_dat_file.py are from a class project in which we used a scintillator + photomultiplier setup to estimate the lifetime of the muon. We did this by having our instruments record only signals occuring within a couple of microseconds of each other, as this situation would, at least some of the time, correspond to a muon entering the scintillator (producing a signal, coming to rest there, and then decaying (producing a second signal). The signals given by the PMT were collected with multiple simultaneous methods that produced different types of data files.

The file xrfhisto.py is from work I did on on an X-ray fluorescence spectrometer in the Summer of 2023. At the time, the data readout software for the spectrometer wasn't working correctly. To temporarily deal with this, I wrote Python code that read the .hyperc files created by the spectrometer and converted them into histograms plottable in matplotlib.
