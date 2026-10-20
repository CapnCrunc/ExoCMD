
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import warnings
import astropy as ap
from astroquery import open_exoplanet_catalogue as oec
from astroquery.open_exoplanet_catalogue import findvalue
warnings.filterwarnings("ignore", category=RuntimeWarning)

fname = 'JWST spectra\GJ_341_b_3.19698_5705_1.tbl'
data = np.genfromtxt(fname, skip_header = 13, unpack = True)
wavelength = np.array(data[0])
depth = np.array(data[2])
depth_pos_err = np.array(data[3])
depth_neg_err = np.array(np.abs(data[4]))

plt.errorbar(wavelength, depth, yerr = [depth_pos_err, depth_neg_err], fmt = 'o')