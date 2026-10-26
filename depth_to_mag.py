import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure
from matplotlib.collections import LineCollection
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import warnings
import astropy as ap
from astroquery import open_exoplanet_catalogue as oec
from astroquery.open_exoplanet_catalogue import findvalue
import glob
import pandas as pd
warnings.filterwarnings("ignore", category=RuntimeWarning)

def depth_to_appmag(depth):
    return -2.5*np.log10(depth/100)

def dta_err_prop(depth, derr):
    return ((-2.5)/(depth*np.log(10)))*derr

def calculations(Bands, Fname, interval, outfilename):
    output = {}
    output.setdefault('Planet Name', [])
    for fname in glob.glob(Fname + '/*'):
        pname = fname.replace(Fname, '')
        slash = pname[0]
        pname = pname.replace(slash, '')
        output['Planet Name'].append(pname)
        for i, b in enumerate(Bands):
            output.setdefault(str(b) + 'mag', [])
            output.setdefault(str(b) + 'mag_pos_err', [])
            output.setdefault(str(b) + 'mag_neg_err', [])
            df = pd.read_csv(fname)
            Wavelength = np.array(df['Wavelength(microns)'])
            Depth = np.array(df['Depth(%)'])
            PDerr = np.array(df['Positive Depth Error'])
            NDerr = np.array(df['Negative Depth Error'])
            Wrange = []
            Drange = []
            PDErange = []
            NDErange = []
            for n, w in enumerate(Wavelength):
                if (b - interval[i]) < w < (b + interval[i]):
                    Wrange.append(w)
                    Drange.append(Depth[n])
                    PDErange.append(PDerr[n])
                    NDErange.append(NDerr[n])
            Wrange = np.array(Wrange)
            Drange = np.array(Drange)
            PDErange = np.array(PDErange)
            NDErange = np.array(NDErange)
            if len(Drange) == 0:
                appmag = np.nan
                pamerr = np.nan
                namerr = np.nan
                output[str(b) + 'mag'].append(appmag)
                output[str(b) + 'mag_pos_err'].append(pamerr)
                output[str(b) + 'mag_neg_err'].append(namerr)
                continue
            elif len(Drange) == 1:
                appmag = depth_to_appmag(Drange)
                pamerr = np.abs(dta_err_prop(Drange, PDErange))
                namerr = np.abs(dta_err_prop(Drange, NDErange))
                output[str(b) + 'mag'].append(float(appmag))
                output[str(b) + 'mag_pos_err'].append(float(pamerr))
                output[str(b) + 'mag_neg_err'].append(float(namerr))
                continue
            else:
                IntDepth = np.trapz(y = Drange, x = Wrange)
                IPDerr = np.sqrt(np.sum(PDErange**2)/(len(PDErange)))
                INDerr = np.sqrt(np.sum(NDErange**2)/(len(NDErange)))
                appmag = depth_to_appmag(IntDepth)
                pamerr = np.abs(dta_err_prop(IntDepth, IPDerr))
                namerr = np.abs(dta_err_prop(IntDepth, INDerr))
                output[str(b) + 'mag'].append(appmag)
                output[str(b) + 'mag_pos_err'].append(pamerr)
                output[str(b) + 'mag_neg_err'].append(namerr)
    df2 = pd.DataFrame(output)
    df2.to_csv('depth to mag output/' + outfilename)