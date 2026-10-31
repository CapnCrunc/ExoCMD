import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure
from matplotlib.collections import LineCollection
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import warnings
import astropy as ap
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive as nea
import glob
import pandas as pd
warnings.filterwarnings("ignore", category=RuntimeWarning)

def depth_to_absmag(pname, depth, type):
    if type == 'depth':
        m_s = str(nea.query_object(pname, table = 'pscomppars', select = 'sy_vmag')[0][0])
        m_s = float(m_s)
        app = -2.5*np.log10(depth/100) + m_s
    elif type == 'irradiance':
        app = -2.5*np.log10(depth/3.61e-11)
    if pname == 'solar spectra adjusted':
        dist = 4.848e-6
    else:
        dist = str(nea.query_object(pname, table = 'pscomppars', select = 'sy_dist')[0][0])
        dist = dist.replace(' pc', '')
        dist = float(dist)
    return app - 5*np.log10(dist/10)

def dta_err_prop(depth, derr):
    return ((-2.5)/(depth*np.log(10)))*derr

def calculations(Bands, Fname, interval, outfilename, type):
    output = {}
    output.setdefault('Planet Name', [])
    if type == 'irradiance':
        flux_name = 'Irradiance(W/m2/nm)'
    elif type == 'depth':
        flux_name = 'Depth(%)'
    for fname in glob.glob(Fname + '/*'):
        pname = fname.replace(Fname, '')
        slash = pname[0]
        pname = pname.replace(slash, '')
        output['Planet Name'].append(pname)
        pname = pname.replace('_', ' ')
        for i, b in enumerate(Bands):
            output.setdefault(str(b) + 'mag', [])
            output.setdefault(str(b) + 'mag_pos_err', [])
            output.setdefault(str(b) + 'mag_neg_err', [])
            df = pd.read_csv(fname)
            Wavelength = np.array(df['Wavelength(microns)'])
            Depth = np.array(df[flux_name])
            PDerr = np.array(df['Positive Depth Error'])
            NDerr = np.array(df['Negative Depth Error'])
            Wrange = []
            Drange = []
            PDErange = []
            NDErange = []
            for n, w in enumerate(Wavelength):
                if (b - interval[i]) <= w <= (b + interval[i]):
                    Wrange.append(w)
                    Drange.append(Depth[n])
                    PDErange.append(PDerr[n])
                    NDErange.append(NDerr[n])
                else:
                    continue
            Wrange = np.array(Wrange)
            Drange = np.array(Drange)
            PDErange = np.array(PDErange)
            NDErange = np.array(NDErange)
            if len(Drange) == 0:
                absmag = np.nan
                pamerr = np.nan
                namerr = np.nan
                output[str(b) + 'mag'].append(absmag)
                output[str(b) + 'mag_pos_err'].append(pamerr)
                output[str(b) + 'mag_neg_err'].append(namerr)
                continue
            elif len(Drange) == 1:
                absmag = depth_to_absmag(pname, Drange, type)
                pamerr = np.abs(dta_err_prop(Drange, PDErange))
                namerr = np.abs(dta_err_prop(Drange, NDErange))
                output[str(b) + 'mag'].append(float(absmag))
                output[str(b) + 'mag_pos_err'].append(float(pamerr))
                output[str(b) + 'mag_neg_err'].append(float(namerr))
                continue
            else:
                IntDepth = np.trapz(y = Drange, x = Wrange)
                IPDerr = np.sqrt(np.sum(PDErange**2)/(len(PDErange)))
                INDerr = np.sqrt(np.sum(NDErange**2)/(len(NDErange)))
                absmag = depth_to_absmag(pname, IntDepth, type)
                pamerr = np.abs(dta_err_prop(IntDepth, IPDerr))
                namerr = np.abs(dta_err_prop(IntDepth, INDerr))
                output[str(b) + 'mag'].append(absmag)
                output[str(b) + 'mag_pos_err'].append(pamerr)
                output[str(b) + 'mag_neg_err'].append(namerr)
    df2 = pd.DataFrame(output)
    df2.to_csv('depth to mag output/' + outfilename)