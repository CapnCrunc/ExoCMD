from astroquery.simbad import Simbad
from astropy import coordinates as coords
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive as nea
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import figure
from matplotlib.collections import LineCollection
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import warnings
import astropy as ap
import glob
import pandas as pd
from difflib import SequenceMatcher
warnings.filterwarnings("ignore", category=RuntimeWarning)

def find_stellar_mag(outname):
    output = {}
    output.setdefault('Stellar Name', [])
    output.setdefault('J mag', [])
    output.setdefault('H mag', [])
    output.setdefault('K mag', [])

    for pname in glob.glob("organised emission spectra output/*"):
        pname = pname.replace('organised emission spectra output', '')
        backslash = pname[0]
        pname = pname.replace(backslash, '')
        pname = pname.replace('_', ' ')
        sname = nea.query_object(pname, table = "pscomppars", select = 'hostname')[0][0]
        dist = find_distance(pname)[0]
        result = Simbad.query_object(sname)[0]
        output['Stellar Name'].append(result['matched_id'])
        output['J mag'].append(app_to_abs(float(result['J']), dist))
        output['H mag'].append(app_to_abs(float(result['H']), dist))
        output['K mag'].append(app_to_abs(float(result['K']), dist))

    df = pd.DataFrame(output)
    df.to_csv('depth to mag output/stellar mags/' + outname + '.csv')

def find_distance(pname):
    vals = nea.query_object(pname, table = 'pscomppars', select = 'sy_dist, sy_disterr1, sy_disterr2')[0]
    dist = str(vals[0])
    posderr = float(vals[1])
    negderr = float(vals[2])
    dist = dist.replace(' pc', '')
    dist = float(dist)
    return dist, posderr, negderr

def app_to_abs(app, dist):
    return app - 5*np.log10(dist/10)
    
def pdepth_to_absmag(pname, depth, type, m_s):
    if type == 'depth':
        app = -2.5*np.log10(depth/100) + m_s
    elif type == 'irradiance':
        app = -2.5*np.log10(depth/3.61e-11)
    if pname == 'solar spectra adjusted':
        dist = 4.848e-6
    else:
        dist = find_distance(pname)[0]
    return app_to_abs(app, dist)

def dta_err_prop(depth, derr, dist, Derr):
    a = (-2.5/(depth*np.log(10)))*derr
    b = (((1*a)*2) + (-5/(dist*np.log(10)))*Derr)
    return b

def calculations(Bands, Fname, interval, outfilename1, outfilename2, type):
    Simbad.server = "simbad.harvard.edu"
    Simbad.add_votable_fields('J', 'H', 'K')
    find_stellar_mag(outfilename2)
    output = {}
    output.setdefault('Planet Name', [])
    if type == 'irradiance':
        flux_name = 'Irradiance(W/m2/nm)'
    elif type == 'depth':
        flux_name = 'Depth(%)'
    sdf = pd.read_csv('depth to mag output/stellar mags/' + outfilename2 + '.csv')
    sname = np.array(sdf['Stellar Name'])
    sJ = np.array(sdf['J mag'])
    sH = np.array(sdf['H mag'])
    sK = np.array(sdf['K mag'])
    for j, fname in enumerate(glob.glob(Fname + '/*')):
        pname = fname.replace(Fname, '')
        slash = pname[0]
        pname = pname.replace(slash, '')
        output['Planet Name'].append(pname)
        pname = pname.replace('_', ' ')
        check = SequenceMatcher(None, sname[j], pname).ratio()
        d = find_distance(pname)
        dist = d[0]
        pderr = d[1]
        nderr = np.abs(d[2])
        if check > 0.7:
            mJ = sJ[j]
            mH = sH[j]
            mK = sK[j]
            sbands = [mJ, mH, mK]
        else:
            sbands = [np.nan, np.nan, np.nan]
            print('planets and stars not in correct order')
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
                absmag = pdepth_to_absmag(pname, Drange, type, sbands[i])
                pamerr = np.abs(dta_err_prop(Drange, PDErange, dist, pderr))
                namerr = np.abs(dta_err_prop(Drange, NDErange, dist, nderr))
                output[str(b) + 'mag'].append(float(absmag))
                output[str(b) + 'mag_pos_err'].append(float(pamerr))
                output[str(b) + 'mag_neg_err'].append(float(namerr))
                continue
            else:
                IntDepth = np.trapz(y = Drange, x = Wrange)
                IPDerr = np.sqrt(np.sum(PDErange**2)/(len(PDErange)))
                INDerr = np.sqrt(np.sum(NDErange**2)/(len(NDErange)))
                absmag = pdepth_to_absmag(pname, IntDepth, type, sbands[i])
                pamerr = np.abs(dta_err_prop(IntDepth, IPDerr, dist, pderr))
                namerr = np.abs(dta_err_prop(IntDepth, INDerr, dist, nderr))
                output[str(b) + 'mag'].append(absmag)
                output[str(b) + 'mag_pos_err'].append(pamerr)
                output[str(b) + 'mag_neg_err'].append(namerr)
    df2 = pd.DataFrame(output)
    df2.to_csv('depth to mag output/planet mags/' + outfilename1)