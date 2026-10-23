
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

#figure(figsize = (20, 20))

    #Pname = ['GJ_341_b', 'GJ_357_b', 'GJ_486_b', 'GJ_1132_b', 'GJ_1214_b', 'GJ_3090_b', 'GJ_3470_b', 'HAT_P_12_b', 'HAT_P_14_b', 'HAT_P_26_b', 
         #'HAT_P_30_b', 'HD_15337_c', 'HD_189733_b', 'HD_209458_b', 'HIP_67522_b', 'K2_18_b', 'KELT_7_b', 'L_98_59_b', 'L_98_59_c', 'L_98_59_d', 
         #'L_168_9_b', 'LHS_475_b', 'LHS_1140_b', 'NGTS_2_b', 'TOI_199_b', 'TOI_260_b', 'TOI_270_b', 'TOI_270_d', 'TOI_421_b', 'TOI_732_c', 
         #'TOI_776_b', 'TOI_776_c', 'TOI_836_b', 'TOI_836_c', 'TOI_1130_b', 'TOI_1685_b', 'TrES_4_b', 'V1298_Tau_b', 'WASP_15_b', 'WASP_39_b', 
         #'WASP_80_b', 'WASP_94_A_b', 'WASP_107_b', 'WASP_121_b', 'WASP_127_b', 'WASP_166_b']

#plots = np.size(Pname)
#fig, axs = plt.subplots(plots)
#fig.suptitle('Planet Spectra')

#for ax in axs.flat:
    #ax.set(xlabel='Wavelength (microns)', ylabel='Transmission Depth (%)')

#for ax in axs.flat:
    #ax.label_outer()

def main_code(Pnames, Fname, Lskip):

    for m, pname in enumerate(Pnames):
        Wavelength = []
        Depth = []
        DepthPosErr = []
        DepthNegErr = []

        for fname in glob.glob(Fname + '/*'):
            check = pname in fname
            if check == True:
                data = np.genfromtxt(fname, skip_header = Lskip, unpack = True)
                wavelength = data[0]
                depth = data[2]
                depth_pos_err = data[3]
                depth_neg_err = np.abs(data[4])
                try:
                    for n, x in enumerate(wavelength):
                        Wavelength.append(x)
                        Depth.append(depth[n])
                        DepthPosErr.append(depth_pos_err[n])
                        DepthNegErr.append(depth_neg_err[n])
                except:
                    Wavelength.append(wavelength)
                    Depth.append(depth)
                    DepthPosErr.append(depth_pos_err)
                    DepthNegErr.append(depth_neg_err)

        Wavelength = np.array(Wavelength)
        Depth = np.array(Depth)
        DepthPosErr = np.array(DepthPosErr)
        DepthNegErr = np.array(DepthNegErr)

        WavelengthPosSort = Wavelength.argsort()
        WavelengthSorted = np.array(Wavelength[WavelengthPosSort])
        DepthSorted = np.array(Depth[WavelengthPosSort])
        DepthPosErrSorted = np.array(DepthPosErr[WavelengthPosSort])
        DepthNegErrSorted = np.array(DepthNegErr[WavelengthPosSort])

        SortedData = {'Wavelength(microns)': WavelengthSorted, 
                      'Transmission Depth(%)': DepthSorted, 
                      'Positive Depth Error': DepthPosErrSorted, 
                      'Negative Depth Error': DepthNegErrSorted}

        df = pd.DataFrame(SortedData)
        df.to_csv('spectra organised data/' + pname, index = False)

        MinWavelength = WavelengthSorted[0]
        MinWTransmission = DepthSorted[0]
        MinWTPErr = DepthPosErrSorted[0]
        MinWTNErr = DepthNegErrSorted[0]
        MaxWavelength = WavelengthSorted[-1]
        MaxWTransmission = DepthSorted[-1]
        MaxWTPErr = DepthPosErrSorted[-1]
        MaxWTNErr = DepthNegErrSorted[-1]

        plt.errorbar(WavelengthSorted, DepthSorted, yerr = np.array([DepthPosErrSorted, DepthNegErrSorted]), fmt = '-o', label = 'Spectrum', 
                     zorder = 1)
        plt.errorbar(MinWavelength, MinWTransmission, yerr = np.array([[MinWTPErr, MinWTNErr]]).T, fmt = 'o', 
                     label = 'Minimum Wavelength:' + str(MinWavelength), zorder = 2)
        plt.errorbar(MaxWavelength, MaxWTransmission, yerr = np.array([[MaxWTPErr, MaxWTNErr]]).T, fmt = 'o', 
                     label = 'Maximum Wavelength:' + str(MaxWavelength), zorder = 2)
        plt.title(pname + ' Spectra')
        plt.xlabel('Wavelength (microns)')
        plt.ylabel('Transmission Depth (%)')
        plt.legend()
        plt.savefig('spectra plots/' + pname)
        plt.show()