
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

cat = oec.get_catalogue()

for planet in cat.findall(".//planet"):
    print(findvalue(planet, "name"))