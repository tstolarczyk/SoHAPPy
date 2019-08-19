# SET GAMMAPY_EXTRA to a decent value before starting
import matplotlib.pyplot as plt
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'

import os
os.environ['GAMMAPY_EXTRA']='D:\CTA\Analyse\gammapy-extra'
print(os.getenv('GAMMAPY_EXTRA'))
print(os.listdir(os.getenv('GAMMAPY_EXTRA')))

import numpy

# Gammapy
from gammapy.spectrum.models import Absorption
from gammapy.scripts import CTAPerf

# Utilties
from config_handler import ConfigHandler
from grb_utils import GammaRayBurst
from plot_utils import PlotGammaRayBurst


#filepath = '/Users/julien/Documents/WorkingDir/Tools/python/grb_paper/data/long_grb_test/LGRB1/'
filepath = '../data/long_grb_test/LGRB3/'

absorption = Absorption.read_builtin('dominguez')

# This is the native IRF in fits format
#south_irf_file = 'D:/CTA/Analyse/IRF-June2018/caldb/data/cta/prod3b/bcf/South_z20_N_100s/irf_file.fits'

# This is a fits from root file suited for point-like sources
south_irf_file = 'irf/CTA-Performance-South-20deg-100s_20170627.fits.gz'

cta_perf = CTAPerf.read(south_irf_file)

grb = GammaRayBurst.from_file(filepath=filepath, absorption=absorption)
grb.run_simulation(cta_perf)
print(grb)

PlotGammaRayBurst.plot_stats_detection(grb, savefig=True, outdir='./out/')
PlotGammaRayBurst.make_gif_from_models(grb, savefig=True, outdir='./out/')
plt.show(block=False)
