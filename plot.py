
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sn

from glob import glob
from os.path import join as pjoin

mpl.rcParams.update({
    "axes.labelsize": 20,
    "axes.titlesize": 18,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
})


def plotResult(result_path,
               save_dir="/Users/mt/workspace/Papers/gp-stat/Figures/Schkade2010/"):

    res = pd.read_csv(result_path, index_col=0)
    piv = res.pivot("postDelibSD", "preDelibSD", "successRate")

    plt.figure()
    ax = sn.heatmap(piv, square=True,
                    cbar_kws={'label': 'Success matching observed'})
    ax.invert_yaxis()

    basename = result_path[5:-4]
    plt.title(basename.replace('_', '\n').replace('=', ' = '))

    if save_dir:
        plt.savefig(pjoin(save_dir, basename + ".pdf"))


def makeSchkadeFigures():

    schkade_data_files = (
        glob(pjoin("data", "Boulder*.csv")) +
        glob(pjoin("data", "COSprings*.csv"))
    )
    # print(schkade_data_files)

    for f in schkade_data_files:
        plotResult(f)


def makeMoscoviciFigures():

    moscovici_data_files = glob(pjoin("data", "Moscovici*.csv"))

    for f in moscovici_data_files:
        plotResult(f,
                   save_dir='/Users/mt/workspace/Papers/gp-stat/Figures/Moscovici1969')
