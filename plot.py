
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sn

from glob import glob
from os.path import join as pjoin

mpl.rcParams.update({
    "axes.labelsize": 14,
    "axes.titlesize": 15,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
})


def plotResult(result_path,
               save_dir="/Users/mt/workspace/Papers/gp-stat/Figures"):

    res = pd.read_csv(result_path, index_col=0)
    piv = res.pivot("postDelibSD", "preDelibSD", "successRate")

    plt.figure()
    ax = sn.heatmap(piv, square=True,
                    cbar_kws={'label': 'Success matching observed'})
    ax.invert_yaxis()

    basename = result_path[5:-4]
    plt.title(basename.replace('_', ', '))

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
