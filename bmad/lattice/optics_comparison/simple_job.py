import os, sys
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt

from fcc_plots import fcc_axes

df_tap = pd.read_csv('./h/data_tapered.csv')
df_ref = pd.read_csv('./h/data_reference.csv')

fig1, ax1 = fcc_axes()
ax1.plot(df_tap['S'], df_tap['BETX_BMAD'], label='Bmad')
ax1.plot(df_tap['S'], df_tap['BETX_MADX'], label='MAD-X')
ax1.set_ylabel(r"$\beta_x$ [m]", fontsize=16)
ax1.legend(loc='upper right', fontsize=12, framealpha=1)

plt.show()


