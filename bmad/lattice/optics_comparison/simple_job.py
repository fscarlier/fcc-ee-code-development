import os, sys
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt

from fcc_plots import fcc_axes

df_tap = pd.read_csv('./h/data_tapered.csv')
df_ref = pd.read_csv('./h/data_reference.csv')
df_pyat_ref = pd.read_csv('./h/data_reference_pyat.csv')
df_pyat_tap = pd.read_csv('./h/data_tapered_pyat.csv')

fig1, ax1 = fcc_axes()
ax1.plot(df_tap['S'], (df_tap['BETX_BMAD']-df_ref['BETX_BMAD'])/df_ref['BETX_BMAD'], label='Bmad')
ax1.plot(df_tap['S'], (df_tap['BETX_MADX']-df_ref['BETX_MADX'])/df_ref['BETX_MADX'], label='MAD-X')
ax1.plot(df_pyat_tap['S'], (df_pyat_tap['BETX_PYAT']-df_pyat_ref['BETX_PYAT'])/df_pyat_ref['BETX_PYAT'], label='PyAT')
ax1.set_ylabel(r'$\frac{\Delta \beta_x}{\beta_x^{Ref}}$', fontsize=16)
ax1.legend(loc='upper right', fontsize=12, framealpha=1)

fig2, ax2 = fcc_axes()
ax2.plot(df_tap['S'], df_tap['DX_BMAD']-df_ref['DX_BMAD'], label='Bmad')
ax2.plot(df_tap['S'], df_tap['DX_MADX']-df_ref['DX_MADX'], label='MAD-X')
ax2.plot(df_pyat_tap['S'], df_pyat_tap['DX_PYAT']-df_pyat_ref['DX_PYAT'], label='PyAT')
ax2.set_ylabel(r'$\Delta \eta_x$', fontsize=16)
ax2.legend(loc='upper right', fontsize=12, framealpha=1)

plt.show()
