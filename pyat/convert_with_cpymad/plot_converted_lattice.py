import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
from at.physics import linopt, linopt_rad
from fcc_plots import fcc_axes
from convert_madx_to_pyat import convert_madx_lattice_pyat
# MAD-X imports
from cpymad.madx import Madx

# PyAT imports
import at


def import_at_lattice(lattice_path):
    key = 'ring'
    ring = at.load.matfile.load_mat(lattice_path,key=key)
    ring = at.Lattice(ring)
    return ring


def plot_optics(spos, l, sposc, lc):
    fig1, ax1 = fcc_axes()
    ax1.plot(spos, sposc-spos)

    fig2, ax2 = fcc_axes()
    ax2.plot(spos, (lc.beta[:,0] - l.beta[:,0])/l.beta[:,0])
    ax2.set_ylabel('beta x (conv - ref)/ref')

    fig3, ax3 = fcc_axes()
    ax3.plot(spos, (lc.beta[:,1] - l.beta[:,1])/l.beta[:,1])
    ax3.set_ylabel('beta y (conv - ref)/ref')

    fig4, ax4 = fcc_axes()
    ax4.plot(spos,lc.dispersion[:,0] - l.dispersion[:,0])
    ax4.set_ylabel('dispersion x (conv - ref)')

    fig5, ax5 = fcc_axes()
    ax5.plot(spos, lc.closed_orbit[:,0] - l.closed_orbit[:,0])
    ax5.set_ylabel('orbit x (conv - ref)')

    fig6, ax6 = fcc_axes()
    ax6.plot(spos, lc.closed_orbit[:,2] - l.closed_orbit[:,2])
    ax6.set_ylabel('orbit y (conv - ref)')

    fig7, ax7 = fcc_axes()
    ax7.plot(spos, lc.closed_orbit[:,4] - l.closed_orbit[:,4])
    ax7.set_ylabel('orbit z (conv - ref)')

    fig8, ax8 = fcc_axes()
    ax8.plot(spos, lc.closed_orbit[:,5] - l.closed_orbit[:,5])
    ax8.set_ylabel('orbit pz (conv - ref)')
    plt.show()


def get_indexes(lat):
    dip = at.get_refpts(lat, at.lattice.elements.Dipole)
    quad = at.get_refpts(lat, at.lattice.elements.Quadrupole)
    sext = at.get_refpts(lat, at.lattice.elements.Sextupole)
    idx = np.concatenate([dip, quad, sext])
    idx.sort()
    return idx


def do_main():
    energy = 120e9 
    madx_lattice_path = 'lattice.seq'
    ring = convert_madx_lattice_pyat(madx_lattice_path, energy)
    return ring


if __name__ == '__main__':
    ring_conv = do_main()
    ring = import_at_lattice('fcch_norad.mat')
    ring_conv.radiation_off()
    ring.radiation_off()

    idx = get_indexes(ring)
    idx_conv = get_indexes(ring_conv)

    get_chrom = True
    xy_step = 1.0e-9
    dp_step = 1.0e-6

    l0,q,qp,l = linopt(ring,dp=0,refpts=idx,get_chrom=get_chrom,
                       coupled=False, XYStep=xy_step, DPStep=dp_step)
    spos = ring.get_s_pos(idx)

    l0c,qc,qpc,lc = linopt(ring_conv,dp=0,refpts=idx_conv,get_chrom=get_chrom,
                       coupled=False, XYStep=xy_step, DPStep=dp_step)
    sposc = ring_conv.get_s_pos(idx_conv)
    plot_optics(spos, l, sposc, lc)
