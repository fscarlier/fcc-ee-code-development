import at
import numpy as np
import matplotlib.pyplot as plt
from at.tracking import element_pass, lattice_pass
from at.physics import linopt, fast_ring, linopt_rad
from at.lattice import DConstant, PyElement, check_radiation
from at import load_mat, get_refpts, get_elements, get_value_refpts, set_value_refpts, elements, find_orbit6
from fcc_plots import fcc_axes
import time
import sys

@check_radiation(True)
def average_tapering(ring, multipoles=True, niter=1, **kwargs):
    """
    Delivers a scaled and averaged tapering for a predefined set of dipole families. Uses the result of perfect tapering to average over family set. Modifies the ring.
    PARAMETERS
        ring        individually tapered lattice
    KEYWORDS
        div_number      number of division within the loaded dipole family
        qadrupole=True  scale quadrupole fields
        sextupole=True  scale sextupole fields
        XYStep=1.0e-8   transverse step for numerical computation
        DPStep=1.0E-6   momentum deviation used for computation of orbit6
    """
    xy_step = kwargs.pop('XYStep', DConstant.XYStep)
    dp_step = kwargs.pop('DPStep', DConstant.DPStep)
    div_number = kwargs.pop('div_number')
    
    A = dipole_natural_family(ring)
    for Arcs in A:
        fam = family_segmentation(ring, div_number, Arcs)
        for dip in fam:
            for dipoles in dip:
                k_arc = get_value_refpts(ring, dipoles, 'PolynomB', index=0)
                k0_arc = get_value_refpts(ring, dipoles, 'BendingAngle')/get_value_refpts(ring, dipoles, 'Length')
                factor = np.average(k_arc/k0_arc, weights=None)
                set_value_refpts(ring, dipoles, 'PolynomB', factor*k0_arc, index=0)
                set_value_refpts(ring, dipoles, 'PolynomB', 0.0 , index=1)

    if multipoles:
        mults = get_refpts(ring, elements.Multipole)
        k0 = get_value_refpts(ring, dipoles, 'PolynomB', index=0)
        _, o6 = find_orbit6(ring, refpts=range(len(ring)+1),
                            XYStep=xy_step, DPStep=dp_step)
        dpps = (o6[mults, 4] + o6[mults+1, 4]) / 2
        for dpp, el in zip(dpps, ring[mults]):
            el.PolynomB *= 1+dpp
            el.PolynomA *= 1+dpp
        set_value_refpts(ring, dipoles, 'PolynomB', k0, index=0)
        
def family_segmentation(ring, div_number, dipole_indexes):
    """
    Selected dipole_indexes are divided into 4 groups as defined by the position of IPs and RFs, followed by the additional subdivision into 'div_number' of near-equal subgroups.
    """
    
    IP_indexes = get_refpts(ring, 'IP*')
    RF_indexes = get_refpts(ring, elements.RFCavity)
    
    #print('Number of dipoles per family:', len(dipole_indexes)/(4*div_number))
    #print('Number of dipoles:', len(dipole_indexes))

    ind1 = int(np.argwhere(dipole_indexes < RF_indexes[0])[-1])
    Arc1 = dipole_indexes[0:ind1]
    Arc1 = np.array_split(Arc1, div_number) # allows for non-equal segment splitting
    
    ind2 = int(np.argwhere(dipole_indexes > RF_indexes[19])[0])
    ind3 = int(np.argwhere(dipole_indexes < IP_indexes[1])[-1])
    Arc2 = dipole_indexes[ind2:ind3]
    Arc2 = np.array_split(Arc2, div_number)

    ind4 = int(np.argwhere(dipole_indexes > IP_indexes[2])[0])
    ind5 = int(np.argwhere(dipole_indexes < RF_indexes[20])[-1])
    Arc3 = dipole_indexes[ind4:ind5]
    Arc3 = np.array_split(Arc3, div_number)

    ind6 = int(np.argwhere(dipole_indexes > RF_indexes[-1])[0])
    ind7 = int(np.argwhere(dipole_indexes < IP_indexes[-1])[-1])
    Arc4 = dipole_indexes[ind6:ind7]
    Arc4 = np.array_split(Arc4, div_number)

    Arcs = [Arc1, Arc2, Arc3, Arc4]
    return Arcs
    
def dipole_natural_family(ring):
    """
    Extracts natural dipole families from the lattice and returns their indexes.
    
    PARAMETERS
        ring        individually tapered lattice
    """
    dipole_indexes = get_refpts(ring, elements.Dipole)
    
    # long and short - main arc dipoles
    dipole_B1 = get_refpts(ring, 'B1_*')
    dipole_B1S = get_refpts(ring, 'B1S_*')
    dipole_B1L = get_refpts(ring, 'B1L_*')
    dipole_main = get_refpts(ring, 'B1*')
    
    # dispersion suppressors
    dipole_BDS = get_refpts(ring, 'BDS*')
    # downstream of IP
    dipole_BC = get_refpts(ring, 'BC*')
    # dispersion suppressor for connecting arc
    dipole_BS = get_refpts(ring, 'BS*')
    # dipoles for connecting arc
    dipole_BG1 = get_refpts(ring, 'BG1_*')
    # dispersion suppressor for the upstream short arc (BLC is EMPTY)
    dipole_BL4 = get_refpts(ring, 'BL4*')
    dipole_BLC = get_refpts(ring, 'BLC*')
    # IP upstream dipoles
    dipole_BWL = get_refpts(ring, 'BWL_*')
    dipole_BCxL = get_refpts(ring, 'BC*L_*')
    
    return dipole_B1, dipole_B1S, dipole_B1L

if __name__ == "__main__":
    t0 = time.time()
    xy_step = 1.0e-9
    dp_step = 1.0e-6
    
    path = './Lattices/'
    filename = 'fcch_norad.mat'
    key = 'ring'
    latticef = path+filename

    # no tapering (ring), perfect tapering (ring_rad) and average tapering (ring_av) lattices
    ring = at.load_lattice(latticef, mat_key=key)
    ring = at.Lattice(ring)
    ring.radiation_off()
    l0, q, qp, ll = linopt(ring,dp=0,refpts=range(len(ring)),get_chrom = True, coupled=False, XYStep=xy_step, DPStep=dp_step)
    
    ring_rad = ring.deepcopy()
    ring_rad.radiation_on(quadrupole_pass='auto')
    ring_rad.set_cavity_phase(method='tracking')
    ring_av = ring_rad.deepcopy()
    ring_rad.tapering(niter = 2, multipoles=True, XYStep=xy_step, DPStep=dp_step)
    l0t,qt,qpt,lt = linopt_rad(ring_rad, refpts=range(len(ring_rad)), get_chrom=True, coupled=False, XYStep=xy_step, DPStep=dp_step)

    print('Computing time:', time.time()-t0, 's')
    
    # visualising orbit and PolynomB
    o0a, _ = at.find_orbit6(ring_av, XYStep=xy_step, DPStep=dp_step)
    o6a = np.squeeze(lattice_pass(ring_av, o0a, refpts=range(len(ring))))
    
    o0, _ = find_orbit6(ring_rad, XYStep=xy_step, DPStep=dp_step)
    o6 = np.squeeze(lattice_pass(ring_rad, o0, refpts=range(len(ring))))
    spos = ring_rad.get_s_pos(range(len(ring_rad)))
    
    ring_av.tapering(niter = 2, multipoles=True, XYStep=xy_step, DPStep=dp_step)
    div_number = 20
    average_tapering(ring_av, multipoles=False, XYStep=xy_step, DPStep=dp_step, div_number = div_number)
    
    o0b, _ = at.find_orbit6(ring_av, XYStep=xy_step, DPStep=dp_step)
    o6b = np.squeeze(lattice_pass(ring_av, o0b, refpts=range(len(ring))))
    l0t_av,qt_av,qpt_av,lt_av = linopt_rad(ring_av, refpts=range(len(ring_av)), get_chrom=True, coupled=False, XYStep=xy_step, DPStep=dp_step)
    
    fig1, ax1 = fcc_axes()
    ax1.plot(spos, o6[0]*1e6, color = 'orange', label='Individual tapering')
    ax1.plot(spos, o6a[0]*1e6, color = 'navy', label='No tapering')
    ax1.plot(spos, o6b[0]*1e6, color = 'red', label='Average tapering')
    ax1.set_xlabel('s [m]')
    ax1.set_ylabel(r'x [$\mu$m]')
    ax1.legend(loc = 1, prop={'size': 12})
    plt.show()

    fig2, ax2 = fcc_axes()
    ax2.plot(spos,(lt_av.beta[:,0]-ll.beta[:,0])/ll.beta[:,0], label = 'PyAT', c='navy')
    ax2.set_xlabel('s [m]')
    ax2.set_ylabel(r'$\frac{\Delta \beta_x}{\beta_x^{Ref}}$')
    ax2.legend(loc = 1, prop={'size': 8})
    plt.show()

    fig3, ax3 = fcc_axes()
    ax3.plot(spos, lt_av.dispersion[:,0] - ll.dispersion[:,0], label = 'PyAT', c='navy')
    ax3.set_ylabel(r'$\eta_x^{Average} - \eta_x^{Ref}$')
    ax3.set_xlabel('s [m]')
    ax3.legend(loc = 1, prop={'size': 8})
    plt.show()

    # Emittance, beta and dispersion values at IP1 and IP2
    env = ring_rad.envelope_parameters()
    epsilon_x = env.emittances[0]
    epsilon_y = env.emittances[1]
    env = ring_av.envelope_parameters()
    epsilon_x_av = env.emittances[0]
    epsilon_y_av = env.emittances[1]
    
    print('Emittance_x for individual and average:', epsilon_x, epsilon_x_av)
    
    IP_indexes = get_refpts(ring_av, 'IP*')
    refb1_x = (lt_av.beta[0,0]-ll.beta[0,0])/ll.beta[0,0]
    refb2_x = (lt_av.beta[IP_indexes[1],0]-ll.beta[IP_indexes[1],0])/ll.beta[IP_indexes[1],0]
    print('Betas at IP1 and IP2:', refb1_x, refb2_x)
    
    refd1_x = (lt_av.dispersion[0,0]-ll.dispersion[0,0])
    refd2_x = (lt_av.dispersion[IP_indexes[1],0]-ll.dispersion[IP_indexes[1],0])
    print('Relative dispersions at IP1 and IP2:', refd1_x, refd2_x)
