import at
import numpy as np
import matplotlib.pyplot as plt
from at.tracking import element_pass, lattice_pass
from at.physics import linopt, fast_ring, linopt_rad
from at.lattice import DConstant, PyElement, check_radiation
from at import load_mat, get_refpts, get_elements, get_value_refpts, set_value_refpts, elements, find_orbit6
import time
import sys

@check_radiation(True)
def average_tapering(ring, ring_rad, quadrupole=True, sextupole=True, niter=1, **kwargs):
    """
    Delivers a scaled and averaged tapering for a predefined set of dipole families. Modifies the ring.
    PARAMETERS
        ring            lattice description
        ring_rad        individually tapered lattice
    KEYWORDS
        family          family of dipoles, set of respective indexes
        qadrupole=True  scale quadrupole fields
        sextupole=True  scale sextupole fields
        XYStep=1.0e-8   transverse step for numerical computation
        DPStep=1.0E-6   momentum deviation used for computation of orbit6
    """
    xy_step = kwargs.pop('XYStep', DConstant.XYStep)
    dp_step = kwargs.pop('DPStep', DConstant.DPStep)
    fam = kwargs.pop('family')
    dipole_indexes = get_refpts(ring, elements.Dipole)
    
    for dip in fam:
        for dipoles in dip:
            k_arc = get_value_refpts(ring_rad, dipoles, 'PolynomB', index=0)
            k0_arc =  get_value_refpts(ring, dipoles, 'BendingAngle')/get_value_refpts(ring, dipoles, 'Length')
            factor = np.average(k_arc/k0_arc, weights=None)
            set_value_refpts(ring, dipoles, 'PolynomB', factor*k0_arc, index=0)
            set_value_refpts(ring, dipoles, 'PolynomB', 0.0 , index=1)
    
    if quadrupole:
        quadrupoles = get_refpts(ring, elements.Quadrupole)
        k01 = get_value_refpts(ring, quadrupoles, 'PolynomB', index=1)
        o0, _ = find_orbit6(ring, XYStep=xy_step, DPStep=dp_step)
        o6 = np.squeeze(lattice_pass(ring, o0, refpts=range(len(ring))))
        dpps = (o6[4, quadrupoles] + o6[4, quadrupoles+1]) / 2
        set_value_refpts(ring, quadrupoles, 'PolynomB', k01*(1+dpps), index=1)

    if sextupole:
        sextupoles = get_refpts(ring, elements.Sextupole)
        k02 = get_value_refpts(ring, sextupoles, 'PolynomB', index=2)
        o0, _ = find_orbit6(ring, XYStep=xy_step, DPStep=dp_step)
        o6 = np.squeeze(lattice_pass(ring, o0, refpts=range(len(ring))))
        dpps = (o6[4, sextupoles] + o6[4, sextupoles+1]) / 2
        set_value_refpts(ring, sextupoles, 'PolynomB', k02*(1+dpps), index=2)
        
def family_segmentation(ring):
    """
    Currently divided into 4 groups as defined by the position of IPs and RFs, additional subdivisin in 'div_number' near-equal subgroups.
    """
    
    div_number = 55
    dipole_indexes = get_refpts(ring, elements.Dipole)
    IP_indexes = get_refpts(ring, 'IP*')
    RF_indexes = get_refpts(ring, elements.RFCavity)
    
    print('Number of dipoles per family:', len(dipole_indexes)/(4*div_number))
    print('Number of dipoles:', len(dipole_indexes))

    ind1 = int(np.argwhere(dipole_indexes < RF_indexes[0])[-1])
    Arc1 = dipole_indexes[0:ind1]
    Arc1 = np.array_split(Arc1, div_number) # allows non-equal splitting

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

if __name__ == "__main__":
    # load lattices, step sizes adjusted to avoid systematics on lattice calculations
    t0 = time.time()
    xy_step = 1.0e-9
    dp_step = 1.0e-6

    ring = at.load_lattice('./Lattices/fcch_norad.mat', mat_key='ring')
    ring = at.Lattice(ring)
    
    # optics of the original lattice
    ring_rad = ring.deepcopy()
    ring_rad.radiation_on(quadrupole_pass='auto')
    ring_rad.set_cavity_phase(method='tracking')
    ring_av = ring_rad.deepcopy()
    ring_rad.tapering(niter = 2, quadrupole=True, sextupole=True, XYStep=xy_step, DPStep=dp_step)

    print('Computing time:', time.time()-t0, 's')
    
    # visualising orbit and PolynomB
    o0a, _ = at.find_orbit6(ring_av, XYStep=xy_step, DPStep=dp_step)
    o6a = np.squeeze(lattice_pass(ring_av, o0a, refpts=range(len(ring))))
    
    o0, _ = find_orbit6(ring_rad, XYStep=xy_step, DPStep=dp_step)
    o6 = np.squeeze(lattice_pass(ring_rad, o0, refpts=range(len(ring))))
    spos = ring_rad.get_s_pos(range(len(ring_rad)))
    
    Arcs = family_segmentation(ring_av)
    average_tapering(ring_av, ring_rad,quadrupole=True, sextupole=True, family = Arcs, XYStep=xy_step, DPStep=dp_step)
    o0b, _ = at.find_orbit6(ring_av, XYStep=xy_step, DPStep=dp_step)
    o6b = np.squeeze(lattice_pass(ring_av, o0b, refpts=range(len(ring))))
    
    plt.figure()
    plt.plot(spos, o6[0], color = 'orange', label='Individual taper')
    plt.plot(spos, o6a[0], color = 'navy', label='No taper')
    plt.plot(spos, o6b[0], color = 'red', label='Average taper')
    plt.xlabel('s [m]')
    plt.ylabel('x [m]')
    plt.legend(loc = 1)
    plt.show()
    
    dipole_indexes = get_refpts(ring_rad, elements.Dipole)
    sdip = ring_rad.get_s_pos(dipole_indexes)
    b = get_value_refpts(ring_rad, dipole_indexes, 'BendingAngle')
    l = get_value_refpts(ring_rad, dipole_indexes, 'Length')
    B0 = get_value_refpts(ring_rad, dipole_indexes, 'PolynomB', index = 0)
    
    b_av = get_value_refpts(ring_av, dipole_indexes, 'BendingAngle')
    l_av = get_value_refpts(ring_av, dipole_indexes, 'Length')
    B0_av = get_value_refpts(ring_av, dipole_indexes, 'PolynomB', index = 0)
    
    plt.plot(sdip, B0/(b/l), color='navy', label='Individual')
    plt.plot(sdip, B0_av/(b_av/l_av), color='orange', label = 'Average')
    plt.xlabel('s [m]')
    plt.ylabel(r'$\frac{B_0}{Angle / Length}$')
    plt.legend(loc=1)
    plt.show()
