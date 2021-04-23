def average_tapering(ring, ring_rad, quadrupole=True, sextupole=True, **kwargs):
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
    
    for dipoles in fam:
        k_arc = get_value_refpts(ring_rad, dipoles, 'PolynomB')[:,0]
        k0_arc =  get_value_refpts(ring, dipoles, 'BendingAngle')/get_value_refpts(ring, dipoles, 'Length')
        factor = np.average(k_arc/k0_arc, weights=None)
        set_value_refpts(ring, dipoles, 'PolynomB', factor*k0_arc, index=0)
    
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
    Not a smart way of doing it, update needed. Currently divided into 'div_number' near-equal groups.
    """
    div_number = 100
    dipole_indexes = get_refpts(ring, elements.Dipole)
    Arcs = np.array_split(dipole_indexes, div_number) # allows non-equal splitting
    
    return Arcs
    
import at
import numpy as np
import matplotlib.pyplot as plt
from at.tracking import element_pass, lattice_pass
from at.physics import linopt, fast_ring, linopt_rad
from at.lattice import DConstant, PyElement
from at import load_mat, get_refpts, get_elements, get_value_refpts, set_value_refpts, elements
import time

if __name__ == "__main__":
    # load lattices, step sizes adjusted to avoid systematics on lattice calculations
    t0 = time.time()
    xy_step = 1.0e-9
    dp_step = 1.0e-6

    ring = at.load_lattice('./Lattices/fcch_norad.mat', mat_key='ring')
    ring = at.Lattice(ring)
    ring.radiation_off()
    spos = ring.get_s_pos(range(len(ring)))

    ring_rad = ring.deepcopy()
    ring_rad.radiation_on(quadrupole_pass='auto')
    ring_rad.set_cavity_phase(method='tracking')
    ring_av = ring_rad.deepcopy()
    ring_rad.tapering(niter = 2, quadrupole=True, sextupole=True, XYStep=xy_step, DPStep=dp_step)
    
    # orbit
    o0a, _ = at.find_orbit6(ring_av, XYStep=xy_step, DPStep=dp_step)
    o6a = np.squeeze(lattice_pass(ring_av, o0a, refpts=range(len(ring))))
    
    o0, _ = find_orbit6(ring_rad, XYStep=xy_step, DPStep=dp_step)
    o6 = np.squeeze(lattice_pass(ring_rad, o0, refpts=range(len(ring))))
    spos = ring_rad.get_s_pos(range(len(ring_rad)))
    
    Arcs = family_segmentation(ring_av)
    average_tapering(ring_av, ring_rad,quadrupole=True, sextupole=True, family = Arcs, XYStep=xy_step, DPStep=dp_step)
    o0b, _ = at.find_orbit6(ring_av, XYStep=xy_step, DPStep=dp_step)
    o6b = np.squeeze(lattice_pass(ring_av, o0b, refpts=range(len(ring))))
    
    print('Computing time:', time.time()-t0, 's')
    
    # visualise
    plt.figure()
    plt.plot(spos, o6[0], color = 'orange', label='Individual taper')
    plt.plot(spos, o6a[0], color = 'navy', label='No taper')
    plt.plot(spos, o6b[0], color = 'red', label='Average taper')
    plt.xlabel('s [m]')
    plt.ylabel('x [m]')
    plt.legend(loc = 1)
    plt.show()
