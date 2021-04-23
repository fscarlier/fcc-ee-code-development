import os, sys
import at
import numpy as np
from cpymad.madx import Madx


def initiate_cpymad(lattice, energy):
    madx = Madx()
    madx.call(file=lattice)
    madx.input('SET, FORMAT="25.20e";')
    madx.command.beam(particle='electron', energy=energy)
    madx.use('l000013')
    return madx


def convert_marker(ele):
    marker = at.lattice.elements.Marker(ele.name.replace('.', '_').upper())
    return marker


def convert_drift(ele):
    drift = at.lattice.elements.Drift(ele.name.replace('.', '_').upper(), ele.l)
    return drift


def convert_dipole(ele):
    try:
        arc_length = (ele.angle*ele.l)/(2*np.sin(ele.angle/2.))
        dipole = at.lattice.elements.Dipole(ele.name.replace('.', '_').upper(), arc_length, bending_angle=ele.angle)
        return dipole
    except AttributeError:
        print('Attribte angle not present in {}'.format(ele.name))


def convert_quadrupole(ele):
    try:
        quad = at.lattice.elements.Quadrupole(ele.name.replace('.', '_').upper(), ele.l, PolynomB=[0.0,ele.k1], PolynomA=[0.0,ele.k1s])
        return quad
    except AttributeError:
        print('Attribte k1 or k1s not present in {}'.format(ele.name))


def convert_sextupole(ele):
    try:
        sext = at.lattice.elements.Sextupole(ele.name.replace('.', '_').upper(), ele.l, PolynomB=[0.0,0.0,ele.k2], PolynomA=[0.0,0.0,ele.k2s])
        return sext
    except AttributeError:
        print('Attribte k2 or k2s not present in {}'.format(ele.name))


def convert_rfcavity(ele, energy):
    try:
        rfcavity = at.lattice.elements.RFCavity(ele.name.replace('.', '_').upper(), ele.l, voltage=ele.volt, frequency=ele.freq*1e6, harmonic_number=ele.harmon, energy=energy)
        return rfcavity
    except AttributeError:
        print('Wrong attribute for element {}'.format(ele.name))


def convert_madx_lattice_pyat(lattice_path, energy):
    madx = initiate_cpymad(lattice_path, energy)
    sequences = []
    for key in madx.sequence:
        seq = madx.sequence[key].elements
        pyat_seq = []
        previous_at = 0
        previous_l = 0
        drift_count = 0
        for ele in seq:
            if ele.base_type.name == 'marker':
                marker = convert_marker(ele)
                pyat_seq.append(marker)
            elif ele.base_type.name == 'drift':
                drift = convert_drift(ele)
                pyat_seq.append(drift)
            elif ele.base_type.name == 'rbend':
                dipole = convert_dipole(ele)
                pyat_seq.append(dipole)
            elif ele.base_type.name == 'quadrupole':
                quadrupole = convert_quadrupole(ele)
                pyat_seq.append(quadrupole)
            elif ele.base_type.name == 'sextupole':
                sextupole = convert_sextupole(ele)
                pyat_seq.append(sextupole)
            elif ele.base_type.name == 'rfcavity':
                rfcavity = convert_rfcavity(ele, energy)
                pyat_seq.append(rfcavity)
            else:
                print('This element with base type ({}) is not converted: {}'.format(ele.base_type, ele.name))
            
            if ele.base_type.name == 'rbend':
                arc_length = (ele.angle*ele.l)/(2*np.sin(ele.angle/2.))
                drift_length = (ele.at - previous_at) - (arc_length/2. + previous_l/2.)
            else: 
                drift_length = (ele.at - previous_at) - (ele.l/2. + previous_l/2.)
            
            if drift_length != 0.0:
                drift = at.lattice.elements.Drift('Drift{}'.format(drift_count), drift_length)
                pyat_seq.insert(-1,drift)
                drift_count += 1

            previous_at = ele.at
            if ele.base_type.name == 'rbend':
                arc_length = (ele.angle*ele.l)/(2*np.sin(ele.angle/2.))
                previous_l = arc_length
            else: 
                previous_l = ele.l

        #TODO: Return multiple lattices when multiple keys are present
        ring = at.Lattice(pyat_seq, key=key, energy=energy*1e9)
        return ring

