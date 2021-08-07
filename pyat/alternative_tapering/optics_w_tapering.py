import at
from at import load_mat
from at.physics import linopt, linopt_rad
import matplotlib.pyplot as plt
import numpy as np

path = './Lattices/'
filename = 'fcch_norad.mat'
key = 'ring'
latticef = path+filename
ring = at.load.matfile.load_mat(latticef, mat_key=key)
ring = at.Lattice(ring)

xy_step = 1.0e-9
dp_step = 1.0e-6

# optics of the original lattice
ring.radiation_off()
l0, q, qp, l = linopt(ring,dp=0,refpts=range(len(ring)),get_chrom = True, coupled=False, XYStep=xy_step, DPStep=dp_step)

# optics of tapered lattice
ring.radiation_on(quadrupole_pass='auto')
ring.set_cavity_phase()
ring.tapering(niter = 2, quadrupole=True, sextupole=True, XYStep=xy_step, DPStep=dp_step)
l0t,qt,qpt,lt = linopt_rad(ring, refpts=range(len(ring)), get_chrom=True, coupled=False, XYStep=xy_step, DPStep=dp_step)

# s along the track
spos = ring.get_s_pos(range(len(ring)))

# optics comparison
plt.figure(figsize=(12,6))
plt.subplot(121)
plt.plot(spos,(lt.beta[:,0]-l.beta[:,0])/l.beta[:,0], label='pyAT', c='blue')
plt.xlabel('s [m]')
plt.ylabel(r'$\frac{\Delta \beta_x}{\beta_x^Ref}$')
plt.legend(loc = 1, prop={'size': 8})

plt.subplot(122)
plt.plot(spos,(lt.beta[:,1]-l.beta[:,1])/l.beta[:,1], label='pyAT', c='blue')
plt.xlabel('s [m]')
plt.ylabel(r'$\frac{\Delta \beta_y}{\beta_y^Ref}$')
plt.legend(loc = 1, prop={'size': 8})
plt.show()

plt.figure(figsize=(12,6))
plt.plot(spos, lt.dispersion[:,0] - l.dispersion[:,0], label='pyAT', c='blue')
plt.xlabel('s [m]')
plt.ylabel(r'$\Delta \eta_x$')
plt.legend(loc = 1, prop={'size': 8})
plt.show()
