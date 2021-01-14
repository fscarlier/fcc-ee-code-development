import at
import os
import matplotlib.pyplot as plt
import numpy as np
import time
import at.plot.specific

path = '/machfs/swhite/FCC/madx_lattices/'
filename = 'fcch_norad.mat'
key = 'ring'
latticef = path+filename
ring = at.load.matfile.load_mat(latticef,key=key)
ring = at.Lattice(ring)

l0, q, qp, ld = at.linopt(ring,refpts=range(len(ring)),coupled=False)
ring.plot_beta()
plt.show()
