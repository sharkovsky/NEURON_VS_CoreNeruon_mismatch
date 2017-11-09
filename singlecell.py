from __future__ import print_function
import sys
import os
import numpy as np

os.environ['HOC_LIBRARY_PATH'] = str( os.path.join( os.getcwd(), 'L5_TTPC1_cADpyr232_1' ) )

from neuron import h, gui
h.load_file("import3d.hoc")
h.load_file("morphology.hoc")
h.load_file("biophysics.hoc")
h.load_file("template.hoc")

pc = h.ParallelContext()
rank = int(pc.id())
nhost = int(pc.nhost())

pc.nthread(1,0)

simulate = False
write_crn_files = True

# Reset NEURON
#pc.gid_clear()
#syndict = None
#cell = None
#stimlist = None

h.Random().Random123_globalindex(123456 + rank)
np.random.seed(42 + rank)

h.cvode.cache_efficient(1)

junction_potential = -14.0
h.tstop = 50.
h.celsius = 34.0

cell = neuron.h.cADpyr232_L5_TTPC1_0fb1ca4724(0) #0 means no synapses :/
pc.set_gid2node(1, rank)
nc = h.NetCon(cell.soma[0](1)._ref_v, None, sec=cell.soma[0])
nc.threshold = -20.
pc.cell(1, nc)

tot_nseg = 0
tot_nsec = 0
for sec in h.allsec():
    tot_nseg += sec.nseg
    tot_nsec += 1
print( 'Before calling geom_nseg' )
print( 'tot nsec:', tot_nsec, 'tot nseg:', tot_nseg)

cell.geom_nseg()

tot_nseg = 0
tot_nsec = 0
for sec in h.allsec():
    tot_nseg += sec.nseg
    tot_nsec += 1
print( 'After calling geom_nseg' )
print( 'tot nsec:', tot_nsec, 'tot nseg:', tot_nseg)



tvec = h.Vector()
idvec = h.Vector()
pc.spike_record(-1,tvec,idvec)

pc.set_maxstep(10)
h.stdinit()

if write_crn_files:
    bbcorewrite_folder = 'coredat/single_cell'
    if rank == 0:
        if not os.path.exists(bbcorewrite_folder):
            os.makedirs(bbcorewrite_folder)

    pc.barrier()
    pc.nrnbbcore_write(bbcorewrite_folder)

h.quit()
