import numpy as np
from cpymad.madx import Madx

import xtrack as xt
import xpart as xp

mad = Madx()
mad.call('../../test_data/hllhc15_noerrors_nobb/sequence.madx')
mad.use('lhcb1')

tw_mad_no_coupling = mad.twiss(ripken=True).dframe()

# introduce coupling
mad.sequence.lhcb1.expanded_elements[7].ksl = [0,1e-4]

tw_mad_coupling = mad.twiss(ripken=True).dframe()

line = xt.Line.from_madx_sequence(mad.sequence.lhcb1)
line.particle_ref = xp.Particles(p0c=7000e9, mass0=xp.PROTON_MASS_EV)

tracker = line.build_tracker()

tw = tracker.twiss()

Ws = np.array(tw.W_matrix)

bety1 = tw.bety1
betx2 = tw.betx2

cmin = tw.c_minus

assert np.isclose(cmin, mad.table.summ.dqmin[0], rtol=0, atol=1e-5)

import matplotlib.pyplot as plt
plt.close('all')
plt.figure(1)
sp1 = plt.subplot(211)
plt.plot(tw.s, tw.bety1, label='bety1')
plt.plot(tw_mad_coupling.s, tw_mad_coupling.beta21, '--')
plt.ylabel(r'$\beta_{1,y}$')
plt.subplot(212, sharex=sp1)
plt.plot(tw.s, tw.betx2, label='betx2')
plt.plot(tw_mad_coupling.s, tw_mad_coupling.beta12, '--')
plt.ylabel(r'$\beta_{2,x}$')
plt.suptitle(r'Xsuite: $C^{-}$'
             f" = {cmin:.2e} "
             r"MAD-X: $C^{-}$ = "
             f"{mad.table.summ.dqmin[0]:.2e}"
             )
plt.show()