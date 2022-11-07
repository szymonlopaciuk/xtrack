import json

# TODO:
# - Put lag on the stable slope
# - Assert that calculated enekick is smaller than voltage at each cavity
# - Check 6d kick
# - Assert no collective
# - Assert no ions

# Some considerations:
# - I observe the right tune on v on the real tracker and not with the one with the
# eneloss of the closed orbit. --> forget, works only when orbit is zero
# On the real tracker I see that the beta beating is introduced by the cavities
# not by the multipoles --> forget, works only when orbit is zero

import numpy as np
from scipy.constants import c as clight
import xtrack as xt

with open('line_no_radiation.json', 'r') as f:
    line = xt.Line.from_dict(json.load(f))

line[3].knl[0] += 1e-6
line[3].ksl[0] += 1e-6

#line_no_rad = line.copy()

line_df = line.to_pandas()
multipoles = line_df[line_df['element_type'] == 'Multipole']
cavities = line_df[line_df['element_type'] == 'Cavity'].copy()

# save voltages
cavities['voltage'] = [cc.voltage for cc in cavities.element.values]
cavities['frequency'] = [cc.frequency for cc in cavities.element.values]
cavities['eneloss_partitioning'] = cavities['voltage'] / cavities['voltage'].sum()


# set voltages to zero
for cc in cavities.element.values:
    cc.voltage = 0

tracker = xt.Tracker(line = line)
tw_no_rad = tracker.twiss(method='4d')

p_test = tw_no_rad.particle_on_co.copy()
tracker.configure_radiation(mode='mean')

tracker.track(p_test, turn_by_turn_monitor='ONE_TURN_EBE')
mon = tracker.record_last_track

n_cavities = len(cavities)

tracker_taper = xt.Tracker(line = line, extra_headers=["#define XTRACK_MULTIPOLE_TAPER"])

import matplotlib.pyplot as plt
plt.close('all')

rtot_eneloss = 1e-10

# Put all cavities on crest and at zero frequency
for cc in cavities.element.values:
    cc.lag = 90
    cc.frequency = 0

while True:
    p_test = tw_no_rad.particle_on_co.copy()
    tracker_taper.configure_radiation(mode='mean')
    tracker_taper.track(p_test, turn_by_turn_monitor='ONE_TURN_EBE')
    mon = tracker_taper.record_last_track

    eloss = -(mon.ptau[0, -1] - mon.ptau[0, 0]) * p_test.p0c[0]
    print(f"Energy loss: {eloss:.3f} eV")

    if eloss < p_test.energy0[0]*rtot_eneloss:
        break

    for ii in cavities.index:
        cc = cavities.loc[ii, 'element']
        eneloss_partitioning = cavities.loc[ii, 'eneloss_partitioning']
        cc.voltage += eloss * eneloss_partitioning

    plt.plot(mon.s.T, mon.ptau.T)

delta_beta_corr = mon.delta[0, :]

i_multipoles = multipoles.index.values
delta_taper = ((mon.delta[0,:][i_multipoles+1] + mon.delta[0,:][i_multipoles]) / 2)
for nn, dd in zip(multipoles['name'].values, delta_taper):
    line[nn].knl *= (1 + dd)
    line[nn].ksl *= (1 + dd)

beta0 = p_test.beta0[0]
v_ratio = []
for icav in cavities.index:
    v_ratio.append(cavities.loc[icav, 'element'].voltage / cavities.loc[icav, 'voltage'])
    inst_phase = np.arcsin(cavities.loc[icav, 'element'].voltage / cavities.loc[icav, 'voltage'])
    freq = cavities.loc[icav, 'frequency']

    zeta = mon.zeta[0, icav]
    lag = 360.*(inst_phase / (2*np.pi) - freq*zeta/beta0/clight)
    lag = 180. - lag # we are above transition

    cavities.loc[icav, 'element'].lag = lag
    cavities.loc[icav, 'element'].frequency = freq
    cavities.loc[icav, 'element'].voltage = cavities.loc[icav, 'voltage']

tracker.configure_radiation(mode='mean')
tw_not_symplectic = tracker.twiss(method='6d', matrix_stability_tol=0.5,
                    eneloss_and_damping=True) # Completely wrong in y when
                                              # closed orbit is not zero

tracker_sympl = xt.Tracker(line = line, extra_headers=["#define XSUITE_SYNRAD_SAME_AS_FIRST"])
tracker_sympl.configure_radiation(mode='mean')
for ee in line.elements:
    if hasattr(ee, 'rescale_pxpy'):
        ee.rescale_pxpy = 1
tw = tracker_sympl.twiss(method='6d', matrix_stability_tol=0.5)



print('Non sympltectic tracker:')
print(f'Tune error =  error_qx: {abs(tw_not_symplectic.qx - tw_no_rad.qx):.3e} error_qy: {abs(tw_not_symplectic.qy - tw_no_rad.qy):.3e}')
print('Sympltectic tracker:')
print(f'Tune error =  error_qx: {abs(tw.qx - tw_no_rad.qx):.3e} error_qy: {abs(tw.qy - tw_no_rad.qy):.3e}')
plt.figure(2)

plt.subplot(2,1,1)
plt.plot(tw_no_rad.s, tw.betx/tw_no_rad.betx - 1)
#tw.betx *= (1 + delta_beta_corr)
#plt.plot(tw_no_rad.s, tw.betx/tw_no_rad.betx - 1)
plt.ylabel(r'$\Delta \beta_x / \beta_x$')

plt.subplot(2,1,2)
plt.plot(tw_no_rad.s, tw.bety/tw_no_rad.bety - 1)
#tw.bety *= (1 + delta_beta_corr)
#plt.plot(tw_no_rad.s, tw.bety/tw_no_rad.bety - 1)
plt.ylabel(r'$\Delta \beta_y / \beta_y$')

plt.figure(10)
plt.subplot(2,1,1)
plt.plot(tw_no_rad.s, tw_no_rad.x, 'k')
plt.plot(tw_no_rad.s, tw_not_symplectic.x, 'b')

plt.subplot(2,1,2)
plt.plot(tw_no_rad.s, tw_no_rad.y, 'k')
plt.plot(tw_no_rad.s, tw_not_symplectic.y, 'b')


plt.show()