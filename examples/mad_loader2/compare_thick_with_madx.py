import xtrack as xt
import xpart as xp
from cpymad.madx import Madx
from xtrack.mad_loader import MadLoader, TeapotSlicing
import matplotlib.pyplot as plt


# Make thin using madx
mad1 = Madx()
mad1.call('lhc_sequence.madx')
mad1.call('lhc_optics.madx')
mad1.beam()
mad1.sequence.lhcb1.use()

tw_mad = mad1.twiss()

# Do the same in xsuite
mad2 = Madx()
mad2.call('lhc_sequence.madx')
mad2.call('lhc_optics.madx')
mad2.beam()
mad2.sequence.lhcb1.use()

ml = MadLoader(mad2.sequence.lhcb1, enable_slicing=False)

ml.slicing_strategies = [
    ml.make_slicing_strategy(
        name_regex=r'(mqt|mqtli|mqtlh)\..*',
        slicing_strategy=TeapotSlicing(2),
    ),
    ml.make_slicing_strategy(
        name_regex=r'(mbx|mbrb|mbrc|mbrs|mbh|mqwa|mqwb|mqy|mqm|mqmc|mqml)\..*',
        slicing_strategy=TeapotSlicing(4),
    ),
    ml.make_slicing_strategy(
        madx_type='mqxb',
        slicing_strategy=TeapotSlicing(16),
    ),
    ml.make_slicing_strategy(
        madx_type='mqxa',
        slicing_strategy=TeapotSlicing(16),
    ),
    ml.make_slicing_strategy(
        madx_type='mq',
        slicing_strategy=TeapotSlicing(2),
    ),
    ml.make_slicing_strategy(
        madx_type='mb',
        slicing_strategy=TeapotSlicing(2),
    ),
    ml.make_slicing_strategy(TeapotSlicing(1)),  # Default catch-all as in MAD-X
]
line_xt = ml.make_line()
line_xt.cycle(name_first_element='ip3')
line_xt.particle_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, q0=1, energy0=7e12)

line_xt.build_tracker()
line_xt.optimize_for_tracking()
tw_xt = line_xt.twiss(method='4d', strengths=True)


print(f'madx tunes = {tw_mad.summary.q1, tw_mad.summary.q2}')
print(f'xt tunes = {tw_xt.qx, tw_xt.qy}')



def plot_coord_for_elements(coord):
    fig, ax = plt.subplots()
    ax.plot(tw_xt.s, tw_xt[coord], '.-', color='r')
    ax.plot(tw_mad.s, tw_mad[coord], '.-', color='b')

    for i, txt in enumerate(tw_xt.name):
        if txt.startswith('drift'):
            continue
        ax.annotate(txt, (tw_xt.s[i], tw_xt[coord][i]), color='r')

    for i, txt in enumerate(tw_mad.name):
        if txt.startswith('drift'):
            continue
        ax.annotate(txt, (tw_mad.s[i], tw_mad[coord][i]), color='b')

    fig.show()


# Compare x along the lines thinned in madx and xsuite:
# plot_coord_for_elements('x')
