import xtrack as xt

collider = xt.Multiline.from_json(
    '../../test_data/hllhc15_thick/hllhc15_collider_thick.json')
collider.build_trackers()

nrj = 7000.
scale = 23348.89927*0.9
scmin = 0.03*7000./nrj
qtlimitx28 = 1.0*225.0/scale
qtlimitx15 = 1.0*205.0/scale
qtlimit2 = 1.0*160.0/scale
qtlimit3 = 1.0*200.0/scale
qtlimit4 = 1.0*125.0/scale
qtlimit5 = 1.0*120.0/scale
qtlimit6 = 1.0*90.0/scale

collider.vars.vary_default.update({
    'kq4.r8b1':    {'step': 1.0E-6, 'limits': ( qtlimit2*scmin, qtlimit2)},
    'kq5.r8b1':    {'step': 1.0E-6, 'limits': (-qtlimit2, qtlimit2*scmin)},
    'kq6.r8b1':    {'step': 1.0E-6, 'limits': ( qtlimit2*scmin, qtlimit2)},
    'kq7.r8b1':    {'step': 1.0E-6, 'limits': (-qtlimit3, qtlimit3*scmin)},
    'kq8.r8b1':    {'step': 1.0E-6, 'limits': ( qtlimit3*scmin, qtlimit3)},
    'kq9.r8b1':    {'step': 1.0E-6, 'limits': (-qtlimit3, qtlimit3*scmin)},
    'kq10.r8b1':   {'step': 1.0E-6, 'limits': ( qtlimit3*scmin, qtlimit3)},
    'kqtl11.r8b1': {'step': 1.0E-6, 'limits': (-qtlimit4, qtlimit4)},
    'kqt12.r8b1':  {'step': 1.0E-6, 'limits': (-qtlimit5, qtlimit5)},
    'kqt13.r8b1':  {'step': 1.0E-6, 'limits': (-qtlimit5, qtlimit5)},
    'kq4.l8b1':    {'step': 1.0E-6, 'limits': (-qtlimit2, qtlimit2*scmin)},
    'kq5.l8b1':    {'step': 1.0E-6, 'limits': ( qtlimit2*scmin, qtlimit2)},
    'kq6.l8b1':    {'step': 1.0E-6, 'limits': (-qtlimit3, qtlimit3*scmin)},
    'kq7.l8b1':    {'step': 1.0E-6, 'limits': ( qtlimit3*scmin, qtlimit3)},
    'kq8.l8b1':    {'step': 1.0E-6, 'limits': (-qtlimit3, qtlimit3*scmin)},
    'kq9.l8b1':    {'step': 1.0E-6, 'limits': ( qtlimit3*scmin, qtlimit3)},
    'kq10.l8b1':   {'step': 1.0E-6, 'limits': (-qtlimit3, qtlimit3*scmin)},
    'kqtl11.l8b1': {'step': 1.0E-6, 'limits': (-qtlimit4, qtlimit4)},
    'kqt12.l8b1':  {'step': 1.0E-6, 'limits': (-qtlimit5, qtlimit5)},
    'kqt13.l8b1':  {'step': 1.0E-6, 'limits': (-qtlimit5, qtlimit5)},
    'kq4.r8b2':    {'step': 1.0E-6, 'limits': (-qtlimit2, qtlimit2*scmin)},
    'kq5.r8b2':    {'step': 1.0E-6, 'limits': ( qtlimit2*scmin, qtlimit2)},
    'kq6.r8b2':    {'step': 1.0E-6, 'limits': (-qtlimit2, qtlimit2*scmin)},
    'kq7.r8b2':    {'step': 1.0E-6, 'limits': ( qtlimit3*scmin, qtlimit3)},
    'kq8.r8b2':    {'step': 1.0E-6, 'limits': (-qtlimit3, qtlimit3*scmin)},
    'kq9.r8b2':    {'step': 1.0E-6, 'limits': ( qtlimit3*scmin, qtlimit3)},
    'kq10.r8b2':   {'step': 1.0E-6, 'limits': (-qtlimit3, qtlimit3*scmin)},
    'kqtl11.r8b2': {'step': 1.0E-6, 'limits': (-qtlimit4, qtlimit4)},
    'kqt12.r8b2':  {'step': 1.0E-6, 'limits': (-qtlimit5, qtlimit5)},
    'kqt13.r8b2':  {'step': 1.0E-6, 'limits': (-qtlimit5, qtlimit5)},
    'kq5.l8b2':    {'step': 1.0E-6, 'limits': (-qtlimit2, qtlimit2*scmin)},
    'kq4.l8b2':    {'step': 1.0E-6, 'limits': ( qtlimit2*scmin, qtlimit2)},
    'kq6.l8b2':    {'step': 1.0E-6, 'limits': ( qtlimit2*scmin, qtlimit2)},
    'kq7.l8b2':    {'step': 1.0E-6, 'limits': (-qtlimit3, qtlimit3*scmin)},
    'kq8.l8b2':    {'step': 1.0E-6, 'limits': ( qtlimit3*scmin, qtlimit3)},
    'kq9.l8b2':    {'step': 1.0E-6, 'limits': (-qtlimit3, qtlimit3*scmin)},
    'kq10.l8b2':   {'step': 1.0E-6, 'limits': ( qtlimit3*scmin, qtlimit3)},
    'kqtl11.l8b2': {'step': 1.0E-6, 'limits': (-qtlimit4, qtlimit4)},
    'kqt12.l8b2':  {'step': 1.0E-6, 'limits': (-qtlimit5, qtlimit5)},
    'kqt13.l8b2':  {'step': 1.0E-6, 'limits': (-qtlimit5, qtlimit5)},
})
mux_b1_target = 3.0985199176526272
muy_b1_target = 2.7863674079923726

collider.varval['kq6.l8b1'] *= 1.1
collider.varval['kq6.r8b1'] *= 1.1

tab_boundary_right = collider.lhcb1.twiss(
    ele_start='ip8', ele_stop='ip1.l1',
    twiss_init=xt.TwissInit(element_name='ip1.l1', line=collider.lhcb1,
                            betx=0.15, bety=0.15))
tab_boundary_left = collider.lhcb1.twiss(
    ele_start='ip5', ele_stop='ip8',
    twiss_init=xt.TwissInit(element_name='ip5', line=collider.lhcb1,
                            betx=0.15, bety=0.15))

opt = collider[f'lhcb1'].match(
    default_tol={None: 1e-7, 'betx': 1e-6, 'bety': 1e-6},
    solve=False,
    ele_start=f's.ds.l8.b1', ele_stop=f'e.ds.r8.b1',
    # Left boundary
    twiss_init='preserve_start', table_for_twiss_init=tab_boundary_left,
    targets=[
        xt.Target('alfx', 0, at='ip8'),
        xt.Target('alfy', 0, at='ip8'),
        xt.Target('betx', 1.5, at='ip8'),
        xt.Target('bety', 1.5, at='ip8'),
        xt.Target('dx', 0, at='ip8'),
        xt.Target('dpx', 0, at='ip8'),
        xt.TargetList(('betx', 'bety', 'alfx', 'alfy', 'dx', 'dpx'),
                value=tab_boundary_right, at=f'e.ds.r8.b1',
                tag='stage2'),
        xt.TargetRelPhaseAdvance('mux', mux_b1_target),
        xt.TargetRelPhaseAdvance('muy', muy_b1_target),
    ],
    vary=[
        xt.VaryList([
            f'kq6.l8b1', f'kq7.l8b1',
            f'kq8.l8b1', f'kq9.l8b1', f'kq10.l8b1', f'kqtl11.l8b1',
            f'kqt12.l8b1', f'kqt13.l8b1']
        ),
        xt.VaryList([f'kq4.l8b1', f'kq5.l8b1'], tag='stage1'),
        xt.VaryList([
            f'kq4.r8b1', f'kq5.r8b1', f'kq6.r8b1', f'kq7.r8b1',
            f'kq8.r8b1', f'kq9.r8b1', f'kq10.r8b1', f'kqtl11.r8b1',
            f'kqt12.r8b1', f'kqt13.r8b1'],
            tag='stage2')
    ]
)

opt.disable_targets(tag=['stage1', 'stage2'])
opt.disable_vary(tag=['stage2'])
opt.solve()

opt.enable_vary(tag='stage1')
opt.solve()

opt.enable_targets(tag='stage2')
opt.enable_vary(tag='stage2')
opt.solve()