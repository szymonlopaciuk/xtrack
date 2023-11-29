import pathlib
import pytest

import numpy as np
import xtrack as xt
import xobjects as xo

test_data_folder = pathlib.Path(
        __file__).parent.joinpath('../test_data').absolute()

configurations = [
    {
        'wiggler_on': False,
        'vertical_orbit_distortion': False,
        'tilt_machine_by_90_degrees': False,
        'check_against_tracking': False,
    },
    {
        'wiggler_on': False,
        'vertical_orbit_distortion': False,
        'tilt_machine_by_90_degrees': True,
        'check_against_tracking': False,
    },
    {
        'wiggler_on': False,
        'vertical_orbit_distortion': True,
        'tilt_machine_by_90_degrees': False,
        'check_against_tracking': True,
    },
    {
        'wiggler_on': False,
        'vertical_orbit_distortion': True,
        'tilt_machine_by_90_degrees': True,
        'check_against_tracking': False,
    },
    {
        'wiggler_on': True,
        'vertical_orbit_distortion': False,
        'tilt_machine_by_90_degrees': False,
        'check_against_tracking': False,
    },
    {
        'wiggler_on': True,
        'vertical_orbit_distortion': False,
        'tilt_machine_by_90_degrees': True,
        'check_against_tracking': True,
    },
]


@pytest.mark.parametrize('conf', configurations)
def test_eq_emitt(conf):

    test_context = xo.context_default # On GPU this is too slow to run routinely

    print('===============================')
    print(conf)
    print('===============================')

    tilt_machine_by_90_degrees = conf['tilt_machine_by_90_degrees']
    wiggler_on = conf['wiggler_on']
    vertical_orbit_distortion = conf['vertical_orbit_distortion']

    line = xt.Line.from_json(test_data_folder / 'fcc_ee/fccee_h_thin.json')
    line.build_tracker()

    print('Done building tracker')

    if wiggler_on:
        line.vars['on_wiggler_v'] = 0.4

    if vertical_orbit_distortion:
        line['mwi.e5rg..0'].ksl[0] = 2e-7

    # Make sure there is no vertical bend nor skew element
    for ee in line.elements:
        if isinstance(ee, xt.Multipole):
            assert np.all(ee.ksl[1:] == 0)

    if tilt_machine_by_90_degrees:

        tw_before_tilt = line.twiss()

        # Bring the machine to the vertical plane
        for ee in line.elements:
            if isinstance(ee, xt.Multipole):
                knl = ee.knl.copy()
                ksl = ee.ksl.copy()
                hxl = ee.hxl
                hyl = ee.hyl

                ee.hxl = -hyl
                ee.hyl = -hxl

                ee.knl[0] = -ksl[0]
                ee.ksl[0] = -knl[0]
                if len(knl) > 1:
                    ee.knl[1] = -knl[1]
                    ee.ksl[1] = 0
                if len(knl) > 2:
                    ee.knl[2] = 0
                    ee.ksl[2] = knl[2]

            if isinstance(ee, xt.DipoleEdge):
                ee._r21, ee._r43 = ee._r43, ee._r21

        tw_after_tilt = line.twiss()

        assert np.isclose(tw_after_tilt.qy, tw_before_tilt.qx, rtol=0, atol=1e-8)
        assert np.isclose(tw_after_tilt.qx, tw_before_tilt.qy, rtol=0, atol=1e-8)
        assert np.isclose(tw_after_tilt.dqy, tw_before_tilt.dqx, rtol=0, atol=1e-4)
        assert np.isclose(tw_after_tilt.dqx, tw_before_tilt.dqy, rtol=0, atol=1e-4)

        assert np.allclose(tw_after_tilt.bety, tw_before_tilt.betx, rtol=3e-5, atol=0)
        assert np.allclose(tw_after_tilt.betx, tw_before_tilt.bety, rtol=3e-5, atol=0)

        assert np.allclose(tw_after_tilt.y, tw_before_tilt.x, rtol=0, atol=1e-9)
        assert np.allclose(tw_after_tilt.x, tw_before_tilt.y, rtol=0, atol=1e-9)

        assert np.allclose(tw_after_tilt.dy, tw_before_tilt.dx, rtol=0, atol=5e-6)
        assert np.allclose(tw_after_tilt.dx, tw_before_tilt.dy, rtol=0, atol=5e-6)

    line.configure_radiation(model='mean')
    line.compensate_radiation_energy_loss()

    tw_rad = line.twiss(eneloss_and_damping=True)
    ex = tw_rad.eq_gemitt_x
    ey = tw_rad.eq_gemitt_y
    ez = tw_rad.eq_gemitt_zeta

    # for regression testing
    checked = False
    if not tilt_machine_by_90_degrees and not vertical_orbit_distortion and not wiggler_on:
        assert np.isclose(ex, 7.0592e-10, atol=0,     rtol=1e-4)
        assert np.isclose(ey, 0,          atol=1e-14, rtol=0)
        assert np.isclose(ez, 3.6000e-6,  atol=0,     rtol=1e-4)
        checked = True
    elif tilt_machine_by_90_degrees and not vertical_orbit_distortion and not wiggler_on:
        assert np.isclose(ex, 0,          atol=1e-14, rtol=0)
        assert np.isclose(ey, 7.0592e-10, atol=0,     rtol=1e-4)
        assert np.isclose(ez, 3.6000e-6,  atol=0,     rtol=1e-4)
        checked = True
    elif not tilt_machine_by_90_degrees and not vertical_orbit_distortion and wiggler_on:
        assert np.isclose(ex, 6.9954e-10, atol=0,     rtol=1e-4)
        assert np.isclose(ey, 5.8575e-13, atol=0,     rtol=2e-3)
        assert np.isclose(ez, 3.8595e-6,  atol=0,     rtol=1e-4)
        checked = True
    elif tilt_machine_by_90_degrees and not vertical_orbit_distortion and wiggler_on:
        assert np.isclose(ex, 5.8575e-13, atol=0,     rtol=4e-3)  # Quite large, to be kept in mind
        assert np.isclose(ey, 6.9955e-10, atol=0,     rtol=1e-4)
        assert np.isclose(ez, 3.8595e-6,  atol=0,     rtol=1e-4)
        checked = True
    elif not tilt_machine_by_90_degrees and vertical_orbit_distortion and not wiggler_on:
        assert np.isclose(ex, 7.0576e-10, atol=0,     rtol=1e-4)
        assert np.isclose(ey, 2.5281e-12, atol=0,     rtol=2e-3)
        assert np.isclose(ez, 3.5762e-6,  atol=0,     rtol=1e-4)
        checked = True
    elif tilt_machine_by_90_degrees and vertical_orbit_distortion and not wiggler_on:
        assert np.isclose(ex, 2.5531e-12, atol=0,     rtol=4e-3)
        assert np.isclose(ey, 7.0576e-10, atol=0,     rtol=1e-4)
        assert np.isclose(ez, 3.5763e-6,  atol=0,     rtol=1e-4)
        checked = True
    else:
        raise ValueError('Unknown configuration')

    assert checked

    tw_rad2 = line.twiss(eneloss_and_damping=True, method='6d',
                     radiation_method='full',
                     compute_lattice_functions=False,
                     compute_chromatic_properties=False)

    assert 'x' in tw_rad2
    assert 'betx' not in tw_rad2
    assert 'circumference' in tw_rad2
    assert 'qx' not in tw_rad2
    assert 'dqx' not in tw_rad2

    if not vertical_orbit_distortion: # Known inconsistency to be investigated
        assert np.isclose(tw_rad2.eq_gemitt_x, tw_rad.eq_gemitt_x, atol=1e-14, rtol=1.5e-2)
        assert np.isclose(tw_rad2.eq_gemitt_y, tw_rad.eq_gemitt_y, atol=1e-14, rtol=1.5e-2)
        assert np.isclose(tw_rad2.eq_gemitt_zeta, tw_rad.eq_gemitt_zeta, atol=1e-14, rtol=4e-2)
        assert np.isclose(tw_rad2.eq_nemitt_x, tw_rad.eq_nemitt_x, atol=1e-16, rtol=1.5e-2)
        assert np.isclose(tw_rad2.eq_nemitt_y, tw_rad.eq_nemitt_y, atol=1e-16, rtol=1.5e-2)
        assert np.isclose(tw_rad2.eq_nemitt_zeta, tw_rad.eq_nemitt_zeta, atol=1e-16, rtol=4e-2)

    assert np.isclose(tw_rad.eq_nemitt_x, tw_rad.eq_gemitt_x/(tw_rad.gamma0*tw_rad.beta0), atol=1e-16, rtol=0)
    assert np.isclose(tw_rad.eq_nemitt_y, tw_rad.eq_gemitt_y/(tw_rad.gamma0*tw_rad.beta0), atol=1e-16, rtol=0)
    assert np.isclose(tw_rad.eq_nemitt_zeta, tw_rad.eq_gemitt_zeta/(tw_rad.gamma0*tw_rad.beta0), atol=1e-16, rtol=0)
    assert np.isclose(tw_rad2.eq_nemitt_x, tw_rad2.eq_gemitt_x/(tw_rad2.gamma0*tw_rad2.beta0), atol=1e-16, rtol=0)
    assert np.isclose(tw_rad2.eq_nemitt_y, tw_rad2.eq_gemitt_y/(tw_rad2.gamma0*tw_rad2.beta0), atol=1e-16, rtol=0)
    assert np.isclose(tw_rad2.eq_nemitt_zeta, tw_rad2.eq_gemitt_zeta/(tw_rad2.gamma0*tw_rad2.beta0), atol=1e-16, rtol=0)

    if conf['check_against_tracking']:

        line.discard_tracker()
        line.configure_radiation(model='quantum')
        line.build_tracker(_context=test_context)

        p = line.build_particles(num_particles=30)
        line.track(p, num_turns=400, turn_by_turn_monitor=True, time=True)
        mon = line.record_last_track
        print(f'Tracking time: {line.time_last_track}')

        sigma_x_eq = float(np.sqrt(ex * tw_rad.betx[0] + ey * tw_rad.betx2[0] + (np.std(p.delta) * tw_rad.dx[0])**2))
        sigma_y_eq = float(np.sqrt(ex * tw_rad.bety1[0] + ey * tw_rad.bety[0] + (np.std(p.delta) * tw_rad.dy[0])**2))
        sigma_zeta_eq = float(np.sqrt(ez * tw_rad.betz0))

        sigma_x_track = np.std(mon.x, axis=0)[-200:]
        sigma_y_track = np.std(mon.y, axis=0)[-200:]
        sigma_zeta_track = np.std(mon.zeta, axis=0)[-200:]

        if sigma_x_eq > 1e-8:
            assert np.min(np.abs(sigma_x_track/sigma_x_eq - 1.)) < 0.1
        if sigma_y_eq > 1e-8:
            assert np.min(np.abs(sigma_y_track/sigma_y_eq - 1.)) < 0.1
        assert np.min(np.abs(sigma_zeta_track/sigma_zeta_eq - 1.)) < 0.1

        assert np.isclose(sigma_x_eq, np.mean(sigma_x_track), rtol=0.3, atol=1e-9)
        assert np.isclose(sigma_y_eq, np.mean(sigma_y_track), rtol=0.3, atol=1e-9)
        assert np.isclose(sigma_zeta_eq, np.mean(sigma_zeta_track), rtol=0.3, atol=1e-9)
