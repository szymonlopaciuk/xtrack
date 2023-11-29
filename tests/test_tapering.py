import json
import pathlib
import numpy as np
import xtrack as xt

test_data_folder = pathlib.Path(
        __file__).parent.joinpath('../test_data').absolute()

def test_tapering_and_twiss_with_radiation():

    filename = test_data_folder / 'clic_dr/line_for_taper.json'
    configs = [
        {'radiation_method': None, 'p0_correction': True, 'cavity_preserve_angle': False, 'beta_rtol': 2e-2, 'q_atol': 5e-4},
        {'radiation_method': 'kick_as_co', 'p0_correction': True, 'cavity_preserve_angle': False, 'beta_rtol': 2e-2, 'q_atol': 5e-4},
        {'radiation_method': 'kick_as_co', 'p0_correction': True, 'cavity_preserve_angle': True, 'beta_rtol': 1e-3, 'q_atol': 5e-4},
        {'radiation_method': 'scale_as_co', 'p0_correction': True, 'cavity_preserve_angle': True, 'beta_rtol': 1e-5, 'q_atol': 5e-4},
        {'radiation_method': 'full', 'p0_correction': False, 'cavity_preserve_angle': False, 'beta_rtol': 2e-2, 'q_atol': 5e-4},
        {'radiation_method': 'full', 'p0_correction': True, 'cavity_preserve_angle': False, 'beta_rtol': 2e-2, 'q_atol': 5e-4},
        {'radiation_method': 'full', 'p0_correction': True, 'cavity_preserve_angle': True, 'beta_rtol': 2e-5, 'q_atol': 5e-4},
    ]

    with open(filename, 'r') as f:
        line = xt.Line.from_dict(json.load(f))

    # Initial twiss (no radiation)
    line.configure_radiation(model=None)
    line.build_tracker()
    tw_no_rad = line.twiss(method='4d', freeze_longitudinal=True)

    assert tw_no_rad.radiation_method == None

    # Enable radiation
    line.configure_radiation(model='mean')
    # - Set cavity lags to compensate energy loss
    # - Taper magnet strengths
    line.compensate_radiation_energy_loss()

    for conf in configs:
        print(f'Running test with conf: {conf}')

        line.config.XTRACK_CAVITY_PRESERVE_ANGLE = conf['cavity_preserve_angle']

        extra_kwargs = {}
        if conf['radiation_method'] == 'kick_as_co' and conf['cavity_preserve_angle']:
            extra_kwargs['matrix_stability_tol'] = 0.1

        if conf['radiation_method'] == 'scale_as_co':
            extra_kwargs['matrix_stability_tol'] = 0.1
            extra_kwargs['use_full_inverse'] = True

        print('Twiss with radiation')
        # Twiss(es) with radiation
        tw = line.twiss(radiation_method=conf['radiation_method'],
                        eneloss_and_damping=True, **extra_kwargs)
        print('Done')

        assert tw.radiation_method == (conf['radiation_method'] or 'kick_as_co')

        # Check twiss at_s
        print('Twiss at_s')
        i_ele = len(tw.s)//3
        tws = line.twiss(at_s=tw.s[i_ele],
                        radiation_method=conf['radiation_method'],
                        eneloss_and_damping=True, **extra_kwargs)
        print('Done')

        line.config.XTRACK_CAVITY_PRESERVE_ANGLE = False

        if conf['p0_correction']:
            p0corr = 1 + tw.delta
        else:
            p0corr = 1

        # mask for taperable elemements
        tt = line.get_table().rows[:-1] # remove endpoint
        mask_taperable = (tt.element_type == 'Multipole') | (tt.element_type == 'DipoleEdge')
        assert np.sum(mask_taperable) == 17420

        delta_taper = line.attr['delta_taper']
        assert np.allclose(delta_taper[mask_taperable],
            0.5*(tw.delta[:-1] + tw.delta[1:])[mask_taperable], rtol=0, atol=1e-6)

        assert np.isclose(tw.delta[0], 0, rtol=0, atol=1e-6)
        assert np.isclose(tw.delta[-1], 0, rtol=0, atol=1e-6)

        assert np.isclose(tw.qx, tw_no_rad.qx, rtol=0, atol=conf['q_atol'])
        assert np.isclose(tw.qy, tw_no_rad.qy, rtol=0, atol=conf['q_atol'])

        assert np.isclose(tw.dqx, tw_no_rad.dqx, rtol=0, atol=1.5e-2*tw.qx)
        assert np.isclose(tw.dqy, tw_no_rad.dqy, rtol=0, atol=1.5e-2*tw.qy)

        assert np.allclose(tw.x, tw_no_rad.x, rtol=0, atol=1e-7)
        assert np.allclose(tw.y, tw_no_rad.y, rtol=0, atol=1e-7)

        assert np.allclose(tw.betx*p0corr, tw_no_rad.betx, rtol=conf['beta_rtol'], atol=0)
        assert np.allclose(tw.bety*p0corr, tw_no_rad.bety, rtol=conf['beta_rtol'], atol=0)

        assert np.allclose(tw.dx, tw.dx, rtol=0.0, atol=0.1e-3)

        assert np.allclose(tw.dy, tw.dy, rtol=0.0, atol=0.1e-3)

        assert np.allclose(tw.s[i_ele], tws.s, rtol=0, atol=1e-7)
        assert np.allclose(tw.x[i_ele], tws.x, rtol=0, atol=1e-7)
        assert np.allclose(tw.y[i_ele], tws.y, rtol=0, atol=1e-7)
        assert np.allclose(tw.betx[i_ele], tws.betx, rtol=1e-3, atol=0)
        assert np.allclose(tw.bety[i_ele], tws.bety, rtol=1e-3, atol=0)

        eneloss = tw.eneloss_turn
        assert eneloss/line.particle_ref.energy0 > 0.01
        assert np.isclose(
            line['rf'].voltage*np.sin((line['rf'].lag + line['rf'].lag_taper)/180*np.pi),
            eneloss/4, rtol=1e-5)
        assert np.isclose(
            line['rf1'].voltage*np.sin((line['rf'].lag + line['rf'].lag_taper)/180*np.pi),
            eneloss/4, rtol=1e-5)
        assert np.isclose(
            line['rf2a'].voltage*np.sin((line['rf'].lag + line['rf'].lag_taper)/180*np.pi),
            eneloss/4*0.6, rtol=1e-5)
        assert np.isclose(
            line['rf2b'].voltage*np.sin((line['rf'].lag + line['rf'].lag_taper)/180*np.pi),
            eneloss/4*0.4, rtol=1e-5)
        assert np.isclose(
            line['rf3'].voltage*np.sin((line['rf'].lag + line['rf'].lag_taper)/180*np.pi),
            eneloss/4, rtol=1e-5)
