# copyright ############################### #
# This file is part of the Xtrack Package.  #
# Copyright (c) CERN, 2021.                 #
# ######################################### #
import json
import pathlib
import pytest

import numpy as np
import xobjects as xo
import xtrack as xt
import xpart as xp
from xobjects.test_helpers import for_all_test_contexts

test_data_folder = pathlib.Path(
    __file__).parent.joinpath('../test_data').absolute()


@for_all_test_contexts
def test_simple_collective_line(test_context):
    num_turns = 100
    elements = [xt.Drift(length=2., _context=test_context) for _ in range(5)]
    elements[3].iscollective = True
    line = xt.Line(elements=elements)
    line.reset_s_at_end_turn = False

    particles = xp.Particles(x=[1e-3, 2e-3, 3e-3], p0c=7e12,
            _context=test_context)
    line.build_tracker(_context=test_context)
    line.track(particles, num_turns=num_turns)

    particles.move(_context=xo.ContextCpu())

    assert np.all(particles.at_turn == num_turns)
    assert np.allclose(particles.s, 10 * num_turns, rtol=0, atol=1e-14)



@for_all_test_contexts
def test_ebe_monitor(test_context):
    line = xt.Line(elements=[xt.Multipole(knl=[0, 1.]),
                            xt.Drift(length=0.5),
                            xt.Multipole(knl=[0, -1]),
                            xt.Cavity(frequency=400e7, voltage=6e6),
                            xt.Drift(length=.5),
                            xt.Drift(length=0)])

    line.build_tracker(_context=test_context)

    particles = xp.Particles(x=[1e-3, -2e-3, 5e-3], y=[2e-3, -4e-3, 3e-3],
                            zeta=1e-2, p0c=7e12, mass0=xp.PROTON_MASS_EV,
                            _context=test_context)

    line.track(particles.copy(), turn_by_turn_monitor='ONE_TURN_EBE')

    mon = line.record_last_track

    for ii, ee in enumerate(line.elements):
        for tt, nn in particles.per_particle_vars:
            assert np.all(particles.to_dict()[nn] == getattr(mon, nn)[:, ii])
        ee.track(particles)
        particles.at_element += 1


@for_all_test_contexts
def test_cycle(test_context):
    d0 = xt.Drift()
    c0 = xt.Cavity()
    d1 = xt.Drift()
    r0 = xt.SRotation()
    particle_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, gamma0=1.05)

    for collective in [True, False]:
        line = xt.Line(elements=[d0, c0, d1, r0])
        d1.iscollective = collective

        line.build_tracker(_context=test_context)
        line.particle_ref = particle_ref

        cline_name = line.cycle(name_first_element='e2')
        cline_index = line.cycle(index_first_element=2)

        assert cline_name.tracker is not None
        assert cline_index.tracker is not None

        for cline in [cline_index, cline_name]:
            assert cline.element_names[0] == 'e2'
            assert cline.element_names[1] == 'e3'
            assert cline.element_names[2] == 'e0'
            assert cline.element_names[3] == 'e1'

            assert cline.elements[0] is d1
            assert cline.elements[1] is r0
            assert cline.elements[2] is d0
            assert cline.elements[3] is c0

            assert cline.particle_ref.mass0 == xp.PROTON_MASS_EV
            assert cline.particle_ref.gamma0 == 1.05


@for_all_test_contexts
def test_synrad_configuration(test_context):
    for collective in [False, True]:
        elements = [xt.Multipole(knl=[1]) for _ in range(10)]
        if collective:
            elements[5].iscollective = True
            elements[5].move(_context=test_context)

        line = xt.Line(elements=elements)
        line.configure_radiation(model='mean')
        line.build_tracker(_context=test_context)

        for ee in line.elements:
            assert ee.radiation_flag == 1
        p = xp.Particles(x=[0.01, 0.02], _context=test_context)
        line.track(p)
        p.move(_context=xo.ContextCpu())
        assert np.all(p._rng_s1 + p._rng_s2 + p._rng_s3 + p._rng_s4 == 0)

        line.configure_radiation(model='quantum')
        for ee in line.elements:
            assert ee.radiation_flag == 2
        p = xp.Particles(x=[0.01, 0.02], _context=test_context)
        line.track(p)
        p.move(_context=xo.ContextCpu())
        assert np.all(p._rng_s1 + p._rng_s2 + p._rng_s3 + p._rng_s4 > 0)

        line.configure_radiation(model=None)
        for ee in line.elements:
            assert ee.radiation_flag == 0
        p = xp.Particles(x=[0.01, 0.02], _context=test_context)
        line.track(p)
        p.move(_context=xo.ContextCpu())
        assert np.all(p._rng_s1 + p._rng_s2 + p._rng_s3 + p._rng_s4 > 0)


@for_all_test_contexts
def test_partial_tracking(test_context):
    n_elem = 9
    elements = [ xt.Drift(length=1.) for _ in range(n_elem) ]
    line = xt.Line(elements=elements)
    line.build_tracker(_context=test_context)
    assert not line.iscollective
    particles_init = xp.Particles(_context=test_context,
        x=[1e-3, -2e-3, 5e-3], y=[2e-3, -4e-3, 3e-3],
        zeta=1e-2, p0c=7e12, mass0=xp.PROTON_MASS_EV,
        at_turn=0, at_element=0)

    _default_track(line, particles_init)
    _ele_start_until_end(line, particles_init)
    _ele_start_with_shift(line, particles_init)
    _ele_start_with_shift_more_turns(line, particles_init)
    _ele_stop_from_start(line, particles_init)
    _ele_start_to_ele_stop(line, particles_init)
    _ele_start_to_ele_stop_with_overflow(line, particles_init)


@for_all_test_contexts
def test_partial_tracking_with_collective(test_context):
    n_elem = 9
    elements = [xt.Drift(length=1., _context=test_context) for _ in range(n_elem)]
    # Make some elements collective
    elements[3].iscollective = True
    elements[7].iscollective = True
    line = xt.Line(elements=elements)
    line.build_tracker(_context=test_context)
    assert line.iscollective
    assert len(line.tracker._parts) == 5
    particles_init = xp.Particles(
            _context=test_context,
            x=[1e-3, -2e-3, 5e-3], y=[2e-3, -4e-3, 3e-3],
            zeta=1e-2, p0c=7e12, mass0=xp.PROTON_MASS_EV,
            at_turn=0, at_element=0)

    _default_track(line, particles_init)
    _ele_start_until_end(line, particles_init)
    _ele_start_with_shift(line, particles_init)
    _ele_start_with_shift_more_turns(line, particles_init)
    _ele_stop_from_start(line, particles_init)
    _ele_start_to_ele_stop(line, particles_init)
    _ele_start_to_ele_stop_with_overflow(line, particles_init)


# Track from the start until the end of the first, second, and tenth turn
def _default_track(line, particles_init):
    n_elem = len(line.element_names)
    for turns in [1, 2, 10]:
        expected_end_turn = turns
        expected_end_element = 0
        expected_num_monitor = expected_end_turn if expected_end_element==0 else expected_end_turn+1

        particles = particles_init.copy()
        line.track(particles, num_turns=turns, turn_by_turn_monitor=True)
        check, end_turn, end_element, end_s = _get_at_turn_element(particles)
        assert (check and end_turn==expected_end_turn and end_element==expected_end_element
                    and end_s==expected_end_element)
        assert line.record_last_track.x.shape == (len(particles.x), expected_num_monitor)


# Track, from any ele_start, until the end of the first, second, and tenth turn
def _ele_start_until_end(line, particles_init):
    n_elem = len(line.element_names)
    for turns in [1, 2, 10]:
        for start in range(n_elem):
            expected_end_turn = turns
            expected_end_element = 0
            expected_num_monitor = expected_end_turn if expected_end_element==0 else expected_end_turn+1

            particles = particles_init.copy()
            particles.at_element = start
            particles.s = start
            line.track(particles, num_turns=turns, ele_start=start, turn_by_turn_monitor=True)
            check, end_turn, end_element, end_s = _get_at_turn_element(particles)
            assert (check and end_turn==expected_end_turn and end_element==expected_end_element
                        and end_s==expected_end_element)
            assert line.record_last_track.x.shape==(len(particles.x), expected_num_monitor)


# Track, from any ele_start, any shifts that stay within the first turn
def _ele_start_with_shift(line, particles_init):
    n_elem = len(line.element_names)
    for start in range(n_elem):
        for shift in range(1,n_elem-start):
            expected_end_turn = 0
            expected_end_element = start+shift
            expected_num_monitor = expected_end_turn if expected_end_element==0 else expected_end_turn+1

            particles = particles_init.copy()
            particles.at_element = start
            particles.s = start
            line.track(particles, ele_start=start, num_elements=shift, turn_by_turn_monitor=True)
            check, end_turn, end_element, end_s = _get_at_turn_element(particles)
            assert (check and end_turn==expected_end_turn and end_element==expected_end_element
                        and end_s==expected_end_element)
            assert line.record_last_track.x.shape==(len(particles.x), expected_num_monitor)

# Track, from any ele_start, any shifts that are larger than one turn (up to 3 turns)
def _ele_start_with_shift_more_turns(line, particles_init):
    n_elem = len(line.element_names)
    for start in range(n_elem):
        for shift in range(n_elem-start, 3*n_elem+1):
            expected_end_turn = round(np.floor( (start+shift)/n_elem ))
            expected_end_element = start + shift - n_elem*expected_end_turn
            expected_num_monitor = expected_end_turn if expected_end_element==0 else expected_end_turn+1

            particles = particles_init.copy()
            particles.at_element = start
            particles.s = start
            line.track(particles, ele_start=start, num_elements=shift, turn_by_turn_monitor=True)
            check, end_turn, end_element, end_s = _get_at_turn_element(particles)
            assert (check and end_turn==expected_end_turn and end_element==expected_end_element
                        and end_s==expected_end_element)
            assert line.record_last_track.x.shape==(len(particles.x), expected_num_monitor)


# Track from the start until any ele_stop in the first, second, and tenth turn
def _ele_stop_from_start(line, particles_init):
    n_elem = len(line.element_names)
    for turns in [1, 2, 10]:
        for stop in range(1, n_elem):
            expected_end_turn = turns-1
            expected_end_element = stop
            expected_num_monitor = expected_end_turn if expected_end_element==0 else expected_end_turn+1

            particles = particles_init.copy()
            line.track(particles, num_turns=turns, ele_stop=stop, turn_by_turn_monitor=True)
            check, end_turn, end_element, end_s = _get_at_turn_element(particles)
            assert (check and end_turn==expected_end_turn and end_element==expected_end_element
                        and end_s==expected_end_element)
            assert line.record_last_track.x.shape==(len(particles.x), expected_num_monitor)


# Track from any ele_start until any ele_stop that is larger than ele_start
# for one, two, and ten turns
def _ele_start_to_ele_stop(line, particles_init):
    n_elem = len(line.element_names)
    for turns in [1, 2, 10]:
        for start in range(n_elem):
            for stop in range(start+1,n_elem):
                expected_end_turn = turns-1
                expected_end_element = stop
                expected_num_monitor = expected_end_turn if expected_end_element==0 else expected_end_turn+1

                particles = particles_init.copy()
                particles.at_element = start
                particles.s = start
                line.track(particles, num_turns=turns, ele_start=start, ele_stop=stop, turn_by_turn_monitor=True)
                check, end_turn, end_element, end_s = _get_at_turn_element(particles)
                assert (check and end_turn==expected_end_turn and end_element==expected_end_element
                            and end_s==expected_end_element)
                assert line.record_last_track.x.shape==(len(particles.x), expected_num_monitor)


# Track from any ele_start until any ele_stop that is smaller than or equal to ele_start (turn increses by one)
# for one, two, and ten turns
def _ele_start_to_ele_stop_with_overflow(line, particles_init):
    n_elem = len(line.element_names)
    for turns in [1, 2, 10]:
        for start in range(n_elem):
            for stop in range(start+1):
                expected_end_turn = turns
                expected_end_element = stop
                expected_num_monitor = expected_end_turn if expected_end_element==0 else expected_end_turn+1

                particles = particles_init.copy()
                particles.at_element = start
                particles.s = start
                line.track(particles, num_turns=turns, ele_start=start, ele_stop=stop, turn_by_turn_monitor=True)
                check, end_turn, end_element, end_s = _get_at_turn_element(particles)
                assert (check and end_turn==expected_end_turn and end_element==expected_end_element
                            and end_s==expected_end_element)
                assert line.record_last_track.x.shape==(len(particles.x), expected_num_monitor)


# Quick helper function to:
#   1) check that all survived particles are at the same element and turn
#   2) return that element and turn
def _get_at_turn_element(particles):
    part_cpu = particles.copy(_context=xo.ContextCpu())
    at_element = np.unique(part_cpu.at_element[part_cpu.state>0])
    at_turn = np.unique(part_cpu.at_turn[part_cpu.state>0])
    at_s = np.unique(part_cpu.s[part_cpu.state>0])
    all_together = len(at_turn)==1 and len(at_element)==1 and len(at_s)==1
    return all_together, at_turn[0], at_element[0], at_s[0]

def test_tracker_hashable_config():
    line = xt.Line([])
    line.build_tracker()
    line.config.TEST_FLAG_BOOL = True
    line.config.TEST_FLAG_INT = 42
    line.config.TEST_FLAG_FALSE = False
    line.config.ZZZ = 'lorem'
    line.config.AAA = 'ipsum'

    expected = (
        ('AAA', 'ipsum'),
        ('TEST_FLAG_BOOL', True),
        ('TEST_FLAG_INT', 42),
        ('XFIELDS_BB3D_NO_BEAMSTR', True), # active by default
        ('XFIELDS_BB3D_NO_BHABHA', True), # active by default
        ('XTRACK_GLOBAL_XY_LIMIT', 1.0), # active by default
        ('XTRACK_MULTIPOLE_NO_SYNRAD', True), # active by default
        ('ZZZ', 'lorem'),
    )

    assert line.tracker._hashable_config() == expected


def test_tracker_config_to_headers():
    line = xt.Line([])
    line.build_tracker()

    line.config.clear()
    line.config.TEST_FLAG_BOOL = True
    line.config.TEST_FLAG_INT = 42
    line.config.TEST_FLAG_FALSE = False
    line.config.ZZZ = 'lorem'
    line.config.AAA = 'ipsum'

    expected = [
        '#define TEST_FLAG_BOOL',
        '#define TEST_FLAG_INT 42',
        '#define ZZZ lorem',
        '#define AAA ipsum',
    ]

    assert set(line.tracker._config_to_headers()) == set(expected)


@for_all_test_contexts
def test_tracker_config(test_context):
    class TestElement(xt.BeamElement):
        _xofields = {
            'dummy': xo.Float64,
        }
        _extra_c_sources = ["""
            /*gpufun*/
            void TestElement_track_local_particle(
                    TestElementData el,
                    LocalParticle* part0)
            {
                //start_per_particle_block (part0->part)

                    #if TEST_FLAG == 2
                    LocalParticle_set_x(part, 7);
                    #endif

                    #ifdef TEST_FLAG_BOOL
                    LocalParticle_set_y(part, 42);
                    #endif

                //end_per_particle_block
            }
            """]

    test_element = TestElement(_context=test_context)
    line = xt.Line([test_element])
    line.build_tracker(_context=test_context)

    particles = xp.Particles(p0c=1e9, x=[0], y=[0], _context=test_context)

    p = particles.copy()
    line.config.TEST_FLAG = 2
    line.track(p)
    assert p.x[0] == 7.0
    assert p.y[0] == 0.0
    first_kernel, first_data = line.tracker.get_track_kernel_and_data_for_present_config()

    p = particles.copy()
    line.config.TEST_FLAG = False
    line.config.TEST_FLAG_BOOL = True
    line.track(p)
    assert p.x[0] == 0.0
    assert p.y[0] == 42.0
    current_kernel, current_data = line.tracker.get_track_kernel_and_data_for_present_config()
    assert current_kernel is not first_kernel
    assert current_data is not first_data

    line.config.TEST_FLAG = 2
    line.config.TEST_FLAG_BOOL = False
    assert len(line.tracker.track_kernel) == 3 # As line.track_kernel.keys() =
                                          # dict_keys([(), (('TEST_FLAG', 2),), (('TEST_FLAG_BOOL', True),)])
    current_kernel, current_data = line.tracker.get_track_kernel_and_data_for_present_config()
    assert current_kernel is first_kernel
    assert current_data is first_data


@for_all_test_contexts
def test_optimize_for_tracking(test_context):
    fname_line_particles = test_data_folder / 'hllhc15_noerrors_nobb/line_and_particle.json'

    with open(fname_line_particles, 'r') as fid:
        input_data = json.load(fid)

    line = xt.Line.from_dict(input_data['line'])
    line.particle_ref = xp.Particles.from_dict(input_data['particle'])

    line.build_tracker(_context=test_context)

    particles = line.build_particles(
        x_norm=np.linspace(-2, 2, 1000), y_norm=0.1, delta=3e-4,
        nemitt_x=2.5e-6, nemitt_y=2.5e-6)

    p_no_optimized = particles.copy()
    p_optimized = particles.copy()

    num_turns = 10

    line.track(p_no_optimized, num_turns=num_turns, time=True)
    df_before_optimize = line.to_pandas()
    n_markers_before_optimize = (df_before_optimize.element_type == 'Marker').sum()
    assert n_markers_before_optimize > 4 # There are at least the IPs

    line.optimize_for_tracking(keep_markers=True)
    df_optimize_keep_markers = line.to_pandas()
    n_markers_optimize_keep = (df_optimize_keep_markers.element_type == 'Marker').sum()
    assert n_markers_optimize_keep == n_markers_before_optimize

    line.optimize_for_tracking(keep_markers=['ip1', 'ip5'])
    df_optimize_ip15 = line.to_pandas()
    n_markers_optimize_ip15 = (df_optimize_ip15.element_type == 'Marker').sum()
    assert n_markers_optimize_ip15 == 2

    line.optimize_for_tracking()

    assert type(line['mb.b10l3.b1..1']) is xt.SimpleThinBend
    assert type(line['mq.10l3.b1..1']) is xt.SimpleThinQuadrupole

    df_optimize = line.to_pandas()
    n_markers_optimize = (df_optimize.element_type == 'Marker').sum()
    assert n_markers_optimize == 0

    n_multipoles_before_optimize = (df_before_optimize.element_type == 'Multipole').sum()
    n_multipoles_optimize = (df_optimize.element_type == 'Multipole').sum()
    assert n_multipoles_before_optimize > n_multipoles_optimize

    n_drifts_before_optimize = (df_before_optimize.element_type == 'Drift').sum()
    n_drifts_optimize = (df_optimize.element_type == 'Drift').sum()
    assert n_drifts_before_optimize > n_drifts_optimize

    line.track(p_optimized, num_turns=num_turns, time=True)

    p_no_optimized.move(xo.context_default)
    p_optimized.move(xo.context_default)

    assert np.all(p_no_optimized.state == 1)
    assert np.all(p_optimized.state == 1)

    assert np.allclose(p_no_optimized.x, p_optimized.x, rtol=0, atol=1e-14)
    assert np.allclose(p_no_optimized.y, p_optimized.y, rtol=0, atol=1e-14)
    assert np.allclose(p_no_optimized.px, p_optimized.px, rtol=0, atol=1e-14)
    assert np.allclose(p_no_optimized.py, p_optimized.py, rtol=0, atol=1e-14)
    assert np.allclose(p_no_optimized.zeta, p_optimized.zeta, rtol=0, atol=1e-11)
    assert np.allclose(p_no_optimized.delta, p_optimized.delta, rtol=0, atol=1e-14)


@for_all_test_contexts
def test_backtrack_with_flag(test_context):

    line = xt.Line.from_json(test_data_folder /
                'hllhc15_noerrors_nobb/line_w_knobs_and_particle.json')
    line.build_tracker(_context=test_context)

    line.vars['on_crab1'] = -190
    line.vars['on_crab5'] = -190
    line.vars['on_x1'] = 130
    line.vars['on_x5'] = 130

    p = xp.Particles(_context=test_context,
        p0c=7000e9, x=1e-4, px=1e-6, y=2e-4, py=3e-6, zeta=0.01, delta=1e-4)

    line.track(p, turn_by_turn_monitor='ONE_TURN_EBE')
    mon_forward = line.record_last_track

    line.track(p, backtrack=True, turn_by_turn_monitor='ONE_TURN_EBE')
    mon_backtrack = line.record_last_track

    assert np.allclose(mon_forward.x, mon_backtrack.x, rtol=0, atol=1e-10)
    assert np.allclose(mon_forward.y, mon_backtrack.y, rtol=0, atol=1e-10)
    assert np.allclose(mon_forward.px, mon_backtrack.px, rtol=0, atol=1e-10)
    assert np.allclose(mon_forward.py, mon_backtrack.py, rtol=0, atol=1e-10)
    assert np.allclose(mon_forward.zeta, mon_backtrack.zeta, rtol=0, atol=1e-10)
    assert np.allclose(mon_forward.delta, mon_backtrack.delta, rtol=0, atol=1e-10)


@for_all_test_contexts
@pytest.mark.parametrize(
    'with_progress,turns',
    [(True, 300), (True, 317), (7, 523), (1, 21), (10, 10)]
)
@pytest.mark.parametrize('collective', [True, False], ids=['collective', 'non-collective'])
def test_tracking_with_progress(test_context, with_progress, turns, collective):
    elements = [xt.Drift(length=2, _context=test_context) for _ in range(5)]
    elements[3].iscollective = collective
    line = xt.Line(elements=elements)
    line.reset_s_at_end_turn = False

    particles = xp.Particles(x=[1e-3, 2e-3, 3e-3], p0c=7e12, _context=test_context)
    line.build_tracker(_context=test_context)
    line.track(particles, num_turns=turns, with_progress=with_progress)
    particles.move(xo.ContextCpu())

    assert np.all(particles.at_turn == turns)
    assert np.allclose(particles.s, 10 * turns, rtol=0, atol=1e-14)


@for_all_test_contexts
@pytest.mark.parametrize(
    'ele_start,ele_stop,expected_x',
    [
        (None, None, [0, 0.005, 0.010, 0.015, 0.020, 0.025]),
        (None, 3, [0, 0.005, 0.010, 0.015, 0.020, 0.023]),
        (2, None, [0, 0.003, 0.008, 0.013, 0.018, 0.023]),
        (2, 3, [0, 0.003, 0.008, 0.013, 0.018, 0.021]),
        (3, 2, [0, 0.002, 0.007, 0.012, 0.017, 0.022, 0.024]),
    ],
)
@pytest.mark.parametrize('with_progress', [False, True, 1, 2, 3])
def test_tbt_monitor_with_progress(test_context, ele_start, ele_stop, expected_x, with_progress):
    line = xt.Line(elements=[xt.Drift(length=1, _context=test_context)] * 5)
    line.build_tracker(_context=test_context)

    p = xt.Particles(px=0.001)
    line.track(p, num_turns=5, turn_by_turn_monitor=True, with_progress=with_progress, ele_start=ele_start, ele_stop=ele_stop)

    monitor_recorded_x = line.record_last_track.x
    assert monitor_recorded_x.shape == (1, len(expected_x) - 1)

    recorded_x = np.concatenate([monitor_recorded_x[0], p.x])
    assert np.allclose(recorded_x, expected_x, atol=1e-16)
