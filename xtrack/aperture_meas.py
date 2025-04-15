import numpy as np
import xtrack as xt

def measure_aperture_extent(
        line,
        dx=1e-3,
        dy=1e-3,
        x_range=(-0.1, 0.1),
        y_range=(-0.1, 0.1),
):
    x_test = np.arange(x_range[0], x_range[1], dx)
    y_test = np.arange(y_range[0], y_range[1], dy)

    n_x = len(x_test)

    x_probe = np.concatenate([x_test, 0*y_test])
    y_probe = np.concatenate([0*x_test, y_test])

    p = line.build_particles(x=x_probe, y=y_probe)

    with xt.line._preserve_config(line):
        line.freeze_longitudinal()
        line.freeze_vars(['x', 'px', 'y', 'py'])
        line.config.XSUITE_RESTORE_LOSS = True

        line.track(p, turn_by_turn_monitor='ONE_TURN_EBE')
        mon = line.record_last_track

    x_h_aper = mon.x[:n_x, :]
    s_h_aper = mon.s[:n_x, :]
    state_h_aper = mon.state[:n_x, :]
    state_h_aper[:, :-1] = state_h_aper[:, 1:] # due to the way they are logged

    mean_x = 0.5*(x_h_aper[:-1, :] + x_h_aper[1:, :])
    diff_loss_h = np.diff(state_h_aper, axis=0)
    zeros = mean_x * 0
    x_aper_low_mat = np.where(diff_loss_h>0, mean_x, zeros)
    x_aper_low_discrete = x_aper_low_mat.sum(axis=0)
    x_aper_high_mat = np.where(diff_loss_h<0, mean_x, zeros)
    x_aper_high_discrete = x_aper_high_mat.sum(axis=0)

    y_v_aper = mon.y[n_x:, :]
    state_v_aper = mon.state[n_x:, :]
    state_v_aper[:, :-1] = state_v_aper[:, 1:] # due to the way they are logged

    mean_y = 0.5*(y_v_aper[:-1, :] + y_v_aper[1:, :])
    diff_loss_v = np.diff(state_v_aper, axis=0)
    zeros = mean_y * 0
    y_aper_low_mat = np.where(diff_loss_v>0, mean_y, zeros)
    y_aper_low_discrete = y_aper_low_mat.sum(axis=0)
    y_aper_high_mat = np.where(diff_loss_v<0, mean_y, zeros)
    y_aper_high_discrete = y_aper_high_mat.sum(axis=0)

    s_aper = s_h_aper[0, :]

    mask_interp_low_h = x_aper_low_discrete != 0
    x_aper_low = np.interp(s_aper,
                            s_aper[mask_interp_low_h], x_aper_low_discrete[mask_interp_low_h])
    mask_interp_high_h = x_aper_high_discrete != 0
    x_aper_high = np.interp(s_aper,
                            s_aper[mask_interp_high_h], x_aper_high_discrete[mask_interp_high_h])
    x_aper_low_discrete[~mask_interp_low_h] = np.nan
    x_aper_high_discrete[~mask_interp_high_h] = np.nan

    mask_interp_low_v = y_aper_low_discrete != 0
    y_aper_low = np.interp(s_aper,
                            s_aper[mask_interp_low_v], y_aper_low_discrete[mask_interp_low_v])
    mask_interp_high_v = y_aper_high_discrete != 0
    y_aper_high = np.interp(s_aper,
                            s_aper[mask_interp_high_v], y_aper_high_discrete[mask_interp_high_v])
    y_aper_low_discrete[~mask_interp_low_v] = np.nan
    y_aper_high_discrete[~mask_interp_high_v] = np.nan

    # I force the values, for the case in which there are multiple apertures
    # at the same location
    x_aper_low[mask_interp_low_h] = x_aper_low_discrete[mask_interp_low_h]
    x_aper_high[mask_interp_low_h] = x_aper_high_discrete[mask_interp_low_h]
    y_aper_low[mask_interp_high_h] = y_aper_low_discrete[mask_interp_high_h]
    y_aper_high[mask_interp_high_h] = y_aper_high_discrete[mask_interp_high_h]

    # Force nan at end_point
    x_aper_low_discrete[-1] = np.nan
    x_aper_high_discrete[-1] = np.nan
    y_aper_low_discrete[-1] = np.nan
    y_aper_high_discrete[-1] = np.nan

    out = xt.Table({
        'name': np.array(list(line.element_names) + ['_end_point']),
        's': s_aper,
        'x_aper_low': x_aper_low,
        'x_aper_high': x_aper_high,
        'x_aper_low_discrete': x_aper_low_discrete,
        'x_aper_high_discrete': x_aper_high_discrete,
        'y_aper_low': y_aper_low,
        'y_aper_high': y_aper_high,
        'y_aper_low_discrete': y_aper_low_discrete,
        'y_aper_high_discrete': y_aper_high_discrete,
    })

    return out


def measure_aperture_shape(line,
        dr=1e-3,
        d_angle=np.pi / 20,
        r_max=0.1,
):
    r_range = np.linspace(0, r_max, int(r_max / dr))
    theta_range = np.linspace(0, 2*np.pi, int(2 * np.pi / d_angle))
    r, theta = np.meshgrid(r_range, theta_range)

    x_probe = r * np.cos(theta)
    y_probe = r * np.sin(theta)

    polar_shape = x_probe.shape
    x_probe = x_probe.flatten()
    y_probe = y_probe.flatten()

    p = line.build_particles(x=x_probe, y=y_probe)

    with xt.line._preserve_config(line):
        line.freeze_longitudinal()
        line.freeze_vars(['x', 'px', 'y', 'py'])
        line.config.XSUITE_RESTORE_LOSS = True

        line.track(p, turn_by_turn_monitor='ONE_TURN_EBE')
        mon = line.record_last_track

    _, num_elements = mon.s.shape

    # These will have shape (num_elements, num_thetas, num_rs), so e.g.
    # s[3, num_thetas/4] will correspond to the 'spoke' of points going outwards
    # from (0, 0) at 90° at element #3
    s = mon.s[0, :]
    x = mon.x.T.reshape(-1, *polar_shape)
    y = mon.y.T.reshape(-1, *polar_shape)
    state = mon.state.T.reshape(-1, *polar_shape)

    state_diff_loc = np.diff(state, prepend=1) < 0
    # Min assumes that there will be only one state 'transition' per 'spoke'
    x_lost = np.min(np.where(state_diff_loc, x, np.inf), axis=-1)
    y_lost = np.min(np.where(state_diff_loc, y, np.inf), axis=-1)

    x_lost[x_lost == np.inf] = np.nan
    y_lost[y_lost == np.inf] = np.nan

    assert np.all(np.isnan(x_lost) == np.isnan(y_lost))
    nans = np.isnan(x_lost) | np.isnan(y_lost)
    nans_seq = nans[:, 0]
    x_lost_interp = np.array(x_lost)
    np.putmask(x_lost_interp, nans, _multi_interp(s[nans_seq], s[~nans_seq], x_lost[~nans_seq]))
    y_lost_interp = np.array(y_lost)
    np.putmask(y_lost_interp, nans, _multi_interp(s[nans_seq], s[~nans_seq], y_lost[~nans_seq]))

    out = xt.Table({
        'name': np.array(list(line.element_names) + ['_end_point']),
        's': s,
        'polygon_x_discrete': x_lost,
        'polygon_y_discrete': y_lost,
        'polygon_x': x_lost_interp,
        'polygon_y': y_lost_interp,
        'aperture_mask': ~np.isnan(x_lost[:, 0]),
        'all_x': mon.x.T,
        'all_y': mon.y.T,
        'all_state': mon.state.T,
    })

    return out


def _multi_interp(x, sample_x, sample_value):
    # Find indices of the intervals for x
    indices = np.searchsorted(sample_x, x, side="left").clip(1, len(sample_x) - 1)
    left_indices = indices - 1
    right_indices = indices

    # Calculate the slopes and perform the interpolation
    x_left = sample_x[left_indices]
    x_right = sample_x[right_indices]
    value_left = sample_value[left_indices]
    value_right = sample_value[right_indices]

    slopes = (value_right - value_left) / (x_right - x_left)[:, np.newaxis]
    interpolated = value_right + slopes * (x[:, np.newaxis] - x_left[:, np.newaxis])

    return interpolated
