import xtrack as xt
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
from scipy.spatial.transform import Rotation

env = xt.Environment()
env.call('lhc_seq.py')

env.lhcb1.particle_ref = xt.Particles(p0c=6.8e12)
env.lhcb2.particle_ref = xt.Particles(p0c=6.8e12)

env.vars.load_madx('../../test_data/lhc_2024/ats_30cm.madx')

env['on_sep5'] = 0

tw1 = env.lhcb1.twiss4d()
tw2 = env.lhcb2.twiss4d(reverse=True)

sv1 = env.lhcb1.survey()
sv2 = env.lhcb2.survey().reverse()

# Convenience function to compute aperture size and beam sizes
# ============================================================

def get_aperture_size(el):
    if hasattr(el, 'min_x'):
        return el.min_x, el.max_x
    if hasattr(el, 'max_x'):
        return -el.max_x, el.max_x
    return -el.a, el.a


def compute_beam_size(survey, twiss):
    sx = survey.X
    sy = survey.Y
    sz = survey.Z
    theta = survey.theta
    s = twiss.s
    x = twiss.x
    y = twiss.y
    bx = twiss.betx
    by = twiss.bety
    dx = twiss.dx
    dy = twiss.dy
    nemitt_x = 2.5e-6
    nemitt_y = 2.5e-6
    gamma0 = twiss.gamma0
    n_sigmas = 3 # 13.
    sigma_delta = 8e-4

    sigx = n_sigmas * np.sqrt(nemitt_x / gamma0 * bx) + abs(dx) * sigma_delta
    sigy = n_sigmas * np.sqrt(nemitt_y / gamma0 * by) + abs(dy) * sigma_delta

    return s, x, sigx, y, sigy, sx, sy, sz, theta

def ellipse(rxy, rz, beam_xy, beam_z, x, y, z, theta):
    """Make a 3D ellipse.

    Make a 3D ellipse centred at ``(x, y, z)``, with radii ``rx`` and ``rz``, and
    rotated around z-axis by the angle ``theta``. The axes are the traditional
    (matplotlib) axes.

    Parameters
    ----------
    rxy : float
        Radius in the xy-plane.
    rz : float
        z-axis radius.
    beam_xy : float
        Horizontal displacement of the centre before rotation, i.e. along theta.
    beam_z : float
        Vertical displacement of the centre before rotation.
    x : float
        Centre of the ellipse in x.
    y : float
        Centre of the ellipse in y.
    z : float
        Centre of the ellipse in z.
    theta : float
        Angle of rotation around the z-axis.
    """
    ts = np.linspace(0, 2 * np.pi, 20)
    points_xz = np.array([
        (rxy * np.cos(t) + beam_xy, 0, rz * np.sin(t) + beam_z) for t in ts]
    )
    points_xz = Rotation.from_euler('z', theta).apply(points_xz)
    return points_xz + np.tile([x, y, z], (len(ts), 1))


def mesh_from_ellipses(pts):
    num_ellipses, points_per_ellipse, _ = pts.shape
    vertices = pts.reshape(-1, 3)
    num_faces = points_per_ellipse * (num_ellipses - 1) - 1
    faces = np.hstack([
        [4, i, i + 1, points_per_ellipse + i + 1, points_per_ellipse + i]
        for i in range(num_faces)
    ])
    surface = pv.PolyData(vertices, faces)
    return surface


def plot_beam_size(ax, twiss, survey, color):
    element_around = 'ip5'
    s_around = twiss.rows[element_around].s[0]
    section_length = 130
    s_start, s_end = s_around - section_length / 2, s_around + section_length / 2

    sv = survey.rows[s_start:s_end:'s']
    tw = twiss.rows[s_start:s_end:'s']

    s, x, sigx, y, sigy, sx, sy, sz, theta = compute_beam_size(sv, tw)

    pts = np.array([
        ellipse(sigx[i], sigy[i], x[i], y[i], sx[i], sz[i], sy[i], theta[i])
        for i in range(len(sigx))
    ])

    #ax.plot_surface(pts[:, :, 0], pts[:, :, 1], pts[:, :, 2], color=color, alpha=0.5)
    surface = mesh_from_ellipses(pts)
    ax.add_mesh(surface, color=color, opacity=0.5, show_edges=True)

    # centre of the beam
    # ax.plot(
    #     sx + np.cos(theta) * x,
    #     sz + np.sin(theta) * x,
    #     sy + y,
    #     color=color,
    # )

# plt.close('all')
#
# ax = plt.figure().add_subplot(projection='3d')
ax = pv.Plotter()

ax.add_axes(
    line_width=5,
    cone_radius=0.6,
    shaft_length=0.7,
    tip_length=0.3,
    ambient=0.5,
    label_size=(0.4, 0.16),
    xlabel='X',
    ylabel='Z',
    zlabel='Y',
)

ax.set_scale(xscale=1e3, zscale=1e3)
ax.show_bounds(
    show_xaxis=False,
    show_yaxis=True,
    show_zaxis=False,
    show_xlabels=False,
    show_ylabels=True,
    show_zlabels=False,
    ytitle='Z [m]',
    location='origin',
)
title = ax.add_title(f'LHC Beam Envelopes at CMS (3σ, β*=30cm)')
x = title.GetTextProperty()
x.SetFontFamily(4)
x.SetFontFile('/Users/szymonlopaciuk/Library/Fonts/DejaVuSans.ttf')
# ax.set_box_aspect(aspect=(1, 5, 1))
#
# ax.set_title(f'LHC Beam Envelopes at CMS ($\sigma$ = 3, $\\beta^*$ = 30 cm, on_sep5 = {env["on_sep5"]:.1f})')
# ax.set_xlabel('X [m]')
# ax.set_ylabel('Z [m]')
# ax.set_zlabel('Y [m]')
#
# # ax.set_xlim(-6e3, 6e3)
# ax.set_zlim(-6e-3, 6e-3)
#
plot_beam_size(ax, tw1, sv1, color='b')
plot_beam_size(ax, tw2, sv2, color='r')
#
# plt.show()
ax.show()
