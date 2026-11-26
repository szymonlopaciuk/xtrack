// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2025.                 //
// ######################################### //
#ifndef XTRACK_TRACK_MAGNET_EDGE_H
#define XTRACK_TRACK_MAGNET_EDGE_H

#include "xtrack/headers/track.h"
#include "xtrack/beam_elements/elements_src/track_magnet_structs.h"
#include "xtrack/beam_elements/elements_src/edge_flags.h"
#include "xtrack/beam_elements/elements_src/track_dipole_edge_linear.h"
#include "xtrack/beam_elements/elements_src/track_yrotation.h"
#include "xtrack/beam_elements/elements_src/track_wedge.h"
#include "xtrack/beam_elements/elements_src/track_mult_fringe.h"
#include "xtrack/beam_elements/elements_src/track_dipole_fringe.h"
#include "xtrack/beam_elements/elements_src/fringe_convergence.h"


// A magnet face is tracked as a fixed sequence of independent sub-maps, each
// selected by the fringe mask of the face (edge_flags.h), following MAD-NG
// (strex_fringe). Going from the drift into the field the order is
//
//     patch -> vector potential -> y-rotation -> bend fringe -> mult fringe
//         -> Lee-Whiting fringe -> quad wedge -> wedge
//
// and going from the field into the drift it is the reverse. Every sub-map
// is written as the forward map of its face, in the element frame, and
// applies its exact inverse when backtracking, so that backtracking crosses
// the faces in the reverse order and undoes them one by one.
//
// EDGE_LINEAR replaces the dipole and multipole fringes and the face geometry
// (y-rotation, wedges) with the r21/r43 linear edge; the Lee-Whiting fringe,
// which has no linear version, is not affected. A face without EDGE_ACTIVE is
// skipped altogether.
// MagnetFace_prepare turns the mask and the parameters into one do_* switch
// per sub-map (off for zero parameters), and every sub-map returns at once
// when its switch is off.


// Vector potential: picked up when entering the field region, dropped when
// leaving it (we assume a drift outside the magnet).
GPUFUN
void MagnetEdge_vector_potential_single_particle(
    LocalParticle* part,
    const uint8_t is_leaving_magnet,
    xt_float_or_tpsa_arg ksol,
    const double x0_solenoid,
    const double y0_solenoid
) {
    if (is_leaving_magnet) {
        LocalParticle_set_ax(part, 0.);
        LocalParticle_set_ay(part, 0.);
    }
    else {
        LocalParticle_set_ax(
            part, -0.5 * ksol * (LocalParticle_get_y(part) - y0_solenoid)
        );
        LocalParticle_set_ay(
            part, 0.5 * ksol * (LocalParticle_get_x(part) - x0_solenoid)
        );
    }
}


// cosh(f) and sinh(f)/f as functions of m = f^2, regular at m = 0 (series
// below m = 1 on the constant part, so that they also hold for TPSA)
GPUFUN
void MagnetEdge_cosh_sinhc(
    xt_float_or_tpsa_arg m,
    xt_float_or_tpsa* cosh_f,
    xt_float_or_tpsa* sinhc_f
) {
    if (xt_float_or_tpsa_const_part(m) < 1.) {
        xt_float_or_tpsa c = 1., s = 1.;
        for (int n = 12; n >= 1; n--) {
            c = 1. + m * c / ((2. * n - 1.) * (2. * n));
            s = 1. + m * s / ((2. * n) * (2. * n + 1.));
        }
        *cosh_f = c;
        *sinhc_f = s;
    }
    else {
        const xt_float_or_tpsa f = sqrt(m);
        *cosh_f = cosh(f);
        *sinhc_f = sinh(f) / f;
    }
}


// Lee-Whiting soft quadrupole fringe. Sources:
//   - the map: E. Forest, ch. 13.2.3, eq. 13.33;
//   - its use in production: SAD, subroutine tqlfre (bz = 0 branch):
//     https://github.com/KatsOide/SAD/blob/d082c7c5161ea111132e12c0cf9288fb892c23ba/src/tqlfre.f#L18
//   - the implementation followed here, MAD-NG qsad_fringe ([FRINGE2QUAD],
//     marked "unchecked" there):
//     https://github.com/MethodicalAcceleratorDesign/MAD/blob/d3cabd9cdebde62ebedb51bab61ac033b9159489/src/madl_dynmap.mad#L1928
//     with the sign conventions (lw, sdir) of strex_fringe (#L1594).
// MAD-NG evaluates it in the frame rotated by -atan2(k1s, k1)/2, which is
// singular at k1 = k1s = 0. Here the rotation is absorbed into the matrix
// K = [[k1, -k1s], [-k1s, -k1]] (= b2 times the reflection in that frame,
// K^2 = b2^2), so that the exponentials become
// cosh(f) + sinh(f)/f f K/b2, with f = c1 b2, functions of f^2 = c1^2 (k1^2 +
// k1s^2). The map is then regular at zero strength, with TPSA derivatives.
// g1 carries MAD-NG's lw (+1 on the face crossed first, -1 on the other) and g2
// the direction sdir, so that each backtracking call inverts the forward one
// (as in MAD-NG's strex_fringe).
GPUFUN XT_NOINLINE
void MagnetEdge_quad_lw_fringe_single_particle(
    LocalParticle* part,
    xt_float_or_tpsa_arg k1,   // total quadrupole strength
    xt_float_or_tpsa_arg k1s,  // total skew quadrupole strength
    const double g1,
    const double g2
) {
    const double chi = LocalParticle_get_chi(part);
    const double beta0 = LocalParticle_get_beta0(part);

    const xt_float_or_tpsa x = LocalParticle_get_x(part);
    const xt_float_or_tpsa px = LocalParticle_get_px(part);
    const xt_float_or_tpsa y = LocalParticle_get_y(part);
    const xt_float_or_tpsa py = LocalParticle_get_py(part);
    xt_float_or_tpsa t = LocalParticle_get_zeta(part) / beta0;
    const xt_float_or_tpsa pt = LocalParticle_get_ptau(part);

    const xt_float_or_tpsa _pz = 1 / sqrt(1 + 2 / beta0 * pt + POW2(pt));
    const xt_float_or_tpsa dt = (1 / beta0 + pt) * _pz;

    const xt_float_or_tpsa c1 = chi * g1 * _pz;
    const xt_float_or_tpsa c2 = chi * g2 * _pz;

    const xt_float_or_tpsa kk = POW2(k1) + POW2(k1s);
    xt_float_or_tpsa cosh_f = 1., sinhc_f = 1.;
    MagnetEdge_cosh_sinhc(POW2(c1) * kk, &cosh_f, &sinhc_f);

    // K q and K p
    const xt_float_or_tpsa kq_x = k1 * x - k1s * y;
    const xt_float_or_tpsa kq_y = -k1s * x - k1 * y;
    const xt_float_or_tpsa kp_x = k1 * px - k1s * py;
    const xt_float_or_tpsa kp_y = -k1s * px - k1 * py;

    // Time of flight: q.K p, p.K p, p.p
    const xt_float_or_tpsa a_s = cosh_f - POW2(c1) * kk / 2 * sinhc_f;
    const xt_float_or_tpsa sum =
        c1 * (x * kp_x + y * kp_y) + a_s * _pz * c2 * (px * kp_x + py * kp_y)
        + _pz * c1 * c2 * kk * (cosh_f / 2 - sinhc_f) * (POW2(px) + POW2(py));
    t -= dt * sum * _pz;

    LocalParticle_set_x(part, cosh_f * x + sinhc_f * c1 * kq_x + c2 * _pz * kp_x);
    LocalParticle_set_y(part, cosh_f * y + sinhc_f * c1 * kq_y + c2 * _pz * kp_y);
    LocalParticle_set_px(part, cosh_f * px - sinhc_f * c1 * kp_x);
    LocalParticle_set_py(part, cosh_f * py - sinhc_f * c1 * kp_y);
    LocalParticle_set_zeta(part, t * beta0);
}



// The Lee-Whiting fringe at a face inside a solenoid field bz: the map acts on
// the kinetic momenta. Ported from SAD, subroutine tqlfre, bz != 0 branch:
//     https://github.com/KatsOide/SAD/blob/d082c7c5161ea111132e12c0cf9288fb892c23ba/src/tqlfre.f#L83
// SAD solves the implicit relation of the map on the kinetic momenta in closed
// form (the factor f below), so that it stays symplectic in the canonical
// variables; at bz = 0 it is MagnetEdge_quad_lw_fringe_single_particle. Its
// inverse (backtracking) is the same map with g1 and g2 flipped, as there.
// Adapted to xtrack:
//   - SAD's frame is that of an upright quadrupole: a skew gradient rotates
//     the coordinates by theta = atan2(k1s, k1) / 2 (cos and sin from the half
//     angle) to the frame where the gradient is b2 = (k1^2 + k1s^2)^(1/2), and
//     back. An upright quadrupole (k1s = 0) is not rotated, and keeps the sign
//     of k1: regular at k1 = 0, also for TPSA. The vector potential bz/2 (-y, x)
//     does not change under the rotation.
//   - The momenta are canonical in the symmetric gauge of the solenoid faces,
//     centred on the solenoid axis (x0, y0): shifting p by the constant
//     bz/2 (-y0, x0) (canonical) centres it on the magnet axis, as in SAD.
//   - SAD's z is zeta * rvv (conjugate to delta): zeta changes by dz / rvv.
//   - SAD's ak/al f1 |f1| / 24 and ak/al f2 are chi b2 g1 and chi b2 g2, with
//     g1, g2 as for MagnetEdge_quad_lw_fringe_single_particle (signed for the
//     face and the direction of travel).
GPUFUN XT_NOINLINE
void MagnetEdge_quad_lw_fringe_solenoid_single_particle(
    LocalParticle* part,
    xt_float_or_tpsa_arg k1,   // total quadrupole strength
    xt_float_or_tpsa_arg k1s,  // total skew quadrupole strength
    const double g1,
    const double g2,
    xt_float_or_tpsa_arg bz,   // solenoid strength at the face
    const double x0,
    const double y0
) {
    const double chi = LocalParticle_get_chi(part);
    const xt_float_or_tpsa rvv = LocalParticle_get_rvv(part);
    const xt_float_or_tpsa p = 1 + LocalParticle_get_delta(part);
    const xt_float_or_tpsa x_in = LocalParticle_get_x(part);
    const xt_float_or_tpsa y_in = LocalParticle_get_y(part);
    const xt_float_or_tpsa px_in = LocalParticle_get_px(part) - bz / 2 * y0;
    const xt_float_or_tpsa py_in = LocalParticle_get_py(part) + bz / 2 * x0;

    // Frame of the quadrupole: rotation by theta
    xt_float_or_tpsa b2 = k1, c = 1., s = 0.;
    if (!xt_float_or_tpsa_is_zero(k1s)) {
        b2 = sqrt(POW2(k1) + POW2(k1s));
        c = sqrt((1 + k1 / b2) / 2);
        s = k1s / (2 * b2 * c);
    }
    const xt_float_or_tpsa x = c * x_in - s * y_in;
    const xt_float_or_tpsa y = s * x_in + c * y_in;
    const xt_float_or_tpsa px = c * px_in - s * py_in;
    const xt_float_or_tpsa py = s * px_in + c * py_in;

    // SAD tqlfre
    const xt_float_or_tpsa a = chi * b2 * g1 / p;
    const xt_float_or_tpsa bp = chi * b2 * g2 / POW2(p);
    const xt_float_or_tpsa ea = exp(a);
    const xt_float_or_tpsa bzph = bz / 2;
    const xt_float_or_tpsa bb = bp * bzph;
    const xt_float_or_tpsa f = 1 / ((1 - bb) * (1 + bb));
    const xt_float_or_tpsa xf = (ea * x + bp * (px + bzph * y / ea - bb * py)) * f;
    const xt_float_or_tpsa yf = (y / ea - bp * (py - bzph * x * ea - bb * px)) * f;
    const xt_float_or_tpsa pxf = (px + bzph * yf) / ea;
    const xt_float_or_tpsa pyf = (py - bzph * xf) * ea;
    const xt_float_or_tpsa dz =
        -((a * ea * x + bp * ea * (1 + 0.5 * a) * pxf) * pxf
          - (a * y / ea + bp / ea * (1 - 0.5 * a) * pyf) * pyf) / p;
    const xt_float_or_tpsa px_out = pxf - bzph * y;
    const xt_float_or_tpsa py_out = pyf + bzph * x;

    // Back to the frame of the magnet
    LocalParticle_set_x(part, c * xf + s * yf);
    LocalParticle_set_y(part, -s * xf + c * yf);
    LocalParticle_set_px(part, c * px_out + s * py_out + bz / 2 * y0);
    LocalParticle_set_py(part, -s * px_out + c * py_out - bz / 2 * x0);
    LocalParticle_add_to_zeta(part, dz / rvv);
}



// A strength in the plan (a field marked FT, magnet_plans.py), from its
// address: a double, or in TPSA builds the TPSA held in place
// (MagnetPlanTpsaStorage, track_magnet.h)
#ifdef XTRACK_TPSA_TRACK
    #define XT_FT(ptr) (*reinterpret_cast<xt_float_or_tpsa*>(ptr))
#else
    #define XT_FT(ptr) (*(ptr))
#endif
#define XT_FT_PTR(ptr) (&XT_FT(ptr))

// The arrays the plan refers to are held by address: the multipoles (of the
// element or knl_tot of the plan), the radiation record
#define MAGNET_PLAN_ADDR(ptr) ((int64_t)(size_t)(ptr))

GPUFUN XT_INLINE
const xt_float_or_tpsa* MagnetField_knl(const MagnetField field) {
    return (const xt_float_or_tpsa*)(size_t)MagnetField_get_knl(field);
}

GPUFUN XT_INLINE
const xt_float_or_tpsa* MagnetField_ksl(const MagnetField field) {
    return (const xt_float_or_tpsa*)(size_t)MagnetField_get_ksl(field);
}


// Total dipole and quadrupole strengths: k + knl / length
GPUFUN XT_INLINE
void MagnetField_prepare(MagnetField field) {
    const double length = MagnetField_get_length(field);
    const int64_t k_order = MagnetField_get_k_order(field);
    const int64_t kl_order = MagnetField_get_kl_order(field);
    const double knl_scale = MagnetField_get_knl_scale(field);
    const xt_float_or_tpsa* knl = MagnetField_knl(field);
    const xt_float_or_tpsa* ksl = MagnetField_ksl(field);
    xt_float_or_tpsa k0 = 0., k1 = 0., k1s = 0.;
    if (k_order > -1) k0 += XT_FT(MagnetField_getp1_knorm(field, 0));
    if (k_order > 0) {
        k1 += XT_FT(MagnetField_getp1_knorm(field, 1));
        k1s += XT_FT(MagnetField_getp1_kskew(field, 1));
    }
    if (fabs(length) > 1e-10 && kl_order > -1) {
        k0 += knl_scale * knl[0] / length;
        if (kl_order > 0) {
            k1 += knl_scale * knl[1] / length;
            k1s += knl_scale * ksl[1] / length;
        }
    }
    XT_FT(MagnetField_getp_k0(field)) = k0;
    XT_FT(MagnetField_getp_k1(field)) = k1;
    XT_FT(MagnetField_getp_k1s(field)) = k1s;

    uint8_t has_mult = 0;
    for (int64_t ii = 1; ii <= k_order && !has_mult; ii++) {
        has_mult = !xt_float_or_tpsa_is_zero(XT_FT(MagnetField_getp1_knorm(field, ii)))
                   || !xt_float_or_tpsa_is_zero(XT_FT(MagnetField_getp1_kskew(field, ii)));
    }
    for (int64_t ii = 1; ii <= kl_order && !has_mult; ii++) {
        has_mult = !xt_float_or_tpsa_is_zero(knl[ii]) || !xt_float_or_tpsa_is_zero(ksl[ii]);
    }
    MagnetField_set_has_mult(field, has_mult);
}



// Reset the settings of a face: all fringes switched off and no patch
GPUFUN
void MagnetFace_reset(MagnetFace face, const uint8_t is_exit_face) {
    MagnetFace_set_flags(face, 0);
    MagnetFace_set_is_exit_face(face, is_exit_face);
    MagnetFace_set_angle(face, 0);
    MagnetFace_set_angle_fdown(face, 0);
    MagnetFace_set_fint(face, 0);
    MagnetFace_set_hgap(face, 0);
    MagnetFace_set_f1(face, 0);
    MagnetFace_set_f2(face, 0);
    MagnetFace_set_has_ks(face, 0);
    MagnetFace_set_ks(face, 0);
    MagnetFace_set_has_patch(face, 0);
    MagnetFace_set_patch_angle(face, 0);
    MagnetFace_set_patch_dx(face, 0);
}


// Whether a face is active with a sub-map that uses the derived field
// (MagnetField_prepare): the linear edge and the fringes. The patch and the
// vector potential do not.
GPUFUN XT_INLINE
uint8_t MagnetFace_uses_field(const MagnetFace face) {
    const int64_t flags = MagnetFace_get_flags(face);
    return (flags & EDGE_ACTIVE) && (flags & (EDGE_LINEAR | EDGE_BEND | EDGE_MULT | EDGE_QUAD_LW));
}


// Solenoid strength of the vector potential at a face: its own, or the one
// of the field
GPUFUN XT_INLINE
xt_float_or_tpsa MagnetFace_ks(const MagnetFace face, const MagnetField field) {
    if (MagnetFace_get_has_ks(face)) return MagnetFace_get_ks(face);
    return XT_FT(MagnetField_getp_ks(field));
}


// Whether crossing a face applies anything: the face is active and has a
// sub-map (a model other than "suppressed"), a patch or a solenoid field. An
// active face without them is treated like an inactive one. Needs the ks of
// the field. A face with its own ks always sets the vector potential, also to
// zero: the other face of the element can have a different one.
GPUFUN XT_INLINE
uint8_t MagnetFace_needs_work(const MagnetFace face, const MagnetField field) {
    const int64_t flags = MagnetFace_get_flags(face);
    if (!(flags & EDGE_ACTIVE)) return 0;
    return (flags & (EDGE_LINEAR | EDGE_BEND | EDGE_MULT | EDGE_QUAD_LW))
           || MagnetFace_get_has_patch(face) || MagnetFace_get_has_ks(face)
           || xt_float_or_tpsa_const_part(XT_FT(MagnetField_getp_ks(field))) != 0.;
}


// Derive the sub-maps of a face for the direction of travel. The field must
// have been prepared (MagnetField_prepare) if MagnetFace_uses_field.
GPUFUN XT_INLINE
void MagnetFace_prepare(MagnetFace face, MagnetField field, const uint8_t backtrack) {
    const int64_t flags = MagnetFace_get_flags(face);
    const uint8_t is_exit_face = MagnetFace_get_is_exit_face(face);
    const double angle = MagnetFace_get_angle(face);
    // Rotations and wedges are inverted by flipping the angle sign
    const double rot_sign = backtrack ? 1. : -1.;

    MagnetFace_set_backtrack(face, backtrack);
    const uint8_t has_work = MagnetFace_needs_work(face, field);
    MagnetFace_set_has_work(face, has_work);
    if (!has_work) {
        MagnetFace_set_do_patch(face, 0);
        MagnetFace_set_do_vector_potential(face, 0);
        MagnetFace_set_do_linear(face, 0);
        MagnetFace_set_do_rotate(face, 0);
        MagnetFace_set_do_bend(face, 0);
        MagnetFace_set_do_mult(face, 0);
        MagnetFace_set_do_quad_lw(face, 0);
        MagnetFace_set_do_quad_wedge(face, 0);
        MagnetFace_set_do_wedge(face, 0);
        return;
    }

    // Patch
    const uint8_t has_patch = MagnetFace_get_has_patch(face);
    const double patch_angle = MagnetFace_get_patch_angle(face);
    MagnetFace_set_do_patch(face, has_patch);
    double patch_sin = 0, patch_cos = 1;
    if (has_patch && fabs(patch_angle) > EDGE_ANGLE_EPS) {
        patch_sin = sin(patch_angle);
        patch_cos = cos(patch_angle);
    }
    MagnetFace_set_patch_sin(face, rot_sign * patch_sin);
    MagnetFace_set_patch_cos(face, patch_cos);
    MagnetFace_set_patch_tan(face, rot_sign * patch_sin / patch_cos);

    MagnetFace_set_do_vector_potential(face, 1);

    // Linear edge: a kick, inverted by flipping it
    const uint8_t do_linear = (flags & EDGE_LINEAR) != 0;
    MagnetFace_set_do_linear(face, do_linear);
    if (do_linear) {
        xt_float_or_tpsa* r21 = XT_FT_PTR(MagnetField_getp1_r21(field, is_exit_face));
        xt_float_or_tpsa* r43 = XT_FT_PTR(MagnetField_getp1_r43(field, is_exit_face));
        compute_dipole_edge_linear_coefficients(
            XT_FT(MagnetField_getp_k0(field)), angle, MagnetFace_get_angle_fdown(face),
            MagnetFace_get_hgap(face), MagnetFace_get_fint(face), r21, r43
        );
        if (backtrack) {
            *r21 = -*r21;
            *r43 = -*r43;
        }
    }

    // Magnet fringes; no geometric sub-map without a fringe map
    const uint8_t fringe =
        !(flags & EDGE_LINEAR) && (flags & (EDGE_BEND | EDGE_MULT | EDGE_QUAD_LW));
    const int64_t k_order = MagnetField_get_k_order(field);
    // A fringe of a vanishing field is the identity: skipped (the dipole
    // fringe uses the threshold of DipoleFringe_track_single_particle)
    const uint8_t do_bend = fringe && (flags & EDGE_BEND)
                            && !xt_float_or_tpsa_is_negligible(XT_FT(MagnetField_getp_k0(field)), EDGE_K0_EPS);
    const uint8_t do_mult = fringe && (flags & EDGE_MULT) && MagnetField_get_has_mult(field);
    const uint8_t do_rotate = fringe && fabs(angle) > EDGE_ANGLE_EPS;
    MagnetFace_set_do_bend(face, do_bend);
    MagnetFace_set_do_mult(face, do_mult);
    MagnetFace_set_do_rotate(face, do_rotate);
    MagnetFace_set_do_wedge(face, do_rotate && (k_order >= 0));
    MagnetFace_set_do_quad_wedge(face, do_mult && do_rotate && (k_order >= 1));
    // Soft quadrupole fringe (Lee-Whiting), with or without the hard-edge one,
    // and with the linear edge too
    const double f1 = MagnetFace_get_f1(face);
    const double f2 = MagnetFace_get_f2(face);
    MagnetFace_set_do_quad_lw(
        face,
        (flags & EDGE_QUAD_LW) && (f1 != 0 || f2 != 0)
            && !(xt_float_or_tpsa_is_zero(XT_FT(MagnetField_getp_k1(field)))
                 && xt_float_or_tpsa_is_zero(XT_FT(MagnetField_getp_k1s(field))))
    );

    double sin_ = 0, cos_ = 1, tan_ = 0;
    if (do_rotate) {
        sin_ = sin(angle);
        cos_ = cos(angle);
        tan_ = tan(angle);
    }
    MagnetFace_set_rot_sin(face, rot_sign * sin_);
    MagnetFace_set_rot_cos(face, cos_);
    MagnetFace_set_rot_tan(face, rot_sign * tan_);
    MagnetFace_set_wedge_angle(face, rot_sign * angle);

    // +1 on the face crossed first, -1 on the other (MAD-NG lw)
    const uint8_t entering = (is_exit_face == backtrack);
    MagnetFace_set_quad_lw_g1(face, (entering ? 1. : -1.) * (fabs(f1) * f1 / -24));
    MagnetFace_set_quad_lw_g2(face, f2 * (backtrack ? -1. : 1.));
}


// The sub-maps of a face, each the forward map of the face or, when
// backtracking, its inverse (the signs are set in MagnetFace_prepare).
// Inlined: they take the face by handle; the fringe maps they call take
// scalars and are not inlined (XT_INLINE, XT_NOINLINE: headers/track.h).

GPUFUN XT_INLINE
void MagnetFace_patch_single_particle(
    LocalParticle* part,
    const MagnetFace f,
    const uint8_t entering
) {
    if (!MagnetFace_get_do_patch(f)) return;
    if (entering) {
        YRotation_single_particle(
            part, MagnetFace_get_patch_sin(f), MagnetFace_get_patch_cos(f),
            MagnetFace_get_patch_tan(f)
        );
        LocalParticle_add_to_x(part, MagnetFace_get_patch_dx(f));
    }
    else {
        LocalParticle_add_to_x(part, -MagnetFace_get_patch_dx(f));
        YRotation_single_particle(
            part, MagnetFace_get_patch_sin(f), MagnetFace_get_patch_cos(f),
            MagnetFace_get_patch_tan(f)
        );
    }
}

GPUFUN XT_INLINE
void MagnetFace_vector_potential_single_particle(
    LocalParticle* part,
    const MagnetFace f,
    const MagnetField field,
    const uint8_t entering
) {
    if (!MagnetFace_get_do_vector_potential(f)) return;
    // A face with its own ks (VariableSolenoid) is not the end of the field:
    // the vector potential of its ks is set when leaving too
    MagnetEdge_vector_potential_single_particle(
        part,
        /*is_leaving_magnet*/ !entering && !MagnetFace_get_has_ks(f),
        MagnetFace_ks(f, field),
        MagnetField_get_x0_solenoid(field),
        MagnetField_get_y0_solenoid(field)
    );
}

GPUFUN XT_INLINE
void MagnetFace_linear_single_particle(
    LocalParticle* part,
    const MagnetFace f,
    const MagnetField field,
    double taper
) {
    if (!MagnetFace_get_do_linear(f)) return;
    const uint8_t is_exit_face = MagnetFace_get_is_exit_face(f);
    xt_float_or_tpsa r21 = XT_FT(MagnetField_getp1_r21(field, is_exit_face));
    xt_float_or_tpsa r43 = XT_FT(MagnetField_getp1_r43(field, is_exit_face));
    if (taper != 1.) {
        // Tapered: r43 is not linear in k0, both are recomputed
        compute_dipole_edge_linear_coefficients(
            taper * XT_FT(MagnetField_getp_k0(field)), MagnetFace_get_angle(f),
            MagnetFace_get_angle_fdown(f), MagnetFace_get_hgap(f), MagnetFace_get_fint(f),
            &r21, &r43
        );
        if (MagnetFace_get_backtrack(f)) {
            r21 = -r21;
            r43 = -r43;
        }
    }
    DipoleEdgeLinear_single_particle(part, r21, r43);
}

GPUFUN XT_INLINE
void MagnetFace_rotation_single_particle(LocalParticle* part, const MagnetFace f) {
    if (!MagnetFace_get_do_rotate(f)) return;
    YRotation_single_particle(
        part, MagnetFace_get_rot_sin(f), MagnetFace_get_rot_cos(f), MagnetFace_get_rot_tan(f)
    );
}

GPUFUN XT_INLINE
void MagnetFace_bend_fringe_single_particle(
    LocalParticle* part,
    const MagnetFace f,
    const MagnetField field,
    double taper
) {
    if (!MagnetFace_get_do_bend(f)) return;
    // The fringe kick changes sign between the two faces, a property of the
    // face itself
    if (MagnetFace_get_is_exit_face(f)) {
        DipoleFringe_track_single_particle(
            part, MagnetFace_get_fint(f), MagnetFace_get_hgap(f),
            -taper * XT_FT(MagnetField_getp_k0(field)), MagnetFace_get_backtrack(f)
        );
    }
    else {
        DipoleFringe_track_single_particle(
            part, MagnetFace_get_fint(f), MagnetFace_get_hgap(f),
            taper * XT_FT(MagnetField_getp_k0(field)), MagnetFace_get_backtrack(f)
        );
    }
}

GPUFUN XT_INLINE
void MagnetFace_mult_fringe_single_particle(
    LocalParticle* part,
    const MagnetFace f,
    const MagnetField field,
    double taper
) {
    if (!MagnetFace_get_do_mult(f)) return;
    MultFringe_track_single_particle(
        part,
        XT_FT_PTR(MagnetField_getp1_knorm(field, 0)),
        XT_FT_PTR(MagnetField_getp1_kskew(field, 0)),
        MagnetField_get_k_order(field),
        MagnetField_knl(field),
        MagnetField_ksl(field),
        MagnetField_get_kl_order(field),
        MagnetField_get_knl_scale(field),
        MagnetField_get_length(field),
        taper,
        MagnetFace_get_is_exit_face(f),
        /* min_order */ 1
    );
}

GPUFUN XT_INLINE
void MagnetFace_quad_lw_fringe_single_particle(
    LocalParticle* part,
    const MagnetFace f,
    const MagnetField field,
    double taper
) {
    if (!MagnetFace_get_do_quad_lw(f)) return;
    const xt_float_or_tpsa bz = MagnetFace_ks(f, field);
    if (xt_float_or_tpsa_is_zero(bz)) {
        MagnetEdge_quad_lw_fringe_single_particle(
            part, taper * XT_FT(MagnetField_getp_k1(field)), taper * XT_FT(MagnetField_getp_k1s(field)),
            MagnetFace_get_quad_lw_g1(f), MagnetFace_get_quad_lw_g2(f)
        );
    }
    else {
        MagnetEdge_quad_lw_fringe_solenoid_single_particle(
            part, taper * XT_FT(MagnetField_getp_k1(field)), taper * XT_FT(MagnetField_getp_k1s(field)),
            MagnetFace_get_quad_lw_g1(f), MagnetFace_get_quad_lw_g2(f), bz,
            MagnetField_get_x0_solenoid(field), MagnetField_get_y0_solenoid(field)
        );
    }
}

GPUFUN XT_INLINE
void MagnetFace_quad_wedge_single_particle(
    LocalParticle* part,
    const MagnetFace f,
    const MagnetField field,
    double taper
) {
    if (!MagnetFace_get_do_quad_wedge(f)) return;
    Quad_wedge_single_particle(
        part, MagnetFace_get_wedge_angle(f), taper * XT_FT(MagnetField_getp1_knorm(field, 1))
    );
}

GPUFUN XT_INLINE
void MagnetFace_wedge_single_particle(
    LocalParticle* part,
    const MagnetFace f,
    const MagnetField field,
    double taper
) {
    if (!MagnetFace_get_do_wedge(f)) return;
    // sin, cos and tan of wedge_angle are those of the face rotation
    Wedge_trig_single_particle(
        part, MagnetFace_get_rot_sin(f), MagnetFace_get_rot_cos(f), MagnetFace_get_rot_tan(f),
        taper * XT_FT(MagnetField_getp1_knorm(field, 0))
    );
}


// Cross a face from the drift into the field (the order at the top of this
// file)
GPUFUN XT_INLINE
void MagnetFace_enter_single_particle(
    LocalParticle* part,
    const MagnetFace f,
    const MagnetField field,
    double taper  // factor of the strengths (tapering, 1 otherwise)
) {
    MagnetFace_patch_single_particle(part, f, 1);
    MagnetFace_vector_potential_single_particle(part, f, field, 1);
    MagnetFace_linear_single_particle(part, f, field, taper);
    MagnetFace_rotation_single_particle(part, f);
    MagnetFace_bend_fringe_single_particle(part, f, field, taper);
    MagnetFace_mult_fringe_single_particle(part, f, field, taper);
    MagnetFace_quad_lw_fringe_single_particle(part, f, field, taper);
    MagnetFace_quad_wedge_single_particle(part, f, field, taper);
    MagnetFace_wedge_single_particle(part, f, field, taper);
}


// Cross a face from the field into the drift: the same sub-maps in reverse
GPUFUN XT_INLINE
void MagnetFace_leave_single_particle(
    LocalParticle* part,
    const MagnetFace f,
    const MagnetField field,
    double taper  // factor of the strengths (tapering, 1 otherwise)
) {
    MagnetFace_wedge_single_particle(part, f, field, taper);
    MagnetFace_quad_wedge_single_particle(part, f, field, taper);
    MagnetFace_quad_lw_fringe_single_particle(part, f, field, taper);
    MagnetFace_mult_fringe_single_particle(part, f, field, taper);
    MagnetFace_bend_fringe_single_particle(part, f, field, taper);
    MagnetFace_rotation_single_particle(part, f);
    MagnetFace_linear_single_particle(part, f, field, taper);
    MagnetFace_vector_potential_single_particle(part, f, field, 0);
    MagnetFace_patch_single_particle(part, f, 0);
}

#endif // XTRACK_TRACK_MAGNET_EDGE_H
