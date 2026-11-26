// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2023.                 //
// ######################################### //
#ifndef XTRACK_TRACK_MAGNET_H
#define XTRACK_TRACK_MAGNET_H

#include "xtrack/headers/track.h"
#include "xtrack/beam_elements/elements_src/track_magnet_structs.h"
#include "xtrack/beam_elements/elements_src/track_magnet_kick.h"
#include "xtrack/beam_elements/elements_src/track_magnet_drift.h"
#include "xtrack/beam_elements/elements_src/track_magnet_edge.h"
#include "xtrack/beam_elements/elements_src/edge_flags.h"
#include "xtrack/beam_elements/elements_src/track_magnet_configure.h"
#include "xtrack/beam_elements/elements_src/integrator.h"
#include "xtrack/beam_elements/elements_src/track_misalignments.h"

#ifndef XTRACK_MULTIPOLE_NO_SYNRAD
#include "xtrack/beam_elements/elements_src/track_magnet_radiation.h"
#endif

// Order of knl + rel_ref_strength * knl_rel: the order of knl, or the highest
// nonzero relative component above it
GPUFUN
int64_t magnet_multipole_order(
    const int64_t order,
    const int64_t order_rel,
    GPUGLMEM const double* knl_rel,
    GPUGLMEM const double* ksl_rel
) {
    if (knl_rel == NULL) return order;
    for (int64_t ii = order_rel; ii > order; ii--) {
        if (knl_rel[ii] != 0 || ksl_rel[ii] != 0) return ii;
    }
    return order;
}

// Whether knl_rel or ksl_rel has a nonzero component (if not, order_tot is
// the order of knl)
GPUFUN
uint8_t magnet_has_rel(
    const int64_t order_rel,
    GPUGLMEM const double* knl_rel,
    GPUGLMEM const double* ksl_rel
) {
    if (knl_rel == NULL) return 0;
    for (int64_t ii = 0; ii <= order_rel; ii++) {
        if (knl_rel[ii] != 0 || ksl_rel[ii] != 0) return 1;
    }
    return 0;
}

// knl + rel_ref_strength * knl_rel (same for ksl), orders 0 ... order_tot,
// into knl_out, ksl_out
GPUFUN
void magnet_total_multipoles(
    const int64_t order,
    GPUGLMEM const double* knl,
    GPUGLMEM const double* ksl,
    const int64_t order_rel,
    GPUGLMEM const double* knl_rel,
    GPUGLMEM const double* ksl_rel,
    xt_float_or_tpsa_arg rel_ref_strength,
    const int64_t order_tot,
    xt_float_or_tpsa* knl_out,
    xt_float_or_tpsa* ksl_out
) {
    for (int64_t ii = 0; ii <= order_tot; ii++) {
        knl_out[ii] = ii <= order ? knl[ii] : 0.;
        ksl_out[ii] = ii <= order ? ksl[ii] : 0.;
    }
    if (knl_rel == NULL) return;
    for (int64_t ii = 0; ii <= order_rel && ii <= order_tot; ii++) {
        if (knl_rel[ii] != 0) knl_out[ii] += rel_ref_strength * knl_rel[ii];
        if (ksl_rel[ii] != 0) ksl_out[ii] += rel_ref_strength * ksl_rel[ii];
    }
}

#ifndef VSWAP
    #define VSWAP(a, b)     \
        {                   \
            double tmp = a; \
            a = b;          \
            b = tmp;        \
        }
#endif


// The radiation record of a body (held by address in the plan)
GPUFUN XT_INLINE
MagnetRadiationRecord MagnetBody_radiation_record(const MagnetBody body) {
    return (MagnetRadiationRecord)(size_t)MagnetBody_get_radiation_record(body);
}


// Apply the drift of a body over dlength, with taper the factor of the
// strengths (tapering, 1 otherwise)
GPUFUN XT_INLINE
void MagnetDrift_apply(
    LocalParticle* part,
    const MagnetDrift drift,
    double length,
    double taper
) {
    track_magnet_drift_single_particle(
        part,
        length,
        taper * XT_FT(MagnetDrift_getp_k0(drift)),
        taper * XT_FT(MagnetDrift_getp_k1(drift)),
        XT_FT(MagnetDrift_getp_ks(drift)),
        MagnetDrift_get_h(drift),
        MagnetDrift_get_x0_solenoid(drift),
        MagnetDrift_get_y0_solenoid(drift),
        MagnetDrift_get_model(drift)
    );
}


// Apply a prepared kick, with kick_weight the fraction of the body it covers
// and taper the factor of the strengths (tapering, 1 otherwise); its
// multipoles are those of the field
GPUFUN XT_INLINE
void MagnetKick_apply(
    LocalParticle* part,
    const MagnetKick kick,
    const MagnetField field,
    double kick_weight,
    double taper
) {
    track_magnet_kick_single_particle(
        part,
        kick_weight,
        taper,
        MagnetKick_get_length(kick),
        MagnetKick_get_order(kick),
        MagnetKick_get_inv_factorial_order(kick),
        MagnetField_knl(field),
        MagnetField_ksl(field),
        MagnetKick_get_factor_knl_ksl(kick),
        XT_FT(MagnetKick_getp1_knl_main(kick, 0)),
        XT_FT(MagnetKick_getp1_knl_main(kick, 1)),
        XT_FT(MagnetKick_getp1_knl_main(kick, 2)),
        XT_FT(MagnetKick_getp1_knl_main(kick, 3)),
        XT_FT(MagnetKick_getp1_ksl_main(kick, 0)),
        XT_FT(MagnetKick_getp1_ksl_main(kick, 1)),
        XT_FT(MagnetKick_getp1_ksl_main(kick, 2)),
        XT_FT(MagnetKick_getp1_ksl_main(kick, 3)),
        MagnetKick_get_h(kick),
        MagnetKick_get_hxl(kick),
        MagnetKick_get_htot(kick),
        XT_FT(MagnetKick_getp_k0_h_coeff(kick)),
        XT_FT(MagnetKick_getp_k1_h_coeff(kick)),
        MagnetKick_get_rot_frame(kick),
        MagnetKick_get_do_drift_correction(kick),
        MagnetKick_get_is_empty(kick),
        MagnetKick_get_has_main(kick),
        MagnetKick_get_has_curvature(kick)
    );
}


#if !defined(XTRACK_TPSA_TRACK) && !defined(XTRACK_MULTIPOLE_NO_SYNRAD)
// Spin and radiation over a step of length ll of the body, from the
// coordinates before it (old_*) and after it (the particle).
GPUFUN XT_NOINLINE
void magnet_body_spin_radiation(
    LocalParticle* part,
    const MagnetBody body,
    const MagnetField field,
    double const taper,
    double const ll,
    double const old_x,
    double const old_y,
    double const old_zeta,
    double const old_kin_px,
    double const old_kin_py,
    double* dp_record_exit,
    double* dpx_record_exit,
    double* dpy_record_exit
) {
    const MagnetDrift drift = MagnetBody_getp_drift(body);
    double Bx_T, By_T, Bz_T;
    double const new_kin_px = LocalParticle_get_px(part) - LocalParticle_get_ax(part);
    double const new_kin_py = LocalParticle_get_py(part) - LocalParticle_get_ay(part);
    double const mean_kin_px = 0.5 * (old_kin_px + new_kin_px);
    double const mean_kin_py = 0.5 * (old_kin_py + new_kin_py);
    evaluate_field_from_strengths(
        LocalParticle_get_p0c(part),
        LocalParticle_get_q0(part),
        0.5 * (old_x + LocalParticle_get_x(part)),
        0.5 * (old_y + LocalParticle_get_y(part)),
        MagnetBody_get_length(body),
        MagnetField_get_kl_order(field),
        MagnetBody_get_inv_factorial_order(body),
        MagnetField_knl(field),
        MagnetField_ksl(field),
        taper * MagnetBody_get_factor_knl_ksl(body),
        taper * XT_FT(MagnetField_getp1_knorm(field, 0)),
        taper * XT_FT(MagnetField_getp1_knorm(field, 1)),
        taper * XT_FT(MagnetField_getp1_knorm(field, 2)),
        taper * XT_FT(MagnetField_getp1_knorm(field, 3)),
        taper * XT_FT(MagnetField_getp1_kskew(field, 0)),
        taper * XT_FT(MagnetField_getp1_kskew(field, 1)),
        taper * XT_FT(MagnetField_getp1_kskew(field, 2)),
        taper * XT_FT(MagnetField_getp1_kskew(field, 3)),
        XT_FT(MagnetDrift_getp_ks(drift)),
        MagnetBody_get_dks_ds(body),
        MagnetDrift_get_x0_solenoid(drift),
        MagnetDrift_get_y0_solenoid(drift),
        &Bx_T,
        &By_T,
        &Bz_T
    );
    double const dzeta = LocalParticle_get_zeta(part) - old_zeta;
    double const l_path = LocalParticle_get_rvv(part) * (ll - dzeta);
    magnet_spin(part, Bx_T, By_T, Bz_T, MagnetBody_get_h_radiation(body), ll, l_path);
    int64_t const radiation_flag = MagnetBody_get_radiation_flag(body);
    if (radiation_flag) {
        double const B_perp_T = compute_b_perp_mod(
            mean_kin_px, mean_kin_py, LocalParticle_get_delta(part), Bx_T, By_T, Bz_T
        );
        magnet_radiation(
            part, B_perp_T, ll, l_path, radiation_flag, MagnetBody_radiation_record(body),
            dp_record_exit, dpx_record_exit, dpy_record_exit
        );
    }
}
#endif


GPUFUN XT_INLINE
void track_magnet_body_single_particle(
    LocalParticle* part,
    const MagnetBody body,
    const MagnetField field,
    double taper,  // factor of the strengths (tapering, 1 otherwise)
    double* dp_record_exit,
    double* dpx_record_exit,
    double* dpy_record_exit
) {
    const double length = MagnetBody_get_length(body);

// The maps for RUN_INTEGRATOR
#define MAGNET_KICK(part, weight)                                              \
    do {                                                                       \
        if (!MagnetBody_get_kick_is_noop(body))                                \
            MagnetKick_apply(part, MagnetBody_getp_kick(body), field, (weight), taper); \
    } while (0)
#define MAGNET_DRIFT(part, dlength) \
    MagnetDrift_apply(part, MagnetBody_getp_drift(body), (dlength), taper)

#if defined(XTRACK_TPSA_TRACK) || defined(XTRACK_MULTIPOLE_NO_SYNRAD)
    #define WITH_RADIATION(ll, code) \
        {                            \
            code;                    \
        }
    #else
        // Spin and radiation of a step, after it: one shared copy (not
        // inlined), not one per integrator step and magnet class
        #define WITH_RADIATION(ll, code) \
        { \
            double const old_x = LocalParticle_get_x(part); \
            double const old_y = LocalParticle_get_y(part); \
            double const old_zeta = LocalParticle_get_zeta(part); \
            double const old_kin_px = LocalParticle_get_px(part) - LocalParticle_get_ax(part); \
            double const old_kin_py = LocalParticle_get_py(part) - LocalParticle_get_ay(part); \
            code; \
            if (length > 0) { /* spin always, radiation if radiation_flag */ \
                magnet_body_spin_radiation( \
                    part, body, field, taper, (ll), old_x, old_y, old_zeta, old_kin_px, \
                    old_kin_py, dp_record_exit, dpx_record_exit, dpy_record_exit \
                ); \
            } \
        }
    #endif

    if (MagnetBody_get_drift_only(body)) {
        WITH_RADIATION(length,
            MAGNET_DRIFT(part, length);
        )
    }
    else{

        RUN_INTEGRATOR(
            MagnetBody_get_integrator(body), length, MagnetBody_get_num_kicks(body), part,
            MAGNET_DRIFT, MAGNET_KICK, WITH_RADIATION
        );

    }

    #undef MAGNET_KICK
    #undef MAGNET_DRIFT
    #undef WITH_RADIATION
}


// Straight-body RBend: the body is a straight line and the curvature is in
// the patches of the faces. At the entry face the reference orbit makes the
// angle theta_entry = (angle - angle_diff) / 2 with the body and is offset by
// dx_entry from it, at the exit face theta_exit = (angle + angle_diff) / 2 and
// dx_exit. The body is shifted by -rbend_shift, and by half the sagitta if
// compensate_sagitta. The linear fringe is the expansion around the angle of
// the orbit at the face, hence the feed-down angles.
GPUFUN
void rbend_straight_body_patches(
    const double angle,
    const double h,
    const double length_straight,
    const double angle_diff,
    const double rbend_shift,
    const int64_t compensate_sagitta,
    MagnetFace entry_face,
    MagnetFace exit_face
) {
    const double theta_entry = (angle - angle_diff) / 2.0;
    const double theta_exit = (angle + angle_diff) / 2.0;
    double cos_theta_entry = 1., sin_theta_exit = 0., cos_theta_exit = 1.;
    if (fabs(theta_entry) > 1e-10) {
        cos_theta_entry = cos(theta_entry);
    }
    if (fabs(theta_exit) > 1e-10) {
        sin_theta_exit = sin(theta_exit);
        cos_theta_exit = cos(theta_exit);
    }

    double dx_mid = 0. - rbend_shift;
    if (compensate_sagitta && fabs(angle) > 1e-10) {
        // Half the sagitta, via 1 - cos(u) = 2*sin(u/2)^2 to avoid the
        // cancellation that the large 1/h would amplify.
        double const sin_rbqa = sin(angle / 4.);
        dx_mid += sin_rbqa * sin_rbqa / h;
    }

    double dx_entry = dx_mid;
    double dx_exit = dx_mid;
    if (fabs(angle) > 1e-10) {
        double const px0_entry = sin(theta_entry);
        double const px0_mid = px0_entry - h * length_straight / 2;
        double const sqrt_mid = sqrt(1 - px0_mid * px0_mid);
        // Rationalised forms of (1/h)*(sqrt_mid - cos_theta_entry) and
        // (1/h)*(cos_theta_exit - sqrt_mid), so that h cancels exactly.
        dx_entry -=
            0.5 * length_straight * (px0_entry + px0_mid) / (sqrt_mid + cos_theta_entry);
        double const term_exit_over_h =
            2 * cos(angle / 2.) * sin(angle_diff / 2.) / h + length_straight / 2.;
        dx_exit -=
            term_exit_over_h * (px0_mid + sin_theta_exit) / (cos_theta_exit + sqrt_mid);
    }

    MagnetFace_set_angle_fdown(entry_face, MagnetFace_get_angle_fdown(entry_face) + theta_entry);
    MagnetFace_set_angle_fdown(exit_face, MagnetFace_get_angle_fdown(exit_face) + theta_exit);
    MagnetFace_set_has_patch(entry_face, 1);
    MagnetFace_set_patch_angle(entry_face, theta_entry);
    MagnetFace_set_patch_dx(entry_face, dx_entry);
    MagnetFace_set_has_patch(exit_face, 1);
    MagnetFace_set_patch_angle(exit_face, theta_exit);
    MagnetFace_set_patch_dx(exit_face, dx_exit);
}


// Prepare the kick of a body (the per-element part of
// track_magnet_kick_single_particle) from the field and the split of the
// tracking model
GPUFUN
void MagnetKick_prepare(
    MagnetKick kick,
    const MagnetField field,
    const MagnetSplit* split,
    double length,
    double inv_factorial_order,
    double factor_knl_ksl,
    double hxl,
    uint8_t is_empty
) {
    const int64_t order = MagnetField_get_kl_order(field);
    const xt_float_or_tpsa* knl = MagnetField_knl(field);

    MagnetKick_set_length(kick, length);
    MagnetKick_set_order(kick, order);
    MagnetKick_set_inv_factorial_order(kick, inv_factorial_order);
    MagnetKick_set_factor_knl_ksl(kick, factor_knl_ksl);
    // The kick part of k0, k1 (configure_tracking_model) and the other
    // strengths, times the length
    XT_FT(MagnetKick_getp1_knl_main(kick, 0)) = split->k0_kick * length;
    XT_FT(MagnetKick_getp1_knl_main(kick, 1)) = split->k1_kick * length;
    XT_FT(MagnetKick_getp1_knl_main(kick, 2)) = XT_FT(MagnetField_getp1_knorm(field, 2)) * length;
    XT_FT(MagnetKick_getp1_knl_main(kick, 3)) = XT_FT(MagnetField_getp1_knorm(field, 3)) * length;
    for (int ii = 0; ii <= 3; ii++) {
        XT_FT(MagnetKick_getp1_ksl_main(kick, ii)) = XT_FT(MagnetField_getp1_kskew(field, ii)) * length;
    }
    const double h = split->h_kick;
    double htot = h;
    if (length != 0) {
        htot += hxl / length;
    }
    MagnetKick_set_h(kick, h);
    MagnetKick_set_hxl(kick, hxl);
    MagnetKick_set_htot(kick, htot);
    // k0h correction can be computed from this term in the hamiltonian
    // H = 1/2 h k0 x^2
    // (see MAD 8 physics manual, eq. 5.15, and apply Hamilton's eq. dp/ds = -dH/dx)
    xt_float_or_tpsa k0l_mult = 0.0;
    if (order >= 0) {
        k0l_mult = knl[0] * factor_knl_ksl;
    }
    XT_FT(MagnetKick_getp_k0_h_coeff(kick)) = split->k0_h_correction * length + k0l_mult;
    // k1h correction can be computed from this term in the hamiltonian
    // H = 1/3 hk1 x^3 - 1/2 hk1 xy^2
    // (see MAD 8 physics manual, eq. 5.15, and apply Hamilton's eq. dp/ds = -dH/dx)
    xt_float_or_tpsa k1l_mult = 0.0;
    if (order >= 1) {
        k1l_mult = knl[1] * factor_knl_ksl;
    }
    XT_FT(MagnetKick_getp_k1_h_coeff(kick)) = split->k1_h_correction * length + k1l_mult;
    MagnetKick_set_rot_frame(kick, split->kick_rot_frame);
    MagnetKick_set_do_drift_correction(kick, split->do_drift_correction);
    MagnetKick_set_is_empty(kick, is_empty);
    uint8_t has_main = 0;
    for (int ii = 0; ii <= 3; ii++) {
        has_main |= !xt_float_or_tpsa_is_zero(XT_FT(MagnetKick_getp1_knl_main(kick, ii)))
                    || !xt_float_or_tpsa_is_zero(XT_FT(MagnetKick_getp1_ksl_main(kick, ii)));
    }
    MagnetKick_set_has_main(kick, has_main);
    MagnetKick_set_has_curvature(kick, h != 0 || hxl != 0);
}


// Derive the integration of the body from its settings and the field of the
// element (the strengths of the faces, see track_magnet_edge.h; they need not
// be prepared).
GPUFUN XT_INLINE
void MagnetBody_prepare(
    MagnetBody body,
    const MagnetField field,
    int8_t model,
    int8_t integrator,
    int64_t num_multipole_kicks,
    double weight,  // fraction of the element (slices)
    double h,
    double hxl
) {
    const double body_length = MagnetBody_get_length(body);
    const int64_t order = MagnetField_get_kl_order(field);
    const xt_float_or_tpsa* knl = MagnetField_knl(field);
    const xt_float_or_tpsa* ksl = MagnetField_ksl(field);
    xt_float_or_tpsa_arg k0 = XT_FT(MagnetField_getp1_knorm(field, 0));
    xt_float_or_tpsa_arg k1 = XT_FT(MagnetField_getp1_knorm(field, 1));
    xt_float_or_tpsa_arg k2 = XT_FT(MagnetField_getp1_knorm(field, 2));
    xt_float_or_tpsa_arg k3 = XT_FT(MagnetField_getp1_knorm(field, 3));
    xt_float_or_tpsa_arg k0s = XT_FT(MagnetField_getp1_kskew(field, 0));
    xt_float_or_tpsa_arg k1s = XT_FT(MagnetField_getp1_kskew(field, 1));
    xt_float_or_tpsa_arg k2s = XT_FT(MagnetField_getp1_kskew(field, 2));
    xt_float_or_tpsa_arg k3s = XT_FT(MagnetField_getp1_kskew(field, 3));
    xt_float_or_tpsa_arg ks = XT_FT(MagnetField_getp_ks(field));

    if (model == MODEL_KICK_ONLY) {
        integrator = INTEGRATOR_UNIFORM;
        num_multipole_kicks = 1;
    }

    // Adjust the number of multipole kicks based on the weight
    if (weight != 1.0 && num_multipole_kicks > 0) {
        num_multipole_kicks = (int64_t)ceil(num_multipole_kicks * weight);
    }

    const MagnetSplit split = configure_tracking_model(model, k0, k1, h, ks);
    const MagnetDrift drift = MagnetBody_getp_drift(body);
    MagnetDrift_set_model(drift, split.drift_model);
    XT_FT(MagnetDrift_getp_k0(drift)) = split.drift_k0;
    XT_FT(MagnetDrift_getp_k1(drift)) = split.drift_k1;
    XT_FT(MagnetDrift_getp_ks(drift)) = split.drift_ks;
    MagnetDrift_set_h(drift, split.drift_h);
    MagnetDrift_set_x0_solenoid(drift, MagnetField_get_x0_solenoid(field));
    MagnetDrift_set_y0_solenoid(drift, MagnetField_get_y0_solenoid(field));

    // Compute the number of kicks for auto mode
    if (num_multipole_kicks == 0) {  // num_multipole_kicks = 0 means auto mode
        uint8_t const kick_inactive =
            kick_is_inactive(order, knl, ksl, k0, k1, k2, k3, k0s, k1s, k2s, k3s, h);

        // The exact matrix model needs a path-correction kick even at zero field.
        // Otherwise, guess the number of active kicks or perform only the drift.
        if (split.do_drift_correction || !kick_inactive) {
            if (fabs(h) < 1e-8) {
                num_multipole_kicks =
                    1;  // straight magnet, one multipole kick in the middle
            }
            else {
                double b_circum = 2 * 3.14159 / fabs(h);
                num_multipole_kicks = fabs(body_length) / b_circum
                                      / 0.5e-3;  // 0.5 mrad per kick (on average)
                if (num_multipole_kicks < 1) {
                    num_multipole_kicks = 1;
                }
            }
        }
    }

    const uint8_t kick_is_empty = (hxl == 0.0)
                          && kick_is_inactive(
                              order,
                              knl,
                              ksl,
                              split.k0_kick,
                              split.k1_kick,
                              k2,
                              k3,
                              k0s,
                              k1s,
                              k2s,
                              k3s,
                              split.h_kick
                          );

    MagnetBody_set_integrator(body, integrator);
    MagnetBody_set_num_kicks(body, num_multipole_kicks);

    // Curvature for the radiation (only for length > 0)
    if (fabs(split.drift_h) > 0) {
        MagnetBody_set_h_radiation(body, split.drift_h);
    }
    else {
        MagnetBody_set_h_radiation(
            body, body_length > 0 ? split.h_kick + hxl / body_length : 0.
        );
    }

    const uint8_t drift_only = num_multipole_kicks == 0
                               && xt_float_or_tpsa_const_part(split.k0_kick) == 0
                               && xt_float_or_tpsa_const_part(split.k1_kick) == 0
                               && split.h_kick == 0;
    // An empty kick still applies the path correction of mat-kick-mat-exact
    const uint8_t kick_is_noop = kick_is_empty && !split.do_drift_correction;
    MagnetBody_set_drift_only(body, drift_only);
    MagnetBody_set_kick_is_noop(body, kick_is_noop);
    if (drift_only || kick_is_noop) return;

    MagnetKick_prepare(
        MagnetBody_getp_kick(body),
        field,
        &split,
        body_length,
        MagnetBody_get_inv_factorial_order(body),
        MagnetBody_get_factor_knl_ksl(body),
        hxl,
        kick_is_empty
    );
}


// The model or integrator of an element: 0 selects the default of the
// element type (default_magnet_config.h)
GPUFUN XT_INLINE
int8_t magnet_or_default(int8_t value, int8_t default_value) {
    return value == 0 ? default_value : value;
}


// The radiation flag of a slice: 10 takes that of the parent element
GPUFUN XT_INLINE
int64_t magnet_radiation_flag(int64_t flag, int64_t parent_flag) {
    return flag == 10 ? parent_flag : flag;
}



// Settings of a magnet that has nothing beyond a straight body: no field,
// no multipoles, default model and integrator, no radiation, not an RBend.
// Fills the given struct in place (no copy: copying a TPSA allocates; the
// TPSA members are only assigned scalars here). Each element then sets the
// fields it has and lists the others under "Kept at default".
GPUFUN
void MagnetSettings_init_default(MagnetSettings* s) {
    s->weight = 1.;
    s->length = 0.;
    s->order = -1;
    s->inv_factorial_order = 1.;
    s->knl = NULL;
    s->ksl = NULL;
    s->order_rel = -1;
    s->knl_rel = NULL;
    s->ksl_rel = NULL;
    s->rel_ref_strength = 0.;
    s->num_multipole_kicks = 0;
    s->model = MODEL_ADAPTIVE;
    s->integrator = INTEGRATOR_ADAPTIVE;
    s->radiation_flag = 0;
    s->radiation_record = NULL;
    s->delta_taper = 0.;
    s->h = 0.;
    s->hxl = 0.;
    s->k0 = 0.;
    s->k1 = 0.;
    s->k2 = 0.;
    s->k3 = 0.;
    s->k0s = 0.;
    s->k1s = 0.;
    s->k2s = 0.;
    s->k3s = 0.;
    s->k_order = MAGNET_K_ORDER;
    s->ks = 0.;
    s->dks_ds = 0.;
    s->x0_solenoid = 0.;
    s->y0_solenoid = 0.;
    s->rbend_model = RBEND_MODEL_NONE;
    s->rbend_compensate_sagitta = 0;
    s->rbend_shift = 0.;
    s->rbend_angle_diff = 0.;
    s->length_straight = 0.;
    s->body_active = 1;
}


// A face of the element (is_exit_face) in the plan: the entry face is crossed
// first when tracking forward, last when backtracking
GPUFUN XT_INLINE
MagnetFace MagnetPlan_face(MagnetPlan plan, const uint8_t is_exit_face, const uint8_t backtrack) {
    return (is_exit_face == backtrack) ? MagnetPlan_getp_first_face(plan)
                                       : MagnetPlan_getp_second_face(plan);
}

// A face of the element, reset, for the element to set before MagnetPlan_prepare
GPUFUN XT_INLINE
MagnetFace MagnetPlan_new_face(MagnetPlan plan, const uint8_t is_exit_face, LocalParticle* part0) {
    const MagnetFace face = MagnetPlan_face(
        plan, is_exit_face, LocalParticle_check_track_flag(part0, XS_FLAG_BACKTRACK)
    );
    MagnetFace_reset(face, is_exit_face);
    return face;
}


// Prepare a magnet, or a slice of it (weight: its fraction of the length),
// for track_magnet_particles: its first face, its body and its second face, in
// the order of the direction of travel. The settings of the faces, in the
// element frame, are already in the plan (MagnetPlan_new_face). The misalignment is set by the caller
// (X_prepare_plan). Reads only the track flags of part0.
// Runs once per element (prepare_magnet_plans), not per particle: not
// inlined, so the kernel holds one copy, not one per magnet class.
GPUFUN XT_NOINLINE
void MagnetPlan_prepare(
    MagnetPlan plan,
    const MagnetSettings* settings,
    LocalParticle* part0
) {
    // Values changed below (RBend, backtracking, tapering) are local copies;
    // the rest is read from settings
    const double weight = settings->weight;
    double length = settings->length;
    double h = settings->h;
    double hxl = settings->hxl;
    const int64_t rbend_model = settings->rbend_model;
    const int64_t radiation_flag = settings->radiation_flag;

    const uint8_t backtrack = LocalParticle_check_track_flag(part0, XS_FLAG_BACKTRACK);

    MagnetTransform_set_active(MagnetPlan_getp_transform(plan), 0);
    MagnetPlan_set_is_thick(plan, 0);

    // Field of the body and the faces (the multipoles are added below)
    const MagnetField field = MagnetPlan_getp_field(plan);
    XT_FT(MagnetField_getp1_knorm(field, 0)) = settings->k0;
    XT_FT(MagnetField_getp1_knorm(field, 1)) = settings->k1;
    XT_FT(MagnetField_getp1_knorm(field, 2)) = settings->k2;
    XT_FT(MagnetField_getp1_knorm(field, 3)) = settings->k3;
    XT_FT(MagnetField_getp1_kskew(field, 0)) = settings->k0s;
    XT_FT(MagnetField_getp1_kskew(field, 1)) = settings->k1s;
    XT_FT(MagnetField_getp1_kskew(field, 2)) = settings->k2s;
    XT_FT(MagnetField_getp1_kskew(field, 3)) = settings->k3s;
    MagnetField_set_k_order(field, settings->k_order);
    XT_FT(MagnetField_getp_ks(field)) = settings->ks;
    MagnetField_set_x0_solenoid(field, settings->x0_solenoid);
    MagnetField_set_y0_solenoid(field, settings->y0_solenoid);

    // Faces, in the element frame; the RBend changes their angles or adds
    // patches
    const MagnetFace entry_face = MagnetPlan_face(plan, 0, backtrack);
    const MagnetFace exit_face = MagnetPlan_face(plan, 1, backtrack);

    // RBend: the faces of a curved body are at the angle of the reference
    // orbit; a straight body follows the chord, and its length and the
    // length of the curved reference path differ
    double length_curved = 0.;
    double const angle = h * length;
    if (rbend_model == RBEND_MODEL_CURVED_BODY){
        MagnetFace_set_angle(
            entry_face,
            MagnetFace_get_angle(entry_face) + (angle - settings->rbend_angle_diff) / 2.0
        );
        MagnetFace_set_angle(
            exit_face,
            MagnetFace_get_angle(exit_face) + (angle + settings->rbend_angle_diff) / 2.0
        );
    }
    else if (rbend_model == RBEND_MODEL_STRAIGHT_BODY){
        rbend_straight_body_patches(
            angle,
            h,
            settings->length_straight,
            settings->rbend_angle_diff,
            settings->rbend_shift,
            settings->rbend_compensate_sagitta,
            entry_face,
            exit_face
        );
        length_curved = length;
        length = settings->length_straight;
        h = 0; // treat magnet as straight
    }

    // Direction of travel: backtracking runs the body backwards; the faces
    // invert their own maps
    double body_length, body_length_curved, knl_scale_body;
    if (backtrack) {
        body_length = -length * weight;
        body_length_curved = -length_curved * weight;
        knl_scale_body = -weight;
        hxl = -hxl;
    }
    else {
        body_length = length * weight;
        body_length_curved = length_curved * weight;
        knl_scale_body = weight;
    }

    // Tapering: the factor of the element, 1 + delta_taper, is folded into the
    // strengths of the plan. While compensating the radiation energy loss the
    // plan holds the nominal strengths and the factor comes from the particle
    // (magnet_taper_factor)
    MagnetPlan_set_taper(plan, 1.);
    MagnetPlan_set_taper_from_delta(plan, 0);
    double taper_const = 1.;
#if !defined(XTRACK_TPSA_TRACK) && !defined(XTRACK_MULTIPOLE_NO_SYNRAD)
    if (radiation_flag) {
        if (LocalParticle_check_track_flag(part0, XS_FLAG_SR_TAPER)) {
            MagnetPlan_set_taper_from_delta(plan, 1);
        }
        else {
            taper_const = 1. + settings->delta_taper;
        }
    }
#endif
    if (taper_const != 1.) {
        for (int ii = 0; ii <= MAGNET_K_ORDER; ii++) {
            XT_FT(MagnetField_getp1_knorm(field, ii)) = taper_const * XT_FT(MagnetField_getp1_knorm(field, ii));
            XT_FT(MagnetField_getp1_kskew(field, ii)) = taper_const * XT_FT(MagnetField_getp1_kskew(field, ii));
        }
        knl_scale_body *= taper_const;
    }

    // Multipoles of the body and the faces
    const int64_t order_tot = magnet_multipole_order(
        settings->order, settings->order_rel, settings->knl_rel, settings->ksl_rel
    );
    double inv_factorial_order_tot = settings->inv_factorial_order;
    for (int64_t ii = settings->order + 1; ii <= order_tot; ii++) {
        inv_factorial_order_tot /= ii;
    }
#ifdef XTRACK_TPSA_TRACK
    // The multipoles are lifted to TPSAs, in the plan
    magnet_total_multipoles(
        settings->order,
        settings->knl,
        settings->ksl,
        settings->order_rel,
        settings->knl_rel,
        settings->ksl_rel,
        settings->rel_ref_strength,
        order_tot,
        XT_FT_PTR(MagnetPlan_getp1_knl_tot(plan, 0)),
        XT_FT_PTR(MagnetPlan_getp1_ksl_tot(plan, 0))
    );
    MagnetField_set_knl(field, MAGNET_PLAN_ADDR(XT_FT_PTR(MagnetPlan_getp1_knl_tot(plan, 0))));
    MagnetField_set_ksl(field, MAGNET_PLAN_ADDR(XT_FT_PTR(MagnetPlan_getp1_ksl_tot(plan, 0))));
#else
    // Without a relative component the element's own arrays are used
    MagnetField_set_knl(field, MAGNET_PLAN_ADDR(settings->knl));
    MagnetField_set_ksl(field, MAGNET_PLAN_ADDR(settings->ksl));
    if (magnet_has_rel(settings->order_rel, settings->knl_rel, settings->ksl_rel)) {
        magnet_total_multipoles(
            settings->order,
            settings->knl,
            settings->ksl,
            settings->order_rel,
            settings->knl_rel,
            settings->ksl_rel,
            settings->rel_ref_strength,
            order_tot,
            XT_FT_PTR(MagnetPlan_getp1_knl_tot(plan, 0)),
            XT_FT_PTR(MagnetPlan_getp1_ksl_tot(plan, 0))
        );
        MagnetField_set_knl(field, MAGNET_PLAN_ADDR(XT_FT_PTR(MagnetPlan_getp1_knl_tot(plan, 0))));
        MagnetField_set_ksl(field, MAGNET_PLAN_ADDR(XT_FT_PTR(MagnetPlan_getp1_ksl_tot(plan, 0))));
    }
#endif
    MagnetField_set_kl_order(field, order_tot);
    MagnetField_set_knl_scale(field, taper_const);
    MagnetField_set_length(field, length);

    // Faces, in the order they are crossed: the particle enters the field at
    // the first and leaves it at the second, in both directions of travel
    const MagnetFace first_face = MagnetPlan_getp_first_face(plan);
    const MagnetFace second_face = MagnetPlan_getp_second_face(plan);
    if (MagnetFace_uses_field(entry_face) || MagnetFace_uses_field(exit_face)) {
        MagnetField_prepare(field);
    }
    MagnetFace_prepare(first_face, field, backtrack);
    MagnetFace_prepare(second_face, field, backtrack);
    // A face that does nothing is skipped as a whole, not sub-map by sub-map
    MagnetPlan_set_first_active(plan, MagnetFace_get_has_work(first_face));
    MagnetPlan_set_second_active(plan, MagnetFace_get_has_work(second_face));

    // Faces only (edge slices)
    const uint8_t body_active = settings->body_active != 0;
    MagnetPlan_set_body_active(plan, body_active);
    if (!body_active) return;

    // Body
    const MagnetBody body = MagnetPlan_getp_body(plan);
    MagnetBody_set_length(body, body_length);
    MagnetBody_set_inv_factorial_order(body, inv_factorial_order_tot);
    MagnetBody_set_factor_knl_ksl(body, knl_scale_body);
    MagnetBody_set_dks_ds(body, backtrack ? -settings->dks_ds : settings->dks_ds);
    MagnetBody_set_radiation_flag(body, radiation_flag);
    MagnetBody_set_radiation_record(body, MAGNET_PLAN_ADDR(settings->radiation_record));
    MagnetBody_prepare(
        body,
        field,
        settings->model,
        settings->integrator,
        settings->num_multipole_kicks,
        weight,
        h,
        hxl
    );

    // Straight body: the path length is measured on the curved reference
    // (not for the kick only model)
    MagnetPlan_set_correct_s(
        plan, rbend_model == RBEND_MODEL_STRAIGHT_BODY && settings->model != MODEL_KICK_ONLY
    );
    MagnetPlan_set_ds_curved(plan, body_length_curved - body_length);
}


// Misalignment on one side of a magnet: before it (at_exit = 0) or after it
// (at_exit = 1), in the direction of travel. Returns 0 if the particles were
// lost (and the magnet must not be tracked).
GPUFUN
int8_t magnet_transform(LocalParticle* part0, const MagnetTransform tr, int8_t at_exit) {
    if (!at_exit && MagnetTransform_get_invalid(tr)) {
        START_PER_PARTICLE_BLOCK(part0, part);
            LocalParticle_set_state(part, XT_INVALID_THIN_SLICE_TRANSFORM);
        END_PER_PARTICLE_BLOCK;
        return 0;
    }
    int8_t const backtrack = LocalParticle_check_track_flag(part0, XS_FLAG_BACKTRACK);
    // Backtracking crosses the exit side first
    uint8_t const entry_side = at_exit == backtrack;
    double const shift_x = MagnetTransform_get_shift_x(tr);
    double const shift_y = MagnetTransform_get_shift_y(tr);
    double const shift_s = MagnetTransform_get_shift_s(tr);
    double const rot_x_rad = MagnetTransform_get_rot_x_rad(tr);
    double const rot_y_rad = MagnetTransform_get_rot_y_rad(tr);
    double const rot_s_rad = MagnetTransform_get_rot_s_rad(tr);
    double const rot_s_rad_no_frame = MagnetTransform_get_rot_s_rad_no_frame(tr);
    double const anchor = MagnetTransform_get_anchor(tr);
    double const length = MagnetTransform_get_length(tr);
    if (MagnetTransform_get_curved(tr)) {
        double const angle = MagnetTransform_get_angle(tr);
        double const h = MagnetTransform_get_h(tr);
        if (entry_side) {
            track_misalignment_entry_curved(
                part0, shift_x, shift_y, shift_s, rot_y_rad, rot_x_rad,
                rot_s_rad_no_frame, anchor, length, angle, h, rot_s_rad, backtrack
            );
        }
        else {
            track_misalignment_exit_curved(
                part0, shift_x, shift_y, shift_s, rot_y_rad, rot_x_rad,
                rot_s_rad_no_frame, anchor, length, angle, h, rot_s_rad, backtrack
            );
        }
    }
    else {
        if (entry_side) {
            track_misalignment_entry_straight(
                part0, shift_x, shift_y, shift_s, rot_y_rad, rot_x_rad,
                rot_s_rad_no_frame, anchor, length, rot_s_rad, backtrack
            );
        }
        else {
            track_misalignment_exit_straight(
                part0, shift_x, shift_y, shift_s, rot_y_rad, rot_x_rad,
                rot_s_rad_no_frame, anchor, length, rot_s_rad, backtrack
            );
        }
    }
    return 1;
}


// The factor of the strengths of a magnet for a particle: 1, or while
// compensating the radiation energy loss 1 + delta of the particle at the
// entrance of the magnet (the plan then holds the nominal strengths)
GPUFUN XT_INLINE
double magnet_taper_factor(LocalParticle* part, const MagnetPlan plan) {
#if !defined(XTRACK_TPSA_TRACK) && !defined(XTRACK_MULTIPOLE_NO_SYNRAD)
    if (MagnetPlan_get_taper_from_delta(plan)) return 1. + LocalParticle_get_delta(part);
#endif
    return MagnetPlan_get_taper(plan);
}


// Track a magnet, or a slice of it, from its plan: the misalignment, the
// first face, the body, the second face. Every magnet class is tracked by
// this one function, called from one place in track_line, where it is
// inlined (faster than a call on GPU, with NVRTC and clang).
GPUFUN XT_INLINE
void track_magnet_particles(LocalParticle* part0, const MagnetPlan plan) {
    const MagnetTransform transform = MagnetPlan_getp_transform(plan);
    const uint8_t transform_active = MagnetTransform_get_active(transform);
    if (transform_active && !magnet_transform(part0, transform, 0)) {
        return;
    }

    const MagnetField field = MagnetPlan_getp_field(plan);
    const MagnetFace first_face = MagnetPlan_getp_first_face(plan);
    const MagnetFace second_face = MagnetPlan_getp_second_face(plan);
    const uint8_t first_active = MagnetPlan_get_first_active(plan);
    const uint8_t second_active = MagnetPlan_get_second_active(plan);
    if (MagnetPlan_get_body_active(plan)) {
        const MagnetBody body = MagnetPlan_getp_body(plan);
        const uint8_t correct_s = MagnetPlan_get_correct_s(plan);
        const double ds_curved = MagnetPlan_get_ds_curved(plan);
        double dp_record_exit, dpx_record_exit, dpy_record_exit;
        START_PER_PARTICLE_BLOCK(part0, part)
            ;
            const double taper = magnet_taper_factor(part, plan);
            if (first_active) MagnetFace_enter_single_particle(part, first_face, field, taper);
            track_magnet_body_single_particle(
                part, body, field, taper, &dp_record_exit, &dpx_record_exit, &dpy_record_exit
            );
            if (correct_s) {
                LocalParticle_add_to_s(part, ds_curved);
                LocalParticle_add_to_zeta(part, ds_curved);
            }
            if (second_active) MagnetFace_leave_single_particle(part, second_face, field, taper);
        END_PER_PARTICLE_BLOCK;
    }
    else {
        // Faces only (edge slices)
        START_PER_PARTICLE_BLOCK(part0, part)
            ;
            const double taper = magnet_taper_factor(part, plan);
            if (first_active) MagnetFace_enter_single_particle(part, first_face, field, taper);
            if (second_active) MagnetFace_leave_single_particle(part, second_face, field, taper);
        END_PER_PARTICLE_BLOCK;
    }

    if (transform_active) {
        magnet_transform(part0, transform, 1);
    }
}

#ifdef XTRACK_TPSA_TRACK
// Call f on the address of every strength of a plan
template <class F>
void MagnetPlan_for_each_strength(MagnetPlan plan, F f) {
    const MagnetField field = MagnetPlan_getp_field(plan);
    const MagnetBody body = MagnetPlan_getp_body(plan);
    const MagnetDrift drift = MagnetBody_getp_drift(body);
    const MagnetKick kick = MagnetBody_getp_kick(body);
    for (int ii = 0; ii <= MAGNET_K_ORDER; ii++) {
        f(MagnetField_getp1_knorm(field, ii));
        f(MagnetField_getp1_kskew(field, ii));
        f(MagnetKick_getp1_knl_main(kick, ii));
        f(MagnetKick_getp1_ksl_main(kick, ii));
    }
    for (int ii = 0; ii < 2; ii++) {
        f(MagnetField_getp1_r21(field, ii));
        f(MagnetField_getp1_r43(field, ii));
    }
    f(MagnetField_getp_ks(field));
    f(MagnetField_getp_k0(field));
    f(MagnetField_getp_k1(field));
    f(MagnetField_getp_k1s(field));
    f(MagnetDrift_getp_k0(drift));
    f(MagnetDrift_getp_k1(drift));
    f(MagnetDrift_getp_ks(drift));
    f(MagnetKick_getp_k0_h_coeff(kick));
    f(MagnetKick_getp_k1_h_coeff(kick));
    for (int ii = 0; ii <= MAGNET_MAX_MULTIPOLE_ORDER; ii++) {
        f(MagnetPlan_getp1_knl_tot(plan, ii));
        f(MagnetPlan_getp1_ksl_tot(plan, ii));
    }
}

// A plan for the time of a TPSA track call, its TPSAs owned here: made after
// the TPSA default scope (xt_tpsa::default_scope), freed before it ends
struct MagnetPlanTpsaStorage {
    alignas(8) char data[MAGNET_PLAN_SIZE] = {};
    MagnetPlanTpsaStorage() {
        MagnetPlan_for_each_strength(plan(), [](double* p) { new (p) xt_float_or_tpsa(0.); });
    }
    ~MagnetPlanTpsaStorage() {
        MagnetPlan_for_each_strength(plan(), [](double* p) { XT_FT(p).~xt_float_or_tpsa(); });
    }
    MagnetPlanTpsaStorage(const MagnetPlanTpsaStorage&) = delete;
    MagnetPlanTpsaStorage& operator=(const MagnetPlanTpsaStorage&) = delete;
    MagnetPlan plan() { return (MagnetPlan)(void*)data; }
};
#endif

#endif
