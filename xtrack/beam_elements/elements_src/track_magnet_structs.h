// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2025.                 //
// ######################################### //

#ifndef XTRACK_TRACK_MAGNET_STRUCTS_H
#define XTRACK_TRACK_MAGNET_STRUCTS_H

#include "xtrack/headers/track.h"
#include "xtrack/beam_elements/elements_src/edge_flags.h"

// The radiation record (SynchrotronRadiationRecordData, a pointer to this
// struct) is typedef'd by xobjects only in kernels with radiation; the
// struct tag works in all, including the drift-only users of this file.
typedef GPUGLMEM struct SynchrotronRadiationRecordData_s* MagnetRadiationRecord;

// Highest multipole order of a magnet (as in MAD-X), checked in Python
#define MAGNET_MAX_MULTIPOLE_ORDER 20

// The plan (MagnetPlan and its parts) is made of xobjects structs
// (magnet_plans.py). Here: the input of a plan and a result of
// configure_tracking_model.


// Settings of a magnet (or of a slice of it) for track_magnet_particles,
// filled by the element; the faces are set in the plan (MagnetPlan_new_face).
typedef struct MagnetSettings_s {
    double weight;
    double length;
    int64_t order;
    double inv_factorial_order;
    GPUGLMEM const double* knl;
    GPUGLMEM const double* ksl;
    int64_t order_rel;
    GPUGLMEM const double* knl_rel;
    GPUGLMEM const double* ksl_rel;
    xt_float_or_tpsa rel_ref_strength;  // knl_rel, ksl_rel: relative to it
    int64_t num_multipole_kicks;
    int8_t model;
    int8_t integrator;
    int64_t radiation_flag;
    MagnetRadiationRecord radiation_record;
    double delta_taper;
    double h;
    double hxl;
    xt_float_or_tpsa k0;
    xt_float_or_tpsa k1;
    xt_float_or_tpsa k2;
    xt_float_or_tpsa k3;
    xt_float_or_tpsa k0s;
    xt_float_or_tpsa k1s;
    xt_float_or_tpsa k2s;
    xt_float_or_tpsa k3s;
    int64_t k_order;  // order of k0 ... k3 used (MAGNET_K_ORDER, or lower)
    xt_float_or_tpsa ks;
    double dks_ds;
    double x0_solenoid;
    double y0_solenoid;
    int64_t rbend_model;  // RBEND_MODEL_*
    int64_t rbend_compensate_sagitta;
    double rbend_shift;
    double rbend_angle_diff;
    double length_straight;
    int64_t body_active;
#ifdef XTRACK_TPSA_TRACK
    // The TPSA type has no default constructor
    MagnetSettings_s()
        : rel_ref_strength(0.), k0(0.), k1(0.), k2(0.), k3(0.), k0s(0.), k1s(0.), k2s(0.),
          k3s(0.), ks(0.) {}
#endif
} MagnetSettings;


// How a magnet body is integrated: the split of the field between the drift
// and the kicks of the integrator, for a tracking model
typedef struct MagnetSplit_s {
    int8_t drift_model;
    xt_float_or_tpsa drift_k0, drift_k1, drift_ks;
    double drift_h;
    xt_float_or_tpsa k0_kick, k1_kick;
    double h_kick;
    xt_float_or_tpsa k0_h_correction, k1_h_correction;  // curvature terms of the kick
    int8_t kick_rot_frame;
    uint8_t do_drift_correction;  // exact path correction in the kick (model 9)
#ifdef XTRACK_TPSA_TRACK
    // The TPSA type has no default constructor
    MagnetSplit_s()
        : drift_k0(0.), drift_k1(0.), drift_ks(0.), k0_kick(0.), k1_kick(0.),
          k0_h_correction(0.), k1_h_correction(0.) {}
#endif
} MagnetSplit;

#endif // XTRACK_TRACK_MAGNET_STRUCTS_H
