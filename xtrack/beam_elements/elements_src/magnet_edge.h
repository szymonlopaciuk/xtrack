// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2025.                 //
// ######################################### //

#ifndef XTRACK_MAGNET_EDGE_H
#define XTRACK_MAGNET_EDGE_H

#include "xtrack/headers/track.h"
#include "xtrack/beam_elements/elements_src/track_magnet.h"

// The face as a magnet of faces only, as a plan for track_magnet_particles
// (the misalignment is added by MagnetEdge_prepare_plan,
// track_local_particle_with_transformations.h): the field of the face and
// the face, at the entry or the exit, without a body
GPUFUN
void MagnetEdge_prepare_magnet(
        MagnetEdgeData el,
        MagnetPlan plan,
        LocalParticle* part0
) {
    const uint8_t is_exit = MagnetEdgeData_get_is_exit(el) != 0;
    const int64_t k_order = MagnetEdgeData_get_k_order(el);
    GPUGLMEM const double* kn = MagnetEdgeData_getp1_kn(el, 0);
    GPUGLMEM const double* ks = MagnetEdgeData_getp1_ks(el, 0);

    // The other face does nothing (reset)
    MagnetFace entry_face = MagnetPlan_new_face(plan, 0, part0);
    MagnetFace exit_face = MagnetPlan_new_face(plan, 1, part0);
    MagnetFace face = is_exit ? exit_face : entry_face;
    MagnetFace_set_flags(face, MagnetEdgeData_get__flags(el));
    MagnetFace_set_angle(face, MagnetEdgeData_get_face_angle(el));
    MagnetFace_set_angle_fdown(face, MagnetEdgeData_get_face_angle_feed_down(el));
    MagnetFace_set_fint(face, MagnetEdgeData_get_fringe_integral(el));
    MagnetFace_set_hgap(face, MagnetEdgeData_get_half_gap(el));
    MagnetFace_set_f1(face, MagnetEdgeData_get_qint1(el));
    MagnetFace_set_f2(face, MagnetEdgeData_get_qint2(el));

    // Settings aligned by hand: not formatted
    // clang-format off
    MagnetSettings magnet;
    MagnetSettings_init_default(&magnet);
    magnet.length                   = MagnetEdgeData_get_length(el);
    magnet.order                    = MagnetEdgeData_get_kl_order(el);
    magnet.knl                      = MagnetEdgeData_getp1_knl(el, 0);
    magnet.ksl                      = MagnetEdgeData_getp1_ksl(el, 0);
    magnet.k_order                  = k_order <= MAGNET_K_ORDER ? k_order : MAGNET_K_ORDER;
    if (k_order >= 0) { magnet.k0 = kn[0]; magnet.k0s = ks[0]; }
    if (k_order >= 1) { magnet.k1 = kn[1]; magnet.k1s = ks[1]; }
    if (k_order >= 2) { magnet.k2 = kn[2]; magnet.k2s = ks[2]; }
    if (k_order >= 3) { magnet.k3 = kn[3]; magnet.k3s = ks[3]; }
    magnet.body_active              = 0;  // faces only
    // ks: not supported yet in MagnetEdge
    // clang-format on
    MagnetPlan_prepare(plan, &magnet, part0);
}

#endif // XTRACK_MAGNET_EDGE_H
