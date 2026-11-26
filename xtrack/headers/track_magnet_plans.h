// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2025.                 //
// ######################################### //
#ifndef XTRACK_TRACK_MAGNET_PLANS_H
#define XTRACK_TRACK_MAGNET_PLANS_H

/* Magnet plans: the C side of xtrack/magnet_plans.py.

   A magnet is tracked in two steps. Once per track call, the kernel
   prepare_magnet_plans reads the attributes of each magnet and stores them,
   ready to use, in its MagnetPlan (an xobjects struct). The track kernel then tracks every
   particle through the plan with track_magnet_particles.

   This header is inlined in the track kernel source right after the switch
   prepare_magnet_plan, which Python generates because it depends on the type
   ids of the build. A build with no magnets has neither.
   In TPSA tracking, only the switch is used: the plans are made at the
   magnets (magnet_plans.py).
*/

#ifndef XTRACK_TPSA_TRACK

/* Fill the plans in the slots [first_plan, first_plan + num_plans), one
   thread per plan. magnet_plan_elements[slot] is the index of the magnet of
   the slot in the line.
*/
/*gpukern*/
void prepare_magnet_plans(
    /*gpuglmem*/ int8_t* buffer,
    ElementRefData elem_ref_data,
    ArrNMagnetPlan magnet_plans,
    /*gpuglmem*/ int64_t* magnet_plan_elements,
    int64_t first_plan,
    int64_t num_plans,
    /*gpuglmem*/ int8_t* io_buffer,
    uint64_t track_flags
) {
    // A plan reads only the track flags of the particle (the direction of
    // travel)
    LocalParticle lpart;
    lpart.io_buffer = io_buffer;
    lpart.track_flags = track_flags;

    // clang-format off: VECTORIZE_OVER is expanded by xobjects, as written
    VECTORIZE_OVER(ii, num_plans);
        const int64_t slot = first_plan + ii;
        const int64_t elem_idx = magnet_plan_elements[slot];
        /*gpuglmem*/ void* el = ElementRefData_member_elements(elem_ref_data, elem_idx);
        const int64_t elem_type = ElementRefData_typeid_elements(elem_ref_data, elem_idx);
        prepare_magnet_plan(el, elem_type, ArrNMagnetPlan_getp1(magnet_plans, slot), &lpart);
    END_VECTORIZE;
    // clang-format on
}

#endif

#endif /* XTRACK_TRACK_MAGNET_PLANS_H */
