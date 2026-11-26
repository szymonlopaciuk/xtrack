/* This file serves as a template for generating per-element tracking functions
   that handle misalignments and transformations. It expects the following macros
   to be defined (as appropriate) before inclusion:

   - ELEMENT_NAME: Name of the beam element class
   - ALLOW_ROT_AND_SHIFT: Whether to include transformations handling
   - CURVED: Whether the element has the `angle` and `h` parameters
   - IS_SLICE: Whether to get parameters from the parent element (if slice)
   - IS_THICK: Whether the element has a length
   - IS_THICK_DYNAMIC: Whether the element has an isthick parameter
   - THIN_SLICE_OF_CURVED_ELEMENT: Special handling for thin slices of curved
     elements where we disallow some transformations.
   - MAGNET_PLAN: The element is a magnet, prepared as a MagnetPlan by
     ELEMENT_NAME_prepare_magnet and tracked by track_magnet_particles.
*/

#ifndef ELEMENT_NAME
    #error "This file must be included with ELEMENT_NAME macro defined."
#endif

#include "xtrack/headers/track.h"
#include "xtrack/headers/particle_states.h"
#include "xtrack/beam_elements/elements_src/track_misalignments.h"

// Shorthand for `ElementClassNameData`
#define ELEMENT_DATA CONCAT(ELEMENT_NAME, Data)

// Retrieve parameters attached to the element object
#define ELEMENT_GET(el, param_name) CONCAT3(ELEMENT_DATA, _get_, param_name)(el)

#ifdef IS_SLICE
    // Retrieve parameters from the parent element
    #define PARENT_GET(el, param_name) CONCAT3(ELEMENT_DATA, _get__parent_, param_name)(el)
    // In case of a slice the parameters need to be retrieved from the parent...
    #define GET_PARAM(el, param_name) PARENT_GET(el, param_name)
    // ...except the weight, which lives in the slice itself
    #define GET_WEIGHT(el) CONCAT(ELEMENT_DATA, _get_weight)(el)
#else
    // Standard retrieval of parameters from the element
    #define GET_PARAM(el, param_name) ELEMENT_GET(el, param_name)
    // Since it's not a slice, weight is always 1.0
    #define GET_WEIGHT(_) 1.0
#endif

// The misalignment functions for curved and straight elements have slightly
// different signatures, so we abstract the differences here, as we will call
// them in multiple places.
#ifdef CURVED
    #define MISALIGN_FUNCTION_ARGS \
        part0, \
        shift_x, \
        shift_y, \
        shift_s, \
        rot_y_rad, \
        rot_x_rad, \
        rot_s_rad_no_frame, \
        anchor, \
        length * weight, \
        angle * weight, \
        h, \
        rot_s_rad, \
        backtrack
    #define TRACK_MISALIGN_ENTRY track_misalignment_entry_curved(MISALIGN_FUNCTION_ARGS)
    #define TRACK_MISALIGN_EXIT track_misalignment_exit_curved(MISALIGN_FUNCTION_ARGS)
#else
    #define MISALIGN_FUNCTION_ARGS \
        part0, \
        shift_x, \
        shift_y, \
        shift_s, \
        rot_y_rad, \
        rot_x_rad, \
        rot_s_rad_no_frame, \
        anchor, \
        length * weight, \
        rot_s_rad, \
        backtrack
    #define TRACK_MISALIGN_ENTRY track_misalignment_entry_straight(MISALIGN_FUNCTION_ARGS)
    #define TRACK_MISALIGN_EXIT track_misalignment_exit_straight(MISALIGN_FUNCTION_ARGS)
#endif


#ifdef MAGNET_PLAN

/* The plan of a magnet: the magnet (ELEMENT_NAME_prepare_magnet) and its
   misalignment, applied by track_magnet_particles. Not inlined: the classes
   prepared in one switch (prepare_magnet_plan) would add up their frames.
*/
GPUFUN XT_NOINLINE
void CONCAT(ELEMENT_NAME, _prepare_plan)(
    ELEMENT_DATA el,
    MagnetPlan plan,
    LocalParticle* part0
) {
    CONCAT(ELEMENT_NAME, _prepare_magnet)(el, plan, part0);

    // Statically thick classes only: the global aperture is checked after them
    #ifdef IS_THICK
        MagnetPlan_set_is_thick(plan, 1);
    #endif

    #ifdef ALLOW_ROT_AND_SHIFT
        const MagnetTransform tr = MagnetPlan_getp_transform(plan);
        const double shift_x = GET_PARAM(el, shift_x);
        const double shift_y = GET_PARAM(el, shift_y);
        const double shift_s = GET_PARAM(el, shift_s);
        const double rot_x_rad = GET_PARAM(el, rot_x_rad);
        const double rot_y_rad = GET_PARAM(el, rot_y_rad);
        const double rot_s_rad = GET_PARAM(el, rot_s_rad);
        const double rot_s_rad_no_frame = GET_PARAM(el, rot_s_rad_no_frame);
        const uint8_t active = NONZERO(shift_x) || NONZERO(shift_y) || NONZERO(shift_s)
                               || NONZERO(rot_x_rad) || NONZERO(rot_y_rad)
                               || NONZERO(rot_s_rad) || NONZERO(rot_s_rad_no_frame);
        MagnetTransform_set_shift_x(tr, shift_x);
        MagnetTransform_set_shift_y(tr, shift_y);
        MagnetTransform_set_shift_s(tr, shift_s);
        MagnetTransform_set_rot_x_rad(tr, rot_x_rad);
        MagnetTransform_set_rot_y_rad(tr, rot_y_rad);
        MagnetTransform_set_rot_s_rad(tr, rot_s_rad);
        MagnetTransform_set_rot_s_rad_no_frame(tr, rot_s_rad_no_frame);
        MagnetTransform_set_active(tr, active);
        if (!active) return;

        const double weight = GET_WEIGHT(el);
        #if defined(IS_THICK_DYNAMIC)
            const double length = GET_PARAM(el, isthick) ? GET_PARAM(el, length) : 0.0;
        #elif defined(IS_THICK)
            const double length = GET_PARAM(el, length);
        #else
            const double length = 0.0;
        #endif
        MagnetTransform_set_length(tr, length * weight);
        #ifdef CURVED
            MagnetTransform_set_curved(tr, 1);
            MagnetTransform_set_angle(tr, GET_PARAM(el, angle) * weight);
            MagnetTransform_set_h(tr, GET_PARAM(el, h));
        #else
            MagnetTransform_set_curved(tr, 0);
            MagnetTransform_set_angle(tr, 0.);
            MagnetTransform_set_h(tr, 0.);
        #endif

        double anchor = GET_PARAM(el, rot_shift_anchor);
        #ifdef IS_SLICE
            anchor -= ELEMENT_GET(el, slice_offset);
        #endif
        MagnetTransform_set_anchor(tr, anchor);

        uint8_t invalid = 0;
        #ifdef THIN_SLICE_OF_CURVED_ELEMENT
            invalid = rot_x_rad != 0.0 || rot_y_rad != 0.0 || rot_s_rad_no_frame != 0.0;
        #endif
        MagnetTransform_set_invalid(tr, invalid);
    #endif
}

#ifndef XTRACK_TPSA_TRACK
/* Tracking of the element on its own (element.track()): its plan is prepared
   by ELEMENT_NAME_prepare_plan_kernel, in a plan of the element, then the
   particles are tracked through it.
*/
GPUFUN
void CONCAT(ELEMENT_NAME, _track_local_particle_with_transformations)(
    ELEMENT_DATA el,
    LocalParticle* part0,
    MagnetPlan plan
) {
    track_magnet_particles(part0, plan);
}

/* Prepare the plan of the element for element.track(): one thread (the
   track flags of element.track() are those of forward tracking). */
GPUKERN
void CONCAT(ELEMENT_NAME, _prepare_plan_kernel)(
    ELEMENT_DATA el,
    MagnetPlan plan,
    GPUGLMEM int8_t* io_buffer
) {
    #if defined(XO_CONTEXT_CUDA)
        if (blockDim.x * blockIdx.x + threadIdx.x != 0) return;
    #elif defined(XO_CONTEXT_CL)
        if (get_global_id(0) != 0) return;
    #endif
    LocalParticle lpart;
    lpart.io_buffer = io_buffer;
    lpart.track_flags = 0;
    CONCAT(ELEMENT_NAME, _prepare_plan)(el, plan, &lpart);
}
#else
/* TPSA tracking of the element on its own: its plan is made here, its TPSAs
   owned by C (MagnetPlanTpsaStorage)
*/
GPUFUN
void CONCAT(ELEMENT_NAME, _track_local_particle_with_transformations)(
    ELEMENT_DATA el,
    LocalParticle* part0
) {
    MagnetPlanTpsaStorage storage;
    const MagnetPlan plan = storage.plan();
    CONCAT(ELEMENT_NAME, _prepare_plan)(el, plan, part0);
    track_magnet_particles(part0, plan);
}
#endif

#else /* MAGNET_PLAN */

/* If transformations are enabled for an element, this applies them on one
   side of it: before the element (at_exit = 0) or after it (at_exit = 1), in
   the direction of travel. It assumes that the transformations are
   (1) possible for the element and (2) non-zero. Returns 0 if the particles
   were lost (and the element must not be tracked).
*/
#ifdef ALLOW_ROT_AND_SHIFT

GPUFUN
#ifdef XO_CONTEXT_CL
    __attribute__((noinline))
#endif
int8_t CONCAT(ELEMENT_NAME, _apply_transformations)(
    ELEMENT_DATA el,
    LocalParticle* part0,
    int8_t at_exit
) {
    // Retrieve the weight
    const double weight = GET_WEIGHT(el);

    // Determine the length
    #if defined(IS_THICK_DYNAMIC)
        const int8_t is_thick = GET_PARAM(el, isthick);
        const double length = is_thick ? GET_PARAM(el, length) : 0.0;
    #elif defined(IS_THICK)
        const double length = GET_PARAM(el, length);
    #else
        const double length = 0.0;
    #endif

    // Retrieve the `angle` and `h` if available
    #ifdef CURVED
        double angle = GET_PARAM(el, angle);
        double h = GET_PARAM(el, h);
    #endif

    // Retrieve misalignment parameters
    double const shift_x = GET_PARAM(el, shift_x);
    double const shift_y = GET_PARAM(el, shift_y);
    double const shift_s = GET_PARAM(el, shift_s);
    double const rot_x_rad = GET_PARAM(el, rot_x_rad);
    double const rot_y_rad = GET_PARAM(el, rot_y_rad);
    double const rot_s_rad = GET_PARAM(el, rot_s_rad);
    double const rot_s_rad_no_frame = GET_PARAM(el, rot_s_rad_no_frame);
    double anchor = GET_PARAM(el, rot_shift_anchor);

    #ifdef THIN_SLICE_OF_CURVED_ELEMENT
        if (!at_exit && (rot_x_rad != 0.0 || rot_y_rad != 0.0 || rot_s_rad_no_frame != 0.0)) {
            START_PER_PARTICLE_BLOCK(part0, part);
                LocalParticle_set_state(part, XT_INVALID_THIN_SLICE_TRANSFORM);
            END_PER_PARTICLE_BLOCK;
            return 0;
        }
    #endif

    int8_t const backtrack = LocalParticle_check_track_flag(part0, XS_FLAG_BACKTRACK);

    #ifdef IS_SLICE
        double const slice_offset = ELEMENT_GET(el, slice_offset);
        anchor = anchor - slice_offset;
    #endif

    // Backtracking crosses the exit side first
    if (at_exit == backtrack) {
        TRACK_MISALIGN_ENTRY;
    } else {
        TRACK_MISALIGN_EXIT;
    }
    return 1;
}

#endif /* ALLOW_ROT_AND_SHIFT */

/* This is the main function defined in this file. It tracks the element in
   its local frame, with the transformations before and after it if they are
   enabled and non-zero. The element is called from one place only, so that
   its code appears once in the kernel.
*/
GPUFUN
#ifdef XO_CONTEXT_CL
    __attribute__((noinline))
#endif
void CONCAT(ELEMENT_NAME, _track_local_particle_with_transformations)(
    ELEMENT_DATA el,
    LocalParticle* part0
) {
    #ifdef ALLOW_ROT_AND_SHIFT
        double const shift_x = GET_PARAM(el, shift_x);
        double const shift_y = GET_PARAM(el, shift_y);
        double const shift_s = GET_PARAM(el, shift_s);
        double const rot_x_rad = GET_PARAM(el, rot_x_rad);
        double const rot_y_rad = GET_PARAM(el, rot_y_rad);
        double const rot_s_rad = GET_PARAM(el, rot_s_rad);
        double const rot_s_rad_no_frame = GET_PARAM(el, rot_s_rad_no_frame);

        const int64_t rot_shift_active =
            NONZERO(shift_x) ||
            NONZERO(shift_y) ||
            NONZERO(shift_s) ||
            NONZERO(rot_x_rad) ||
            NONZERO(rot_y_rad) ||
            NONZERO(rot_s_rad) ||
            NONZERO(rot_s_rad_no_frame);

        if (rot_shift_active && !CONCAT(ELEMENT_NAME, _apply_transformations)(el, part0, 0)) {
            return;
        }
    #endif

    CONCAT(ELEMENT_NAME, _track_local_particle)(el, part0);

    #ifdef ALLOW_ROT_AND_SHIFT
        if (rot_shift_active) {
            CONCAT(ELEMENT_NAME, _apply_transformations)(el, part0, 1);
        }
    #endif
}


#endif /* MAGNET_PLAN */


#undef ELEMENT_DATA
#undef ELEMENT_GET
#undef PARENT_GET
#undef GET_PARAM
#undef GET_WEIGHT
#undef MISALIGN_FUNCTION_ARGS
#undef TRACK_MISALIGN_ENTRY
#undef TRACK_MISALIGN_EXIT
