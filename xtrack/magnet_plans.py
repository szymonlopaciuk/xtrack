# copyright ############################### #
# This file is part of the Xtrack Package.  #
# Copyright (c) CERN, 2025.                 #
# ######################################### #

"""Magnet plans: how the tracker tracks magnets.

What a plan is
--------------
A magnet (Bend, RBend, Quadrupole, ..., and their slices) is tracked in two
steps, both in C:

1. ``<Class>_prepare_plan`` reads the attributes of the element and stores
   them, ready to use, in a ``MagnetPlan`` (an xobjects struct, below): the
   field, the two faces in the order the particles cross them, the body and
   the transforms.
2. ``track_magnet_particles`` (track_magnet.h) tracks the particles through
   the plan. It does not read the element.
   ``element.track()`` prepares the plan of the element with a kernel of one
   thread, then tracks it.

Step 1 runs once per magnet and track call, step 2 once per particle. A class
is a magnet if it sets ``_magnet_plan = True``.

When things happen
------------------
- When the track kernel is built (``Tracker._compile_kernel_from_classes``):
  ``plan_source`` adds the C of the plans to the kernel source, and
  ``track_line_cases`` adds the magnet cases to the element switch of
  ``track_line``. ``attach_prepare_kernel`` keeps the kernel that fills the
  plans with the track kernel of the same build.
- At the first track call of a tracker (``MagnetPlans``): one plan is
  allocated per magnet of the line. Only the positions of the magnets in the
  line are used.
- At every track call (``MagnetPlans.prepare``): the kernel
  ``prepare_magnet_plans`` fills the plans of the magnets the call will cross,
  from the current values of their attributes. A change of an attribute
  between two calls is therefore always used: nothing is cached in Python.
- During tracking (``track_line``): each magnet is tracked from its slot.

Inspecting a plan
-----------------
``Tracker.get_magnet_plan(name)`` prepares the plan of one magnet and returns
it as a nested dict. In TPSA tracking the plan holds TPSAs, owned by C
(``MagnetPlanTpsaStorage``), and is not inspected.

Exceptions
----------
- Tapering: the factor of the element, ``1 + delta_taper``, is folded into
  the strengths of the plan. While compensating the radiation energy loss
  the plan holds the nominal strengths and the factor comes from the energy
  of each particle (``magnet_taper_factor``, track_magnet.h).
- TPSA tracking: the plan holds TPSA values, owned by C, so it is made at the
  magnet, in one plan per track call, in the same branch of ``track_line``.
"""

import numpy as np
import xobjects as xo

from .general import _pkg_root

_PLANS_HEADER = _pkg_root.joinpath("headers/track_magnet_plans.h")

# The plan, as xobjects structs (the C code reads them with the getters of
# xobjects). A field marked FT holds a strength: a double, or in TPSA builds
# the TPSA held in place (XT_FT, track_magnet_edge.h). Arrays the plan refers
# to (multipoles, radiation record) are held by address.
_FT = xo.Float64
_ADDR = xo.Int64


class MagnetField(xo.Struct):
    knorm = _FT[4]  # MAGNET_K_ORDER + 1
    kskew = _FT[4]
    k_order = xo.Int64
    knl = _ADDR
    ksl = _ADDR
    kl_order = xo.Int64
    knl_scale = xo.Float64
    length = xo.Float64
    ks = _FT
    x0_solenoid = xo.Float64
    y0_solenoid = xo.Float64
    # Derived when a face is active (MagnetField_prepare, MagnetFace_prepare)
    k0 = _FT
    k1 = _FT
    k1s = _FT
    r21 = _FT[2]  # linear edge, [is_exit_face]
    r43 = _FT[2]
    has_mult = xo.UInt8


class MagnetFace(xo.Struct):
    # Settings, in the element frame
    flags = xo.Int64  # edge_flags.h
    is_exit_face = xo.UInt8
    angle = xo.Float64
    angle_fdown = xo.Float64
    fint = xo.Float64
    hgap = xo.Float64
    f1 = xo.Float64
    f2 = xo.Float64
    has_ks = xo.UInt8
    ks = xo.Float64
    has_patch = xo.UInt8
    patch_angle = xo.Float64
    patch_dx = xo.Float64
    # Derived for the direction of travel (MagnetFace_prepare)
    backtrack = xo.UInt8
    has_work = xo.UInt8
    do_patch = xo.UInt8
    do_vector_potential = xo.UInt8
    do_linear = xo.UInt8
    do_rotate = xo.UInt8
    do_bend = xo.UInt8
    do_mult = xo.UInt8
    do_quad_lw = xo.UInt8
    do_quad_wedge = xo.UInt8
    do_wedge = xo.UInt8
    patch_sin = xo.Float64
    patch_cos = xo.Float64
    patch_tan = xo.Float64
    rot_sin = xo.Float64
    rot_cos = xo.Float64
    rot_tan = xo.Float64
    wedge_angle = xo.Float64
    quad_lw_g1 = xo.Float64
    quad_lw_g2 = xo.Float64


class MagnetDrift(xo.Struct):
    model = xo.Int8
    k0 = _FT
    k1 = _FT
    ks = _FT
    h = xo.Float64
    x0_solenoid = xo.Float64
    y0_solenoid = xo.Float64


class MagnetKick(xo.Struct):
    length = xo.Float64
    order = xo.Int64
    inv_factorial_order = xo.Float64
    factor_knl_ksl = xo.Float64
    knl_main = _FT[4]
    ksl_main = _FT[4]
    h = xo.Float64
    hxl = xo.Float64
    htot = xo.Float64
    k0_h_coeff = _FT
    k1_h_coeff = _FT
    rot_frame = xo.UInt8
    do_drift_correction = xo.UInt8
    is_empty = xo.UInt8
    has_main = xo.UInt8
    has_curvature = xo.UInt8


class MagnetBody(xo.Struct):
    length = xo.Float64
    integrator = xo.Int8
    num_kicks = xo.Int64
    drift_only = xo.UInt8
    kick_is_noop = xo.UInt8
    drift = MagnetDrift
    kick = MagnetKick
    inv_factorial_order = xo.Float64
    factor_knl_ksl = xo.Float64
    dks_ds = xo.Float64
    h_radiation = xo.Float64
    radiation_flag = xo.Int64
    radiation_record = _ADDR


class MagnetTransform(xo.Struct):
    active = xo.UInt8
    invalid = xo.UInt8
    curved = xo.UInt8
    shift_x = xo.Float64
    shift_y = xo.Float64
    shift_s = xo.Float64
    rot_x_rad = xo.Float64
    rot_y_rad = xo.Float64
    rot_s_rad = xo.Float64
    rot_s_rad_no_frame = xo.Float64
    anchor = xo.Float64
    length = xo.Float64
    angle = xo.Float64
    h = xo.Float64


class MagnetPlan(xo.Struct):
    field = MagnetField
    first_face = MagnetFace  # in the order they are crossed
    second_face = MagnetFace
    body = MagnetBody
    first_active = xo.UInt8
    second_active = xo.UInt8
    body_active = xo.UInt8
    correct_s = xo.UInt8
    ds_curved = xo.Float64
    is_thick = xo.UInt8
    taper_from_delta = xo.UInt8
    taper = xo.Float64
    transform = MagnetTransform
    knl_tot = _FT[21]  # MAGNET_MAX_MULTIPOLE_ORDER + 1
    ksl_tot = _FT[21]


# The size of a plan, for the TPSA plans made in C (MagnetPlanTpsaStorage)
MagnetPlan._extra_c_sources = [f'#define MAGNET_PLAN_SIZE ({MagnetPlan._size})']

MagnetPlanArray = MagnetPlan[:]


def magnet_classes(kernel_element_classes):
    """Return the magnets among the element classes of a kernel.

    Returns a list of ``(type_id, class)``, where ``type_id`` is the index of
    the class in ``kernel_element_classes``, i.e. its type id in the kernel.
    The classes are xobjects structs or, from a prebuilt kernel, the elements.
    """
    return [(ii, cc) for ii, cc in enumerate(kernel_element_classes)
            if vars(getattr(cc, '_DressingClass', cc)).get('_magnet_plan', False)]


def plan_source(magnet_classes):
    """Return the C source of the plans, for the track kernel.

    It is made of two parts:

    - the function ``prepare_magnet_plan``, a switch from the type id of an
      element to the ``<Class>_prepare_plan`` of its class. It is generated
      here, because the type ids depend on the classes of the build;
    - the header headers/track_magnet_plans.h, with the kernel
      ``prepare_magnet_plans``. The header is inlined rather than included, so
      that xobjects specializes its kernel for the context like the rest of
      the source.

    Without magnets, there is no plan source and no prepare kernel.
    """
    if not magnet_classes:
        return ""
    header = _PLANS_HEADER.read_text()

    cases = "\n".join(
        f"""    case {ii}:
        {cc.__name__.replace("Data", "")}_prepare_plan(({cc.__name__})el, plan, part);
        break;"""
        for ii, cc in magnet_classes)
    return f"""
// Fill the plan of the element el, of type id elem_type
GPUFUN
void prepare_magnet_plan(
    /*gpuglmem*/ void* el,
    int64_t elem_type,
    MagnetPlan plan,
    LocalParticle* part
) {{
    switch (elem_type) {{
{cases}
    }}
}}
""" + header


def track_line_cases(magnet_ids):
    """Return the cases of the element switch of ``track_line`` for the
    magnets, given their type ids.

    All magnets take the same branch: get the plan, then track the particles
    through it. Only where the plan comes from differs: its slot, or, in TPSA
    tracking, a plan made at the magnet. The global aperture is checked after
    the statically thick classes (``is_thick``), as in the other
    elements. The indentation matches the generated ``track_line``.
    """
    if not magnet_ids:
        return ""
    case_labels = "\n".join(
        f"                        case {ii}:" for ii in magnet_ids)
    return f"""
{case_labels}
                        {{
#ifdef XTRACK_TPSA_TRACK
                            // The plan holds TPSA values: it is prepared here, in the plan of the call
                            const MagnetPlan plan = magnet_plan_tpsa.plan();
                            prepare_magnet_plan(el, elem_type, plan, &lpart);
#else
                            // The plan prepared for this track call
                            const MagnetPlan plan = ArrNMagnetPlan_getp1(
                                magnet_plans, magnet_plan_slots[elem_idx]);
#endif
                            track_magnet_particles(&lpart, plan);
                            #ifdef XTRACK_GLOBAL_XY_LIMIT
                            if (MagnetPlan_get_is_thick(plan)) global_aperture_check(&lpart);
                            #endif
                            break;
                        }}
"""


def static_thick(cls):
    """Whether the class declares ``isthick = True`` literally (not a property)."""
    return vars(cls._DressingClass).get('_isthick') is True


def track_line_args():
    """Return the arguments that the plans add to ``track_line``: the plans,
    and the slot of each element."""
    return [
        xo.Arg(MagnetPlanArray, name="magnet_plans"),
        xo.Arg(xo.Int64, pointer=True, name="magnet_plan_slots"),
    ]


def prepare_kernel_description(tracker_data_type):
    """Return the description of the kernel ``prepare_magnet_plans``
    (headers/track_magnet_plans.h), which fills the plans."""
    return xo.Kernel(
        c_name='prepare_magnet_plans',
        args=[
            xo.Arg(xo.Int8, pointer=True, name="buffer"),
            xo.Arg(tracker_data_type, name="tracker_data"),
            xo.Arg(MagnetPlanArray, name="magnet_plans"),
            xo.Arg(xo.Int64, pointer=True, name="magnet_plan_elements"),
            xo.Arg(xo.Int64, name="first_plan"),
            xo.Arg(xo.Int64, name="num_plans"),
            xo.Arg(xo.Int8, pointer=True, name="io_buffer"),
            xo.Arg(xo.UInt64, name="track_flags"),
        ],
        n_threads='num_plans',
    )


def attach_prepare_kernel(kernels):
    """Return the track kernel of a build, with the prepare kernel of the
    same build attached as ``prepare_magnet_plans``.

    The two must stay together: the switch in the prepare kernel uses the type
    ids of its build. Track kernels are cached and can be shared between
    trackers (``Tracker.track_kernel``), so the prepare kernel goes with them.
    """
    track_line = kernels['track_line']
    track_line.prepare_magnet_plans = kernels.get('prepare_magnet_plans')
    return track_line


class MagnetPlans:
    """The plans of the magnets of a line, one slot per magnet.

    Made at the first track call of a tracker. It only depends on which
    elements of the line are magnets, which is fixed once the tracker is
    built.

    Attributes:
        plans: the plans (a ``MagnetPlanArray``), on the context.
        element_of_slot: for each slot, the index of its magnet in the line
            (numpy).
        slot_of_element: for each element of the line, its slot, or -1 if it
            is not a magnet (on the context, read by ``track_line``).
    """

    def __init__(self, tracker_data, context):
        is_magnet = [getattr(el, '_magnet_plan', False)
                     for el in tracker_data.elements]

        self.element_of_slot = np.flatnonzero(is_magnet).astype(np.int64)

        slot_of_element = np.full(len(is_magnet), -1, dtype=np.int64)
        slot_of_element[self.element_of_slot] = np.arange(
            len(self.element_of_slot))
        self.slot_of_element = context.nparray_to_context_array(
            slot_of_element)

        num_slots = max(len(self.element_of_slot), 1)  # no empty arrays
        self.plans = MagnetPlanArray(num_slots, _context=context)

        # Read by the prepare kernel; padded, so that it is never empty
        self._element_of_slot_ctx = context.nparray_to_context_array(
            np.concatenate([self.element_of_slot, [0]]))

    def prepare(self, track_kernel, tracker_data, io_buffer, track_flags,
                ele_start, ele_stop):
        """Fill the plans of the magnets in ``[ele_start, ele_stop)`` from the
        current values of their attributes.

        Runs the prepare kernel attached to ``track_kernel``
        (``attach_prepare_kernel``). Does nothing if there is no magnet in the
        range.
        """
        first, stop = np.searchsorted(self.element_of_slot,
                                      [ele_start, ele_stop])
        if stop <= first:
            return
        track_kernel.prepare_magnet_plans(
            buffer=tracker_data._buffer.buffer,
            tracker_data=tracker_data._element_ref_data,
            magnet_plans=self.plans,
            magnet_plan_elements=self._element_of_slot_ctx,
            first_plan=int(first),
            num_plans=int(stop - first),
            io_buffer=io_buffer,
            track_flags=track_flags,
        )

    def read(self, context, element_index):
        """Return the plan of the element ``element_index`` of the line, as it
        is in its slot (a ``MagnetPlan``)."""
        slot = np.searchsorted(self.element_of_slot, element_index)
        if (slot == len(self.element_of_slot)
                or self.element_of_slot[slot] != element_index):
            raise ValueError(f'element {element_index} is not a magnet')
        return self.plans[int(slot)]

    def track_args(self):
        """Return the arguments of ``track_line`` for the plans."""
        return dict(magnet_plans=self.plans,
                    magnet_plan_slots=self.slot_of_element)
