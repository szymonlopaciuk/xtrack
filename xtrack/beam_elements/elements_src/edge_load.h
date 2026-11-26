// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2026.                 //
// ######################################### //
#ifndef XTRACK_EDGE_LOAD_H
#define XTRACK_EDGE_LOAD_H

#include "xtrack/beam_elements/elements_src/edge_flags.h"

// Fill the settings of one face (a MagnetFace of the plan) from the flat edge fields
// of an element. GET is the getter prefix (e.g. BendData_get_, or
// ThinSliceBendData_get__parent_ for a slice), side is entry or exit.
// Wanted, not implemented: an exit face following the entry one, which would
// resolve the exit mask and the NaN exit parameters here (see edge_flags.h).

#define EDGE_FLAGS(GET, el, side)        GET##_edge_##side##_flags(el)
#define EDGE_PARAM(GET, el, side, param) GET##edge_##side##_##param(el)

// Straight magnets: mask, soft dipole fringe (fringe integral, half gap) and
// the integrals of the soft quadrupole fringe (Lee-Whiting), f1 and f2
// (one statement each: they can follow an if)
#define EDGE_LOAD_STRAIGHT(face, GET, el, side) do {            \
    MagnetFace_set_flags(face, EDGE_FLAGS(GET, el, side));      \
    MagnetFace_set_fint(face, EDGE_PARAM(GET, el, side, fint)); \
    MagnetFace_set_hgap(face, EDGE_PARAM(GET, el, side, hgap)); \
    MagnetFace_set_f1(face, EDGE_PARAM(GET, el, side, qint1));  \
    MagnetFace_set_f2(face, EDGE_PARAM(GET, el, side, qint2));  \
} while (0)

// Bends and magnets: as straight magnets, plus face angle and linear feed-down
#define EDGE_LOAD_BEND(face, GET, el, side) do {                              \
    EDGE_LOAD_STRAIGHT(face, GET, el, side);                                  \
    MagnetFace_set_angle(face, EDGE_PARAM(GET, el, side, angle));             \
    MagnetFace_set_angle_fdown(face, EDGE_PARAM(GET, el, side, angle_fdown)); \
} while (0)

#endif  // XTRACK_EDGE_LOAD_H
