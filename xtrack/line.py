import json
import numpy as np

import xobjects as xo
import xpart as xp

from .loader_sixtrack import _expand_struct
from .loader_mad import iter_from_madx_sequence
from .beam_elements import element_classes, Multipole
from . import beam_elements
from .beam_elements import Drift


import logging

log=logging.getLogger(__name__)

def mk_class_namespace(extra_classes):
    try:
       import xfields as xf
       all_classes= element_classes + xf.element_classes + extra_classes
    except ImportError:
        log.warning("Xfields not installed correctly")

    out=AttrDict()
    for cl in all_classes:
        out[cl.__name__]=cl
    return out


_thick_element_types = (beam_elements.Drift, ) #TODO add DriftExact

def _is_drift(element): # can be removed if length is zero
    return isinstance(element, (beam_elements.Drift,) )

def _is_thick(element):
    return  ((hasattr(element, "isthick") and element.isthick) or
             (isinstance(element, _thick_element_types)))




# missing access to particles._m:
deg2rad = np.pi / 180.

class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self



class Line:

    @classmethod
    def from_dict(cls, dct, _context=None, _buffer=None, classes=()):
        class_dict=mk_class_namespace(classes)

        _buffer, _ =xo.get_a_buffer(size=8,context=_context, buffer=_buffer)
        elements = []
        for el in dct["elements"]:
            eltype = class_dict[el["__class__"]]
            eldct=el.copy()
            del eldct['__class__']
            if hasattr(el,'XoStruct'):
               newel = eltype.from_dict(eldct,_buffer=_buffer)
            else:
               newel = eltype.from_dict(eldct)
            elements.append(newel)

        self = cls(elements=elements, element_names=dct['element_names'])
        if 'particle_ref' in dct.keys():
            self.particle_ref = xp.Particles.from_dict(dct['particle_ref'],
                                    _context=_buffer.context)

        if '_var_manager' in dct.keys():
            import xdeps as xd
            import math
            self._var_management = {}
            vardata = dct['_var_management_data']

            manager = xd.Manager()

            _vref = manager.ref(vardata['var_values'],'vars')
            _eref = manager.ref(vardata['mad_elements_dct'], 'mad_elements_dct')
            _fref = manager.ref(math,'f')
            _lref = manager.ref(self.element_dict, 'line_dict')
            manager.reload(dct['_var_manager'])

            self._var_management['manager'] = manager
            self._var_management['data'] = vardata
            self._var_management['lref'] = _lref
            self._var_management['vref'] = _vref
            self._var_management['fref'] = _fref
            self._var_management['eref'] = _eref
            self.vars = _vref

        return self


    @classmethod
    def from_sixinput(cls, sixinput, classes=()):
        class_dict=mk_class_namespace(classes)

        line_data, rest, iconv = _expand_struct(sixinput, convert=class_dict)

        ele_names = [dd[0] for dd in line_data]
        elements = [dd[2] for dd in line_data]

        line = cls(elements=elements, element_names=ele_names)

        other_info = {}
        other_info["rest"] = rest
        other_info["iconv"] = iconv

        line.other_info = other_info

        return line

    @classmethod
    def from_madx_sequence(
        cls,
        sequence,
        classes=(),
        ignored_madtypes=[],
        exact_drift=False,
        drift_threshold=1e-6,
        deferred_expressions=False,
        install_apertures=False,
        apply_madx_errors=False,
    ):

        class_dict=mk_class_namespace(classes)

        elements = []
        element_names = []
        for el_name, el in iter_from_madx_sequence(
            sequence,
            class_dict,
            ignored_madtypes=ignored_madtypes,
            exact_drift=exact_drift,
            drift_threshold=drift_threshold,
            install_apertures=install_apertures,
        ):
            elements.append(el)
            element_names.append(el_name)

        line = cls(elements=elements, element_names=element_names)

        if deferred_expressions:
            mad = sequence._madx

            # Extract all values
            from collections import defaultdict
            import xdeps as xd
            from xdeps.madxutils import MadxEval
            import math

            # Extract globals values from madx
            _var_values=defaultdict(lambda :0)
            for name,par in mad.globals.cmdpar.items():
                _var_values[name]=par.value

            # Extract element values from madx
            _mad_elements_dct={}
            for name,elem in mad.elements.items():
                elemdata={}
                for parname, par in elem.cmdpar.items():
                    elemdata[parname]=par.value
                _mad_elements_dct[name]=elemdata
                _mad_elements_dct[name]['__basetype__'] = elem.base_type.name

            _ref_manager = manager=xd.Manager()
            _vref=manager.ref(_var_values,'vars')
            _eref=manager.ref(_mad_elements_dct,'mad_elements_dct')
            _fref=manager.ref(math,'f')
            madeval=MadxEval(_vref,_fref,_eref).eval

            # Extract expressions from madx globals
            for name,par in mad.globals.cmdpar.items():
                if par.expr is not None:
                    _vref[name]=madeval(par.expr)

            # Extract expressions from madx elements
            for name,elem in mad.elements.items():
                for parname, par in elem.cmdpar.items():
                    if par.expr is not None:
                        if par.dtype==12: # handle lists
                            for ii,ee in enumerate(par.expr):
                                if ee is not None:
                                    _eref[name][parname][ii]=madeval(ee)
                        else:
                            _eref[name][parname]=madeval(par.expr)

            _lref = manager.ref(line.element_dict, 'line_dict')

            for nn, ee in line.element_dict.items():
                if isinstance(ee, beam_elements.Multipole):
                    assert nn in _mad_elements_dct.keys()
                    ee_mad_dct = _mad_elements_dct[nn]
                    ref_knl = line.element_dict[nn].knl.copy()
                    ref_ksl = line.element_dict[nn].ksl.copy()
                    if ee_mad_dct['__basetype__'] == 'hkicker':
                        _lref[nn].knl[0] = -_eref[nn]['kick']
                    elif ee_mad_dct['__basetype__'] == 'vkicker':
                        _lref[nn].ksl[0] = _eref[nn]['kick']
                    elif ee_mad_dct['__basetype__'] == 'multipole':
                        _lref[nn].knl = _eref[nn]['knl']
                        _lref[nn].ksl[0] = _eref[nn]['ksl'][0]
                    elif ee_mad_dct['__basetype__'] in ['tkicker', 'kicker']:
                        if hasattr(ee_mad_dct, 'hkick'):
                            _lref[nn].knl[0] = -_eref[nn]['hkick']
                        if hasattr(ee_mad_dct, 'vkick'):
                            _lref[nn].ksl[0] = _eref[nn]['vkick']
                    else:
                        raise ValueError('???')
                    assert np.allclose(line.element_dict[nn].knl, ref_knl, 1e-18)
                    assert np.allclose(line.element_dict[nn].ksl, ref_ksl, 1e-18)

            line.vars = _vref
            line._var_management = {}
            line._var_management['data'] = {}
            line._var_management['data']['var_values'] = _var_values
            line._var_management['data']['mad_elements_dct'] = _mad_elements_dct

            line._var_management['manager'] = _ref_manager
            line._var_management['lref'] = _lref
            line._var_management['vref'] = _vref
            line._var_management['fref'] = _fref
            line._var_management['eref'] = _eref


        if apply_madx_errors:
            line._apply_madx_errors(sequence)

        return line

    def __init__(self, elements=(), element_names=None, particle_ref=None):
        if isinstance(elements,dict):
            element_dict=elements
            if element_names is None:
                raise ValueError('`element_names must be provided'
                                 ' if `elements` is a dictionary.')
        else:
            if element_names is None:
                element_names = [ f"e{ii}" for ii in range(len(elements))]
            element_dict = dict(zip(element_names, elements))

        self.element_dict=element_dict
        self.element_names=element_names

        self.particle_ref = particle_ref

        self._var_management = None
        self.vars = None

    @property
    def elements(self):
        return tuple([self.element_dict[nn] for nn in self.element_names])

    def filter_elements(self, mask=None, exclude_types_starting_with=None):

        if mask is None:
            assert exclude_types_starting_with is not None

        if exclude_types_starting_with is not None:
            assert mask is None
            mask = [not(ee.__class__.__name__.startswith(exclude_types_starting_with))
                    for ee in self.elements]

        new_elements = []
        assert len(mask) == len(self.elements)
        for ff, ee in zip(mask, self.elements):
            if ff:
                new_elements.append(ee)
            else:
                if _is_thick(ee) and not _is_drift(ee):
                    new_elements.append(Drift(length==ee.length))
                else:
                    new_elements.append(Drift(length=0))

        return self.__class__(elements=new_elements,
                              element_names=self.element_names)

    def cycle(self, index_first_element=None, name_first_element=None):

        if ((index_first_element is not None and name_first_element is not None)
               or (index_first_element is None and name_first_element is None)):
             raise ValueError(
                "Plaese provide either `index_first_element` or `name_first_element`.")

        if name_first_element is not None:
            assert self.element_names.count(name_first_element) == 1, (
                f"name_first_element={name_first_element} occurs more than once!"
            )
            index_first_element = self.element_names.index(name_first_element)

        new_elements = (list(self.elements[index_first_element:])
                        + list(self.elements[:index_first_element]))
        new_element_names = (list(self.element_names[index_first_element:])
                        + list(self.element_names[:index_first_element]))

        return self.__class__(
                         elements=new_elements, element_names=new_element_names)

    def _freeze(self):
        self.element_names = tuple(self.element_names)

    def _frozen_check(self):
        if isinstance(self.element_name, tuple):
            raise ValueError(
                'This action is not allowed as the line is frozen!')

    def __len__(self):
        return len(self.element_names)

    def to_dict(self):
        out = {}
        out["elements"] = [el.to_dict() for el in self.elements]
        out["element_names"] = self.element_names[:]
        if self.particle_ref is not None:
            out['particle_ref'] = self.particle_ref.to_dict()
        if self._var_management is not None:
            out['_var_management_data'] = self._var_management['data']
            out['_var_manager'] = self._var_management['manager'].dump()
        return out

    def insert_element(self, idx, element, name):
        self._frozen_check()
        assert name not in self.element_dict.keys()
        self.element_dict[name] = element
        self.element_names.insert(idx, name)
        return self

    def append_element(self, element, name):
        self._frozen_check()
        assert name not in self.element_dict.keys()
        self.element_dict[name] = element
        self.element_names.append(name)
        return self

    def get_length(self):
        ll = 0
        for ee in self.elements:
            if _is_thick(ee):
                ll += ee.length

        return ll

    def get_s_elements(self, mode="upstream"):

        assert mode in ["upstream", "downstream"]
        s_prev = 0
        s = []
        for ee in self.elements:
            if mode == "upstream":
                s.append(s_prev)
            if _is_thick(ee):
                s_prev += ee.length
            if mode == "downstream":
                s.append(s_prev)
        return s

    def remove_inactive_multipoles(self, inplace=False):

        self._frozen_check()

        newline = Line(elements=[], element_names=[])

        for ee, nn in zip(self.elements, self.element_names):
            if isinstance(ee, (beam_elements.Multipole)):
                aux = [ee.hxl, ee.hyl] + list(ee.knl) + list(ee.ksl)
                if np.sum(np.abs(np.array(aux))) == 0.0:
                    continue
            newline.append_element(ee, nn)

        if inplace:
            self.element_names = newline.element_names
            return self
        else:
            return newline

    def remove_zero_length_drifts(self, inplace=False):

        self._frozen_check()

        newline = Line(elements=[], element_names=[])

        for ee, nn in zip(self.elements, self.element_names):
            if _is_drift(ee):
                if ee.length == 0.0:
                    continue
            newline.append_element(ee, nn)

        if inplace:
            self.element_names = newline.element_names
            return self
        else:
            return newline

    def merge_consecutive_drifts(self, inplace=False):

        self._frozen_check()

        newline = Line(elements=[], element_names=[])

        for ee, nn in zip(self.elements, self.element_names):
            if len(newline.elements) == 0:
                newline.append_element(ee, nn)
                continue

            if _is_drift(ee):
                prev_ee = newline.elements[-1]
                prev_nn = newline.element_names[-1]
                if _is_drift(prev_ee):
                    prev_ee.length += ee.length
                    prev_nn += ('_' + nn)
                    newline.element_names[-1] = prev_nn
                else:
                    newline.append_element(ee, nn)
            else:
                newline.append_element(ee, nn)

        if inplace:
            self.element_dict.update(newline.element_dict)
            self.element_names = newline.element_names
            return self
        else:
            return newline

    def merge_consecutive_multipoles(self, inplace=False):

        self._frozen_check()
        if hasattr(self, '_var_management'):
            raise NotImplementedError('`merge_consecutive_multipoles` not'
                                      ' available when deferred expressions are'
                                      ' used')

        newline = Line(elements=[], element_names=[])

        for ee, nn in zip(self.elements, self.element_names):
            if len(newline.elements) == 0:
                newline.append_element(ee, nn)
                continue

            if isinstance(ee, beam_elements.Multipole):
                prev_ee = newline.elements[-1]
                prev_nn = newline.element_names[-1]
                if (isinstance(prev_ee, beam_elements.Multipole)
                    and prev_ee.hxl==ee.hxl==0 and prev_ee.hyl==ee.hyl==0
                    ):

                    oo=max(len(prev_ee.knl), len(prev_ee.ksl),
                           len(ee.knl), len(ee.ksl))
                    knl=np.zeros(oo,dtype=float)
                    ksl=np.zeros(oo,dtype=float)
                    for ii,kk in enumerate(prev_ee.knl):
                        knl[ii]+=kk
                    for ii,kk in enumerate(ee.knl):
                        knl[ii]+=kk
                    for ii,kk in enumerate(prev_ee.ksl):
                        ksl[ii]+=kk
                    for ii,kk in enumerate(ee.ksl):
                        ksl[ii]+=kk
                    newee = beam_elements.Multipole(
                            knl=knl, ksl=ksl, hxl=prev_ee.hxl, hyl=prev_ee.hyl,
                            length=prev_ee.length,
                            radiation_flag=prev_ee.radiation_flag)
                    prev_nn += ('_' + nn)
                    newline.element_names[-1] = prev_nn
                    newline.elements[-1] = newee
                else:
                    newline.append_element(ee, nn)
            else:
                newline.append_element(ee, nn)

        if inplace:
            self.element_dict.update(newline.element_dict)
            self.element_names = newline.element_names
            return self
        else:
            return newline

    def get_elements_of_type(self, types):
        if not hasattr(types, "__iter__"):
            type_list = [types]
        else:
            type_list = types

        names = []
        elements = []
        for ee, nn in zip(self.elements, self.element_names):
            for tt in type_list:
                if isinstance(ee, tt):
                    names.append(nn)
                    elements.append(ee)

        return elements, names

    def _find_element_ids(self, element_name):
        """Find element_name in this Line instance's
        self.elements_name list. Assumes the names are unique.

        Return index before and after the element, taking into account
        attached _aperture instances (LimitRect, LimitEllipse, ...)
        which would follow the element occurrence in the list.

        Raises IndexError if element_name not found in this Line.
        """
        # will raise error if element not present:
        idx_el = self.element_names.index(element_name)
        try:
            # if aperture marker is present
            idx_after_el = self.element_names.index(element_name + "_aperture") + 1
        except ValueError:
            # if aperture marker is not present
            idx_after_el = idx_el + 1
        return idx_el, idx_after_el

    def _add_offset_error_to(self, element_name, dx=0, dy=0):
        idx_el, idx_after_el = self._find_element_ids(element_name)
        xyshift = beam_elements.XYShift(dx=dx, dy=dy)
        inv_xyshift = beam_elements.XYShift(dx=-dx, dy=-dy)
        self.insert_element(idx_el, xyshift, element_name + "_offset_in")
        self.insert_element(
            idx_after_el + 1, inv_xyshift, element_name + "_offset_out"
        )

    def _add_aperture_offset_error_to(self, element_name, arex=0, arey=0):
        idx_el, idx_after_el = self._find_element_ids(element_name)
        idx_el_aper = idx_after_el - 1
        if not self.element_names[idx_el_aper] == element_name + "_aperture":
            # it is allowed to provide arex/arey without providing an aperture
            print('Info: Element', element_name, ': arex/y provided without aperture -> arex/y ignored')
            return
        xyshift = beam_elements.XYShift(dx=arex, dy=arey)
        inv_xyshift = beam_elements.XYShift(dx=-arex, dy=-arey)
        self.insert_element(idx_el_aper, xyshift, element_name + "_aperture_offset_in")
        self.insert_element(
            idx_after_el + 1, inv_xyshift, element_name + "_aperture_offset_out"
        )

    def _add_tilt_error_to(self, element_name, angle):
        '''Alignment error of transverse rotation around s-axis.
        The element corresponding to the given `element_name`
        gets wrapped by SRotation elements with rotation angle
        `angle`.

        In the case of a thin dipole component, the corresponding
        curvature terms in the Multipole (hxl and hyl) are rotated
        by `angle` as well.
        '''
        idx_el, idx_after_el = self._find_element_ids(element_name)
        element = self.elements[self.element_names.index(element_name)]
        if isinstance(element, beam_elements.Multipole) and (
                element.hxl or element.hyl):
            dpsi = angle * deg2rad

            hxl0 = element.hxl
            hyl0 = element.hyl

            hxl1 = hxl0 * np.cos(dpsi) - hyl0 * np.sin(dpsi)
            hyl1 = hxl0 * np.sin(dpsi) + hyl0 * np.cos(dpsi)

            element.hxl = hxl1
            element.hyl = hyl1
        srot = beam_elements.SRotation(angle=angle)
        inv_srot = beam_elements.SRotation(angle=-angle)
        self.insert_element(idx_el, srot, element_name + "_tilt_in")
        self.insert_element(idx_after_el + 1, inv_srot, element_name + "_tilt_out")

    def _add_multipole_error_to(self, element_name, knl=[], ksl=[]):
        # will raise error if element not present:
        assert element_name in self.element_names
        element_index = self.element_names.index(element_name)
        element = self.elements[element_index]

        new_order = max([len(knl), len(ksl), len(element.knl), len(element.ksl)])
        new_knl = new_order*[0]
        new_ksl = new_order*[0]

        # Original strengths
        for ii, vv in enumerate(element.knl):
            new_knl[ii] += element.knl[ii]
        for ii, vv in enumerate(element.ksl):
            new_ksl[ii] += element.ksl[ii]

        # Errors
        for ii, vv in enumerate(knl):
            new_knl[ii] += knl[ii]
        for ii, vv in enumerate(ksl):
            new_ksl[ii] += ksl[ii]

        new_element = Multipole(knl=new_knl, ksl=new_ksl,
                length=element.length, hxl=element.hxl,
                hyl=element.hyl, radiation_flag=element.radiation_flag)

        self.elements[element_index] = new_element

    def _apply_madx_errors(self, madx_sequence):
        """Applies errors from MAD-X sequence to existing
        elements in this Line instance.

        Return names of MAD-X elements with existing align_errors
        or field_errors which were not found in the elements of
        this Line instance (and thus not treated).

        Example via cpymad:
            madx = cpymad.madx.Madx()

            # (...set up lattice and errors in cpymad...)

            seq = madx.sequence.some_lattice
            line = Line.from_madx_sequence(
                                    seq,
                                    apply_madx_errors=True
                              )
        """
        elements_not_found = []
        for element, element_name in zip(
                madx_sequence.expanded_elements,
                madx_sequence.expanded_element_names()
        ):
            if element_name not in self.element_names:
                if element.align_errors or element.field_errors:
                    elements_not_found.append(element_name)
                    continue

            if element.align_errors:
                # add offset
                dx = element.align_errors.dx
                dy = element.align_errors.dy
                if dx or dy:
                    self._add_offset_error_to(element_name, dx, dy)

                # add tilt
                dpsi = element.align_errors.dpsi
                if dpsi:
                    self._add_tilt_error_to(element_name, angle=dpsi / deg2rad)

                # add aperture-only offset
                arex = element.align_errors.arex
                arey = element.align_errors.arey
                if arex or arey:
                    self._add_aperture_offset_error_to(element_name, arex, arey)

                # check for errors which cannot be treated yet:
                #for error_type in dir(element.align_errors):
                 #   if not error_type[0] == '_' and \
                  #          error_type not in ['dx', 'dy', 'dpsi', 'arex',
                   #                            'arey', 'count', 'index']:
                        #print(
                        #    f'Warning: MAD-X error type "{error_type}"'
                        #    " not implemented yet."
                        #)

            if element.field_errors:
                # add multipole error
                if any(element.field_errors.dkn) or \
                            any(element.field_errors.dks):
                    knl = element.field_errors.dkn
                    ksl = element.field_errors.dks
                    on=np.where(knl)[0]
                    os=np.where(ksl)[0]
                    on = on[-1] if len(on)>0 else 0
                    os = os[-1] if len(os)>0 else 0
                    oo = max(os,on)+1
                    knl = knl[:oo]  # delete trailing zeros
                    ksl = ksl[:oo]  # to keep order low
                    self._add_multipole_error_to(element_name, knl, ksl)

        return elements_not_found


