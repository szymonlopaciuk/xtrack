# copyright ############################### #
# This file is part of the Xtrack Package.  #
# Copyright (c) CERN, 2021.                 #
# ######################################### #

from pathlib import Path
import numpy as np
from functools import partial

import xobjects as xo
import xpart as xp
import xtrack

from xobjects.hybrid_class import _build_xofields_dict

from .general import _pkg_root
from .internal_record import RecordIdentifier, RecordIndex, generate_get_record

start_per_part_block = """
   int64_t const n_part = LocalParticle_get__num_active_particles(part0); //only_for_context cpu_serial cpu_openmp
   #pragma omp parallel for                                       //only_for_context cpu_openmp
   for (int jj=0; jj<n_part; jj+=!!CHUNK_SIZE!!){                 //only_for_context cpu_serial cpu_openmp
    //#pragma omp simd
    for (int iii=0; iii<!!CHUNK_SIZE!!; iii++){                   //only_for_context cpu_serial cpu_openmp
      int const ii = iii+jj;                                      //only_for_context cpu_serial cpu_openmp
      if (ii<n_part){                                             //only_for_context cpu_serial cpu_openmp

        LocalParticle lpart = *part0;//only_for_context cpu_serial cpu_openmp
        LocalParticle* part = &lpart;//only_for_context cpu_serial cpu_openmp
        part->ipart = ii;            //only_for_context cpu_serial cpu_openmp

        LocalParticle* part = part0;//only_for_context opencl cuda
""".replace("!!CHUNK_SIZE!!", "128")

end_part_part_block = """
     } //only_for_context cpu_serial cpu_openmp
    }  //only_for_context cpu_serial cpu_openmp
   }   //only_for_context cpu_serial cpu_openmp
"""

def _handle_per_particle_blocks(sources, local_particle_src):

    if isinstance(sources, str):
        sources = (sources, )
        wasstring = True
    else:
        wasstring = False

    out = []
    for ii, ss in enumerate(sources):
        if isinstance(ss, Path):
            with open(ss, 'r') as fid:
                strss = fid.read()
        else:
            strss = ss

        strss = strss.replace('/*placeholder_for_local_particle_src*/',
                                local_particle_src,
                                )
        if '//start_per_particle_block' in strss:

            lines = strss.splitlines()
            for ill, ll in enumerate(lines):
                if '//start_per_particle_block' in ll:
                    lines[ill] = start_per_part_block
                if '//end_per_particle_block' in ll:
                    lines[ill] = end_part_part_block

            # TODO: this is very dirty, just for check!!!!!
            out.append('\n'.join(lines))
        else:
            out.append(ss)


    if wasstring:
        out = out[0]

    return out

def _generate_per_particle_kernel_from_local_particle_function(
                                                element_name, kernel_name,
                                                local_particle_function_name,
                                                additional_args=[]):

    if len(additional_args) > 0:
        add_to_signature = ", ".join([
            f"{' /*gpuglmem*/ ' if arg.pointer else ''} {arg.get_c_type()} {arg.name}"
                for arg in additional_args]) + ", "
        add_to_call = ", " + ", ".join(f"{arg.name}" for arg in additional_args)

    source = ('''
            /*gpukern*/
            '''
            f'void {kernel_name}(\n'
            f'               {element_name}Data el,\n'
'''
                             ParticlesData particles,
'''
            f'{(add_to_signature if len(additional_args) > 0 else "")}'
'''
                             int64_t flag_increment_at_element,
                /*gpuglmem*/ int8_t* io_buffer){
            LocalParticle lpart;
            lpart.io_buffer = io_buffer;

            int64_t part_id = 0;                    //only_for_context cpu_serial cpu_openmp
            int64_t part_id = blockDim.x * blockIdx.x + threadIdx.x; //only_for_context cuda
            int64_t part_id = get_global_id(0);                    //only_for_context opencl

            int64_t part_capacity = ParticlesData_get__capacity(particles);
            if (part_id<part_capacity){
                Particles_to_LocalParticle(particles, &lpart, part_id);
                if (check_is_active(&lpart)>0){
'''
            f'      {local_particle_function_name}(el, &lpart{(add_to_call if len(additional_args) > 0 else "")});\n'
'''
                }
                if (check_is_active(&lpart)>0 && flag_increment_at_element){
                        increment_at_element(&lpart);
                }
            }
        }
''')
    return source

class MetaBeamElement(xo.MetaHybridClass):

    def __new__(cls, name, bases, data):
        _XoStruct_name = name+'Data'
        particles_class = xp.ParticlesInterface

        # Take xofields from data['_xofields'] or from bases
        xofields = _build_xofields_dict(bases, data)
        data = data.copy()
        data['_xofields'] = xofields

        depends_on = []
        extra_c_source = [
            _pkg_root.joinpath('headers','constants.h'),
            _pkg_root.joinpath('headers','checks.h'),
            _pkg_root.joinpath('headers','particle_states.h')
        ]
        kernels = {}

        # Handle internal record
        if '_internal_record_class' in data.keys():
            data['_xofields']['_internal_record_id'] = RecordIdentifier
            if '_skip_in_to_dict' not in data.keys():
                data['_skip_in_to_dict'] = []
            data['_skip_in_to_dict'].append('_internal_record_id')

            depends_on.append(RecordIndex)
            depends_on.append(data['_internal_record_class']._XoStruct)
            extra_c_source.append(
                generate_get_record(ele_classname=_XoStruct_name,
                    record_classname=data['_internal_record_class']._XoStruct.__name__))

        # Get user-defined source, dependencies and kernels
        if '_extra_c_sources' in data.keys():
            extra_c_source.extend(data['_extra_c_sources'])

        if '_depends_on' in data.keys():
            depends_on.extend(data['_depends_on'])

        if '_kernels' in data.keys():
            kernels.update(data['_kernels'])

        # Generate track kernel
        extra_c_source.append(
            _generate_per_particle_kernel_from_local_particle_function(
                element_name=name, kernel_name=name+'_track_particles',
                local_particle_function_name=name+'_track_local_particle'))

        # Add dependency on Particles class
        depends_on.append(particles_class._XoStruct)

        # Define track kernel
        track_kernel_name = f'{name}_track_particles'
        kernels[track_kernel_name] = xo.Kernel(
                    args=[xo.Arg(xo.ThisClass, name='el'),
                        xo.Arg(particles_class._XoStruct, name='particles'),
                        xo.Arg(xo.Int64, name='flag_increment_at_element'),
                        xo.Arg(xo.Int8, pointer=True, name="io_buffer")]
                    )

        # Generate per-particle kernels
        if '_per_particle_kernels' in data.keys():
            for nn, kk in data['_per_particle_kernels'].items():
                extra_c_source.append(
                    _generate_per_particle_kernel_from_local_particle_function(
                        element_name=name, kernel_name=nn,
                        local_particle_function_name=kk.c_name,
                        additional_args=kk.args))
                if particles_class._XoStruct not in depends_on:
                    depends_on.append(particles_class._XoStruct)

                kernels.update(
                    {nn:
                        xo.Kernel(args=[xo.Arg(xo.ThisClass, name='el'),
                            xo.Arg(particles_class._XoStruct, name='particles')]
                            + kk.args + [
                            xo.Arg(xo.Int64, name='flag_increment_at_element'),
                            xo.Arg(xo.Int8, pointer=True, name="io_buffer")])}
                )


        # Call HybridClass metaclass
        data['_depends_on'] = depends_on
        data['_extra_c_sources'] = extra_c_source
        data['_kernels'] = kernels
        new_class = xo.MetaHybridClass.__new__(cls, name, bases, data)

        # Attach some information to the class
        new_class._track_kernel_name = track_kernel_name
        if '_internal_record_class' in data.keys():
            new_class._XoStruct._internal_record_class = data['_internal_record_class']
            new_class._internal_record_class = data['_internal_record_class']

        # Attach methods corresponding to per-particle kernels
        if '_per_particle_kernels' in data.keys():
            for nn in data['_per_particle_kernels'].keys():
                setattr(new_class, nn, PerParticlePyMethodDescriptor(kernel_name=nn))

        return new_class


class BeamElement(xo.HybridClass, metaclass=MetaBeamElement):

    iscollective = None
    isthick = False
    behaves_like_drift = False
    allow_backtrack = False
    skip_in_loss_location_refinement = False
    _kernels_and_classes = None

    def __init__(self, *args, **kwargs):
        xo.HybridClass.__init__(self, *args, **kwargs)

    def init_pipeline(self, pipeline_manager, name, partners_names=[]):
        self._pipeline_manager = pipeline_manager
        self.name = name
        self.partners_names = partners_names

    def compile_kernels(self, particles_class, *args, **kwargs):
        if 'apply_to_source' not in kwargs.keys():
            kwargs['apply_to_source'] = []
        kwargs['apply_to_source'].append(
            partial(_handle_per_particle_blocks,
                    local_particle_src=particles_class.gen_local_particle_api()))
        xo.HybridClass.compile_kernels(self,
                                       extra_classes=[particles_class._XoStruct],
                                       *args, **kwargs)

    def _store_kernels_and_classes(self, kernels, classes):
        if not self.__class__._kernels_and_classes:
            self.__class__._kernels_and_classes = []

        self.__class__._kernels_and_classes.append(
            (self.context, kernels, classes)
        )

    def _get_kernels_and_classes(self):
        if not self.__class__._kernels_and_classes:
            return None, None

        for context, kernels, classes in self.__class__._kernels_and_classes:
            if context == self.context:
                return kernels, classes
            if type(context) == type(self.context):
                if isinstance(context, xo.ContextCpu):
                    return kernels, classes

        return None, None

    def _track_with_minitracker(self, particles, increment_at_element=False):
        if hasattr(self, 'io_buffer') and self.io_buffer is not None:
            io_buffer = self.io_buffer
        else:
            from xtrack import new_io_buffer
            io_buffer = new_io_buffer(capacity=1, _context=self.context)

        # Restore the kernels if the context is right
        track_kernels, element_classes = self._get_kernels_and_classes()

        from xtrack.line import Line
        line = Line(elements=[self])
        tracker = line.build_tracker(
            particles_class=particles.__class__,
            io_buffer=io_buffer,
            track_kernel=track_kernels,
            element_classes=element_classes,
            compile=False,
        )
        tracker.config.DANGER_SKIP_ACTIVE_CHECK_AND_SWAPS = (not increment_at_element)
        tracker.config.XTRACK_MULTIPOLE_NO_SYNRAD = False
        tracker.skip_end_turn_actions = True

        self._store_kernels_and_classes(tracker.track_kernel, tracker.element_classes)

        tracker.io_buffer = io_buffer
        tracker.track(particles)

    def _track_per_particle(self, particles, increment_at_element=False):
        context = self._buffer.context
        if not hasattr(self, '_track_kernel'):
            if self._track_kernel_name not in context.kernels.keys():
                self.compile_kernels(particles_class=particles.__class__)
            self._track_kernel = context.kernels[self._track_kernel_name]

        if hasattr(self, 'io_buffer') and self.io_buffer is not None:
            io_buffer_arr = self.io_buffer.buffer
        else:
            io_buffer_arr = context.zeros(1, dtype=np.int8)  # dummy

        self._track_kernel.description.n_threads = particles._capacity
        self._track_kernel(el=self._xobject, particles=particles,
                           flag_increment_at_element=increment_at_element,
                           io_buffer=io_buffer_arr)

    def track(self, particles, increment_at_element=False):
        old_start_at_element = particles.start_tracking_at_element
        particles.start_tracking_at_element = -1
        old_is_collective = self.iscollective
        self.iscollective = False
        try:
            self._track_with_minitracker(particles, increment_at_element)
        finally:
            particles.start_tracking_at_element = old_start_at_element
            self.iscollective = old_is_collective


    @property
    def context(self):
        return self._buffer.context

    def _arr2ctx(self, arr):
        ctx = self._buffer.context

        if isinstance(arr, list):
            arr = np.array(arr)

        if np.isscalar(arr):
            if hasattr(arr, 'item'):
                return arr.item()
            else:
                return arr
        elif isinstance(arr, ctx.nplike_array_type):
            return arr
        elif isinstance(arr, np.ndarray):
            return ctx.nparray_to_context_array(arr)
        else:
            raise ValueError("Invalid array type")


class PerParticlePyMethod:

    def __init__(self, kernel_name, element):
        self.kernel_name = kernel_name
        self.element = element

    def __call__(self, particles, increment_at_element=False, **kwargs):
        instance = self.element
        context = instance.context
        if not hasattr(instance, '_track_kernel'):
            if instance._track_kernel_name not in context.kernels.keys():
                instance.compile_kernels(particles_class=particles.__class__)
            instance._track_kernel = context.kernels[instance._track_kernel_name]
        self.kernel = context.kernels[self.kernel_name]

        if hasattr(self.element, 'io_buffer') and self.element.io_buffer is not None:
            io_buffer_arr = self.element.io_buffer.buffer
        else:
            context = self.kernel.context
            io_buffer_arr = context.zeros(1, dtype=np.int8)  # dummy

        self.kernel.description.n_threads = particles._capacity
        self.kernel(el=self.element._xobject,
                    particles=particles,
                    flag_increment_at_element=increment_at_element,
                    io_buffer=io_buffer_arr,
                    **kwargs)


class PerParticlePyMethodDescriptor:
    def __init__(self, kernel_name):
        self.kernel_name = kernel_name

    def __get__(self, instance, owner):
        return PerParticlePyMethod(kernel_name=self.kernel_name,
                                   element=instance)
