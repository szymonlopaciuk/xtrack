import builtins
import math
import operator
import re
import traceback
from collections import defaultdict
from pathlib import Path
from typing import Tuple, Optional

import cython as cy
import numpy as np
import scipy.constants as sc
from cython.cimports.cpython.ref import Py_INCREF, PyObject  # noqa
from cython.cimports.libc.stdint import uintptr_t  # noqa
from cython.cimports.libc.stdio import fopen  # noqa
from cython.cimports.posix.unistd import access, R_OK  # noqa
from cython.cimports.xtrack.sequence import parser as xld  # noqa

import xdeps as xd
import xobjects as xo
import xtrack as xt
from xdeps.refs import CallRef, LiteralExpr
from xtrack.beam_elements import element_classes as xt_element_classes

KEEP_LITERAL_EXPRESSIONS = False


BUILTIN_CONSTANTS = {
    # Supported constants from MAD-X manual
    'pi': math.pi,
    'twopi': 2 * math.pi,
    'degrad': 180 / math.pi,  # °/rad
    'raddeg': math.pi / 180,  # rad/°
    'e': math.e,
    'emass': sc.electron_mass * sc.c**2 / sc.e,  # eV
    'pmass': sc.proton_mass * sc.c**2 / sc.e,  # eV
    'nmass': sc.neutron_mass * sc.c**2 / sc.e,  # eV
    'umass': sc.m_u * sc.c**2 / sc.e,  # eV
    'mumass': sc.value('muon mass') * sc.c**2 / sc.e,  # eV
    'clight': sc.c,  # m/s
    'qelect': sc.value('elementary charge'),  # A * s
    'hbar': sc.hbar,  # eV * s
    'erad': sc.value('classical electron radius'),  # m
    'prad': sc.value('classical electron radius') * (sc.m_e / sc.m_p),  # m
}

AVAILABLE_ELEMENT_CLASSES = {cls.__name__: cls for cls in xt_element_classes}

try:
    from xfields import element_classes as xf_element_classes
    AVAILABLE_ELEMENT_CLASSES.update({cls.__name__: cls for cls in xf_element_classes})
except ModuleNotFoundError:
    pass

try:
    from xcoll import element_classes as xc_element_classes
    AVAILABLE_ELEMENT_CLASSES.update({cls.__name__: cls for cls in xc_element_classes})
except ModuleNotFoundError:
    pass


class ParseLogEntry:
    def __init__(self, location, reason, error=True, context=None):
        self.line_no = location['first_line']
        self.column = location['first_column']
        self.end_line_no = location['last_line']
        self.end_column = location['last_column']
        self.reason = reason
        self.context = context
        self.error = error

    def add_file_line(self, file_lines):
        relevant_lines = file_lines[self.line_no - 1:self.end_line_no]

        def _insert(string, where, what):
            return string[:where] + what + string[where:]

        relevant_lines[-1] = _insert(relevant_lines[-1], self.end_column - 1, '\033[0m')
        relevant_lines[0] = _insert(relevant_lines[0], self.column - 1, '\033[4m')

        self.context = '\n'.join(relevant_lines)

    def __repr__(self):
        if self.error:
            out = '\033[91mError \033[0m'
        else:
            out = '\033[93mWarning \033[0m'
        out += f'on line {self.line_no} column {self.column}: {self.reason}'

        if self.context:
            out += '\n\n'
            out += _indent_string(self.context, indent='    > ')

        return out


class ParseError(Exception):
    LIMIT = 15

    def __init__(self, parse_log, error_type="error", file_lines=None):
        parse_log = parse_log.copy()

        if len(parse_log) > 1:
            message = f'{error_type.title()}s occurred while parsing:\n\n'
        else:
            message = ''

        truncated = parse_log[:self.LIMIT]
        for entry in truncated:
            entry.add_file_line(file_lines)
        message += '\n\n'.join([repr(entry) for entry in truncated])

        if len(parse_log) > self.LIMIT:
            message += (
                f'\n\nTruncated {len(parse_log) - self.LIMIT} more '
                f'{error_type}s...'
            )

        super().__init__(message)


@cy.cfunc
def register_error(scanner: xld.yyscan_t, exception, action, add_context=True, location=None):
    parser = parser_from_scanner(scanner)
    caught_exc_string = '\n'.join(traceback.format_exception(exception))
    caught_exc_string = _indent_string(caught_exc_string, indent='    > ')

    full_error_string = (
        f"While {action} the following error occurred:\n\n"
        f"{caught_exc_string}"
    )
    parser.handle_error(full_error_string, add_context=add_context, location=location)


def _indent_string(string, indent='    '):
    return re.sub('^', indent, string, flags=re.MULTILINE)


@cy.cclass
class Parser:
    scanner = cy.declare(xld.yyscan_t)
    log = cy.declare(list, visibility='public')
    xd_manager = cy.declare(object, visibility='public')
    vars = cy.declare(object, visibility='public')
    var_refs = cy.declare(object, visibility='public')
    elements = cy.declare(dict, visibility='public')
    element_refs = cy.declare(object, visibility='public')
    functions = cy.declare(object, visibility='public')
    func_refs = cy.declare(object, visibility='public')
    lines = cy.declare(dict, visibility='public')
    global_elements = cy.declare(dict, visibility='public')
    _context = cy.declare(object, visibility='public')

    def __init__(self, _context=xo.context_default):
        self.log = []

        self.xd_manager = xd.Manager()

        self.vars = defaultdict(lambda: 0)
        self.vars.default_factory = None
        self.var_refs = self.xd_manager.ref(self.vars, 'vars')

        self.elements = {}
        self.element_refs = self.xd_manager.ref(self.elements, 'element_refs')

        self.functions = xt.line.Functions()
        self.func_refs = self.xd_manager.ref(self.functions, 'f')

        self.lines = {}
        self.global_elements = {}

        self._context = _context
        # So we do something interesting here: the objective is that we want to
        # keep a pointer to the parser object on the C side. However, for
        # reasons yet unknown, in the process of passing the pointer to C and
        # then taking it back, python decreases the reference count leading to
        # the premature death of the object (TODO: figure out why it happens?).
        #
        # It seems that we have two options: manually increasing the reference
        # count here, to inform Python that the object is in C, or to somehow
        # make this reference weak. The first solution is not ideal, as we
        # create a cycle: both the python-side parser and the C-side scanner
        # effectively hold strong references to one another.
        #
        # The implementation of the second solution is not super pretty, but it
        # is legitimate -- we cast to PyObject* to skip reference counting:
        xld.yylex_init_extra(
            cy.cast(cy.pointer(PyObject), self),
            cy.address(self.scanner),
        )

    def __del__(self):
        xld.yylex_destroy(self.scanner)

    def parse_string(self, string):
        xld.yy_scan_string(string.encode(), self.scanner)

        # yy_scan_string doesn't reset the line and column, so we do it manually
        xld.yyset_lineno(1, self.scanner)
        xld.yyset_column(1, self.scanner)

        success = xld.yyparse(self.scanner)

        self._assert_success(success, string)

    def parse_file(self, path):
        c_path = path.encode()

        file_readable = access(c_path, R_OK)
        if file_readable != 0:
            raise OSError(f'File `{path}` does not exist or is not readable.')

        xld.yyset_in(fopen(c_path, "r"), self.scanner)

        success = xld.yyparse(self.scanner)

        self._assert_success(success, Path(path))

    def _assert_success(self, success, input_file_or_string=None):
        errors = [entry for entry in self.log if entry.error]

        def _make_error(error_type='error'):
            input = input_file_or_string
            if isinstance(input_file_or_string, Path):
                input = input_file_or_string.read_text()
            parse_error = ParseError(
                self.log, error_type=error_type, file_lines=input.split('\n'))
            return parse_error

        if success != 0 or errors:
            raise _make_error()

        if self.log:
            parse_warning = _make_error(error_type='warning')
            xt.general._print(str(parse_warning))

        self.log = []  # reset log after all errors have been reported

    def get_line(self, name):
        if not self.lines:
            raise ValueError('No sequence was parsed. Either the input has no '
                             'sequences, or the parser has not been run yet.')

        if name is None:
            try:
                name, = self.lines.keys()
            except ValueError:
                raise ValueError('Cannot unambiguously determine the sequence '
                                 'as no name was provided and there is more '
                                 'than one sequence in parsed input.')

        return self.lines[name]

    def handle_error(self, message, add_context=True, location=None):
        self._handle_parse_log_event(
            message, is_error=True, add_context=add_context, location=location)

    def handle_warning(self, message, add_context=True, location=None):
        self._handle_parse_log_event(
            message, is_error=False, add_context=add_context, location=location)

    def _handle_parse_log_event(self, message, is_error, add_context, location):
        text = xld.yyget_text(self.scanner).decode()
        yylloc_ptr = xld.yyget_lloc(self.scanner)
        location = location or yylloc_ptr[0]
        log_entry = ParseLogEntry(location, message, error=is_error)
        if add_context:
            log_entry.context = text
        self.log.append(log_entry)

    def get_identifier_ref(self, identifier, location):
        try:
            if identifier in self.vars:
                return self.var_refs[identifier]

            if identifier in BUILTIN_CONSTANTS:
                return py_float(self.scanner, BUILTIN_CONSTANTS[identifier])

            self.handle_error(
                f'use of an undefined variable `{identifier}`',
                location=location,
            )
            return np.nan
        except Exception as e:
            register_error(self.scanner, e, 'getting an identifier reference')

    def set_value(self, identifier, value):
        self.var_refs[identifier] = value

    def add_global_element(self, name, parent, args):
        self.global_elements[name] = (parent, args)

    def add_line(self, line_name, elements, params):
        try:
            if elements is None or params is None:
                # A parsing error occurred, and we're in recovery.
                # Let's give up on making the line, it won't be any good anyway.
                return

            self.elements[line_name] = local_elements = {}
            element_names = []

            for el_template in elements:
                if el_template is None:
                    # A parsing error occurred and has been registered already.
                    # Ignore this element.
                    continue

                name, parent, args = el_template

                self.add_element_to_line(line_name, name, parent, args)

                parent_name = getattr(local_elements[name], 'parent_name', None)
                if parent_name and parent_name not in local_elements:
                    self.add_element_to_line(
                        line_name,
                        parent_name,
                        *self.global_elements[parent_name],
                    )

                element_names.append(name)

            line = xt.Line(
                elements=self.elements[line_name],
                element_names=element_names,
            )
            line._var_management = {}
            line._var_management['data'] = {}
            line._var_management['data']['var_values'] = self.vars
            line._var_management['data']['functions'] = self.functions

            line._var_management['manager'] = self.xd_manager
            line._var_management['vref'] = self.var_refs
            line._var_management['lref'] = self.element_refs[line_name]
            line._var_management['fref'] = self.func_refs
            self.lines[line_name] = line
        except Exception as e:
            register_error(self.scanner, e, 'building the sequence')

    def add_element_to_line(self, line_name, name, parent, args):
        local_elements = self.elements[line_name]

        if not args and parent not in AVAILABLE_ELEMENT_CLASSES:  # simply insert a replica
            local_elements[name] = xt.Replica(parent_name=parent)
            return

        element_cls = AVAILABLE_ELEMENT_CLASSES[parent]
        kwargs = {k: self._ref_get_value(v) for k, v in args}
        element = element_cls.from_dict(kwargs, _context=self._context)
        local_elements[name] = element
        element_ref = self.element_refs[line_name][name]

        for k, v in args:
            if isinstance(v, list):
                list_ref = getattr(element_ref, k)
                for i, le in enumerate(v):
                    if xd.refs.is_ref(le):
                        list_ref[i] = le
            elif xd.refs.is_ref(v):
                setattr(element_ref, k, v)

    @staticmethod
    def _ref_get_value(value_or_ref):
        if isinstance(value_or_ref, list):
            return [getattr(elem, '_value', elem) for elem in value_or_ref]
        return getattr(value_or_ref, '_value', value_or_ref)


def parser_from_scanner(yyscanner) -> Parser:
    # Cast back to Python, see comment in Parser.__init__:
    parser = xld.yyget_extra(yyscanner)
    return cy.cast(Parser, parser)


def yyerror(_, yyscanner, message):
    parser = parser_from_scanner(yyscanner)
    parser.handle_error(message.decode())


def py_float(scanner, value):
    return py_numeric(scanner, value)


def py_integer(scanner, value):
    return py_numeric(scanner, value)


@cy.exceptval(check=False)
def py_numeric(scanner, value):
    try:
        if KEEP_LITERAL_EXPRESSIONS:
            return LiteralExpr(value)
        return value
    except Exception as e:
        register_error(scanner, e, f'parsing a numeric value')


@cy.exceptval(check=False)
def py_unary_op(scanner, op_string, value):
    try:
        function = getattr(operator, op_string.decode())
        return function(value)
    except Exception as e:
        register_error(scanner, e, f'parsing a unary operation')


@cy.exceptval(check=False)
def py_binary_op(scanner, op_string, left, right):
    try:
        function = getattr(operator, op_string.decode())
        return function(left, right)
    except Exception as e:
        register_error(scanner, e, f'parsing a binary operation')


@cy.exceptval(check=False)
def py_call_func(scanner, func_name, value):
    try:
        name = func_name.decode()
        if name == 'const':
            return value._value if xd.refs.is_ref(value) else value

        parser = parser_from_scanner(scanner)
        if name not in parser.functions:
            parser_from_scanner(scanner).handle_error(
                f'builtin function `{name}` is unknown',
            )
            return np.nan

        return getattr(parser.func_refs, name)(value)
    except Exception as e:
        register_error(scanner, e, f'parsing a function call')


def py_assign(scanner, identifier, value):
    try:
        return identifier.decode(), value
    except Exception as e:
        register_error(scanner, e, f'parsing an assignment')


def py_arrow(scanner, source_name, field_name):
    try:
        raise NotImplementedError
    except Exception as e:
        register_error(scanner, e, f'parsing the arrow syntax')


def py_identifier_atom(scanner, name, location):
    try:
        parser = parser_from_scanner(scanner)
        return parser.get_identifier_ref(name.decode(), location)
    except Exception as e:
        register_error(
            scanner, e, f'parsing an identifier',
            add_context=True, location=location,
        )


def py_set_value(scanner, assignment, location):
    try:
        identifier, value = assignment
        parser = parser_from_scanner(scanner)
        if identifier in BUILTIN_CONSTANTS:
            parser.handle_warning(
                f"variable `{identifier}` shadows a built-in constant",
                add_context=False,
                location=location,
            )
        parser.set_value(identifier, value)
    except Exception as e:
        register_error(
            scanner, e, 'parsing a deferred assign statement',
            add_context=True, location=location,
        )


def py_make_sequence(scanner, name, args, elements):
    try:
        parser = parser_from_scanner(scanner)
        parser.add_line(name.decode(), elements, dict(args))
    except Exception as e:
        register_error(scanner, e, 'parsing a sequence')


def py_clone(scanner, name, parent, args) -> Optional[Tuple[str, str, dict]]:
    try:
        if name.decode() in AVAILABLE_ELEMENT_CLASSES:
            parser = parser_from_scanner(scanner)
            parser.handle_error(f'the name `{name.decode()}` shadows a built-in type.')
            return None

        return name.decode(), parent.decode(), args
    except Exception as e:
        register_error(scanner, e, 'parsing a clone statement')


def py_clone_global(scanner, clone):
    try:
        if clone is None:  # A parsing error already occurred, recover
            return
        parser = parser_from_scanner(scanner)
        name, parent, args = clone
        parser.add_global_element(name, parent, args)
    except Exception as e:
        register_error(scanner, e, 'parsing a global clone statement')
