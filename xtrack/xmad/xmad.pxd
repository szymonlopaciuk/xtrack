#cython: language_level=3
# This file defines the interface between the C part of the parser and the
# Python/Cython modules. We expose the C constructs we need from flex/bison code
# to Python, and declare the C functions ('yyerror') that are expected by the
# parser, which we then implement in the corresponding Cython 'xmad.py' file.

from libc.stdio cimport FILE
from cpython.ref cimport PyObject

ctypedef PyObject* extra_t

ctypedef fused numeric:
    double
    long

cdef extern from "xmad_tab.h":
    ctypedef struct YYLTYPE:
        int first_line
        int first_column
        int last_line
        int last_column

    int yyparse(yyscan_t)

cdef extern from "xmad_lex.h":
    ctypedef struct YY_BUFFER_STATE:
        pass

    ctypedef struct yyscan_t:
        pass

    void* yylex_init(yyscan_t*)
    void* yylex_init_extra(extra_t, yyscan_t*)
    int yylex_destroy(yyscan_t)

    void yyset_lineno(int, yyscan_t)
    void yyset_column(int, yyscan_t)
    YY_BUFFER_STATE* yy_scan_string(const char*, yyscan_t)

    char* yyget_text(yyscan_t)
    extra_t yyget_extra(yyscan_t)
    YYLTYPE* yyget_lloc(yyscan_t)

    char* yyset_in(FILE*, yyscan_t)
    void yyset_extra(extra_t, yyscan_t)

cdef object parser_from_scanner(yyscan_t yyscanner)

cdef public void yyerror(YYLTYPE* yyllocp, yyscan_t yyscanner, const char* message)
cdef public object py_integer(yyscan_t scanner, long value)
cdef public object py_float(yyscan_t scanner, double value)
cdef public object py_numeric(yyscan_t scanner, numeric value)
cdef public object py_unary_op(yyscan_t scanner, const char* op_string, object value)
cdef public object py_binary_op(yyscan_t scanner, const char* op_string, object left, object right)
cdef public tuple py_eq_value_scalar(yyscan_t scanner, const char* identifier, object value)
cdef public tuple py_eq_defer_scalar(yyscan_t scanner, const char* identifier, object value)
cdef public object py_call_func(yyscan_t scanner, const char* func_name, object value)
cdef public object py_arrow(yyscan_t scanner, const char* source_name, const char* field_name)
cdef public object py_identifier_atom(yyscan_t scanner, const char* name)
cdef public void py_set_defer(yyscan_t scanner, tuple assignment)
cdef public void py_set_value(yyscan_t scanner, tuple assignment)
cdef public void py_make_sequence(yyscan_t scanner, const char* name, list args, list elements)
cdef public object py_clone(yyscan_t scanner, const char* name, const char* parent, list args)
cdef public object py_eq_value_sum(yyscan_t scanner, const char* name, object value)
cdef public object py_eq_defer_sum(yyscan_t scanner, const char* name, object value)
cdef public object py_eq_value_array(yyscan_t scanner, const char* name, list array)
cdef public object py_eq_defer_array(yyscan_t scanner, const char* name, list array)
