Vector math library | Spencer Tipping
Licensed under the terms of the MIT source code license

# Introduction

This module provides various functions for doing vector math. It is more oriented towards generating new values than it is towards solving equations. It represents vectors as arrays that can be
manipulated with the sequence macro library (though optimized functions will be generated for you, so you should rarely need to do this). All operations in this library are nondestructive.

## Interface

This vector library generalizes to an arbitrary number of coordinates, but the compiled code contains no loops. Instead, you instantiate an N-dimensional copy of the library and it compiles specialized
functions. So, for example, the compiled function for componentwise addition in three-dimensional space is vplus(a, b) = [a[0] + b[0], a[1] + b[1], a[2] + b[2]]. Generally you'd combine this with a
using[] macro to eliminate duplicate compilation:

    reflect(a, normal) = a /-vproj/ normal /-vscale/ -2 /-vplus/ a,
    using [caterwaul.vector(3, 'v')]

If you're using multiple dimensions at once, you can customize the prefix to distinguish the functions:

    using [caterwaul.merge({}, caterwaul.vector(3, 'v3'), caterwaul.vector(4, 'v4'))]

    caterwaul.module('vector', 'js_all', function ($) {
    $.vector = generator(base_v, composite_v),

    where [generator(base, composite)(n, prefix) = {} /compiled_base /-$.merge/ composite(compiled_base) /-rename/ prefix -where [compiled_base = base(n)],
           rename(o, prefix)                     = o %k*['#{prefix || ""}#{x}'] -seq,
           compile_function(xs, e)               = tree /!$.compile -se- it.tree /eq.tree -where [tree = '(function (_formals) {return _e})'.qs /~replace/ {_formals: xs, _e: e}],

# Vector functions

Each function is implemented in terms of its structure. Simple componentwise functions are specified by providing an expression to use for each component, where a wildcard 'i' will be replaced by the
index of that coordinate. This expression is then used in a reduction, which is some structure that combines components into a single value. For example, this is the reduction for 'plus':

    plus = reduction(3,                   // <- number of dimensions
                     'a, b'.qs,           // <- formal parameters, quoted as syntax
                     '[x]'.qs,            // <- result expression
                     'x, y'.qs,           // <- binary combination of intermediate values
                     'a[i] + b[i]'.qs)    // <- componentwise combination

             base_v(n)         = capture [plus  = r(n, 'a, b'.qs, '[x]'.qs, 'x, y'.qs, 'a[i] + b[i]'.qs),           times = r(n, 'a, b'.qs, '[x]'.qs, 'x, y'.qs, 'a[i] * b[i]'.qs),
                                          minus = r(n, 'a, b'.qs, '[x]'.qs, 'x, y'.qs, 'a[i] - b[i]'.qs),           scale = r(n, 'a, b'.qs, '[x]'.qs, 'x, y'.qs, 'a[i] * b'.qs),
                                          dot   = r(n, 'a, b'.qs, 'x'.qs, 'x + y'.qs, 'a[i] * b[i]'.qs),            norm  = r(n, 'a'.qs, 'Math.sqrt(x)'.qs, 'x + y'.qs, 'a[i] * a[i]'.qs),

                                          min   = r(n, 'a, b'.qs, '[x]'.qs, 'x, y'.qs, 'Math.min(a[i], b[i])'.qs),  macv  = r(n, 'a, b, c'.qs, '[x]'.qs, 'x, y'.qs, 'a[i] + b[i] * c[i]'.qs),
                                          max   = r(n, 'a, b'.qs, '[x]'.qs, 'x, y'.qs, 'Math.max(a[i], b[i])'.qs),  macs  = r(n, 'a, b, c'.qs, '[x]'.qs, 'x, y'.qs, 'a[i] + b * c[i]'.qs),
                                                                                                                    mixv  = r(n, 'a, b, c'.qs, '[x]'.qs, 'x, y'.qs, 'a[i] * (1 - b[i]) + c[i] * b[i]'.qs),
                                                                                                                    mixs  = r(n, 'a, b, c'.qs, '[x]'.qs, 'x, y'.qs, 'a[i] * (1 - b) + c[i] * b'.qs)]

                                -where [r(n, formals, wrap, fold, each) = formals /-compile_function/ body
                                                                          -where [body = wrap /~replace/ {x: n[n] *[each /~replace/ {i: '#{x}'}] /[fold /~replace/ {x: x0, y: x}] -seq}]],

             composite_v(base) = capture [unit = ref_compile(base, 'a'.qs,    'scale(a, 1.0 / norm(a))'.qs),
                                          proj = ref_compile(base, 'a, b'.qs, 'scale(b, dot(a, b) / dot(b, b))'.qs),
                                          orth = ref_compile(base, 'a, b'.qs, 'minus(a, scale(b, dot(a, b) / dot(b, b)))'.qs)]

                                -where [ref_compile(functions, formals, body) = formals /-compile_function/ new_body
                                                                                -where [new_body = body |~replace| functions %v*[new $.expression_ref(x.tree)] -seq]]]});