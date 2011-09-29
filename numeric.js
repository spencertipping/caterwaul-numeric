// Linear algebra library | Spencer Tipping
// Licensed under the terms of the MIT source code license

// Introduction.
// This module provides various functions for doing vector and matrix algebra. It is more oriented towards generating new values than it is towards solving equations. It represents all matrices
// and vectors as arrays that can be manipulated with the sequence macro library. This representation is not particularly fast when used directly, but the accompanying numerical compiler can
// reduce memory allocation to improve performance significantly. All operations in this library are nondestructive.

// Matrices are represented as nested arrays in row-major order. This means that x[1][2] returns the entry in the second row, third column. All matrices have the invariant that they are properly
// rectangular; that is, all row arrays are the same length.

  // Interface.
//   This vector library generalizes to an arbitrary number of coordinates, but the compiled code contains no loops. Instead, you instantiate an N-dimensional copy of the library and it compiles
//   specialized functions. So, for example, the compiled function for componentwise addition in three-dimensional space is vplus(a, b) = [a[0] + b[0], a[1] + b[1], a[2] + b[2]]. Generally you'd
//   combine this with a using[] macro to eliminate duplicate compilation:

  // | reflect(a, normal) = a /-vproj/ normal /-vscale/ -2 /-vplus/ a,
//     using [caterwaul.linear.vector(3, 'v')]

  // If you're using multiple dimensions at once, you can customize the prefix to distinguish the functions:

  // | using [caterwaul.merge({}, caterwaul.linear.vector(3, 'v3'), caterwaul.linear.vector(4, 'v4'))]

caterwaul.js_all()(function ($) {

// Vector functions.
// Each function is implemented in terms of its structure. Simple componentwise functions are specified by providing an expression to use for each component, where a wildcard 'i' will be replaced
// by the index of that coordinate. This expression is then used in a reduction, which is some structure that combines components into a single value. For example, this is the reduction for
// 'plus':

// | plus = reduction(3,                   // <- number of dimensions
//                    'a, b'.qs,           // <- formal parameters, quoted as syntax
//                    '[x]'.qs,            // <- result expression
//                    'x, y'.qs,           // <- binary combination of intermediate values
//                    'a[i] + b[i]'.qs)    // <- componentwise combination

// Most simple vector functions can be defined this way. Others, however, are better defined in terms of each other; for instance, the 'proj' and 'orth' functions never access the vectors
// directly since they are defined in terms of the dot product and vector-scalar multiplication. The obvious solution is to first create the reduction functions and then define things like 'proj'
// and 'orth' to close over them; however, this is problematic from a compilation perspective since we would need the closure state to know the dimension of 'proj' and 'orth'. (Not problematic if
// we compile it here, but problematic if we push it off to the numerical compiler.)

// Better is to inline the functions using compile-time beta expansion. This follows a fairly standard form; inlining a function f within an expression involves allocating a variable for each
// formal, then replacing the function call with the returned expression (after substituting formal names). So, for instance, if we had this setup:

// | plus(a, b) = [a[0] + b[0], a[1] + b[1], a[2] + b[2]],
//   twice(a)   = a /-plus/ a

// We could inline 'plus' into 'twice' like this:

// | twice = function (a) {
//     var p1 = a;
//     var p2 = a;
//     return [p1[0] + p2[0], p1[1] + p2[1], p1[2] + p2[2]];
//   }

// This code contains no closures or other functional elements, so the numerical compiler will have an easy time inlining it further to eliminate the array allocation. Note, however, that this
// strategy doesn't eliminate unnecessary intermediate array allocations. That's also a job for the numerical compiler.

// Caveat: I don't have inlining implemented yet. Right now I'm just compiling any function names into refs and relying on the runtime to inline things.

// Note that cross products aren't handled by the vector library; these are considered to be matrix functions, and the determinant formula is computed specifically for a given matrix size rather
// than being explicitly generalized.

  // Defining alternative componentwise semantics.
//   An extra parameter, 'field', lets you redefine algebraic field operations. This can be useful if you want to build vectors or matrices over non-scalar data structures. In particular, if you
//   specify this parameter you'll need to provide replacements for +, -, *, and /. Here's an example for working with complex numbers:

  // | my_field = {'+': '{r: x.r + y.r, i: x.i + y.i}'.qs,
//                 '-': '{r: x.r - y.r, i: x.i + y.i}'.qs,
//                 '*': '{r: x.r * x.r - x.i * y.i, i: 2 * x.i * y.i}'.qs,
//                 '/': '{r: (x.r*y.r + x.i*y.i) / (y.r*y.r + y.i*y.i), i: (x.i*y.r - x.r*y.i) / (y.r*y.r + y.i*y.i)}'.qs}

  // These expressions will then replace the real-number field used by default. Note here that the complex conjugate operation has duplicated subexpressions; y.r*y.r + y.i*y.i is computed twice.
//   This won't be a problem in the compiled function because all of the expressions are subject to common subexpression elimination prior to being compiled. (This is the mechanism used to
//   optimize matrix array access as well.) Because of this optimization, it's very important that any side-effects of each subexpression be idempotent and reorderable.

  $.linear = capture [vector(n, prefix, field) = {} /base_functions /-$.merge/ high_level_v(base_functions) /-rename/ prefix -where [base_functions = base_v(n, field || scalar_ops)],
                      matrix(n, prefix, field) = {} /base_functions /-$.merge/ high_level_m(base_functions) /-rename/ prefix -where [base_functions = base_m(n, field || scalar_ops)]],

  where [rename(o, prefix)  = o %k*['#{prefix || ""}#{x}'] -seq,

         base_v(n)          = capture [plus  = r(n, 'a, b'.qs, '[x]'.qs, 'x, y'.qs, 'a[i] + b[i]'.qs),  times = r(n, 'a, b'.qs, '[x]'.qs, 'x, y'.qs, 'a[i] * b[i]'.qs),
                                       minus = r(n, 'a, b'.qs, '[x]'.qs, 'x, y'.qs, 'a[i] - b[i]'.qs),  scale = r(n, 'a, b'.qs, '[x]'.qs, 'x, y'.qs, 'a[i] * b'.qs),
                                       dot   = r(n, 'a, b'.qs, 'x'.qs, 'x + y'.qs, 'a[i] * b[i]'.qs),   norm  = r(n, 'a'.qs, 'Math.sqrt(x)'.qs, 'x + y'.qs, 'a[i] * a[i]'.qs)]

                      -where [r(n, formals, wrap, fold, each) = f /!$.compile -se [it.body = body, it.formals = formals]
                                                        -where [body = wrap /~replace/ {x: n[n] *[each.replace({i: '#{x}'})] /[fold.replace({x: x0, y: x})] -seq},
                                                                f    = '(function (_formals) {return _e})'.qs /~replace/ {_formals: formals, _e: body}]],

         high_level_v(base) = capture [unit = ref_compile(base, 'a'.qs,    'scale(a, 1 / norm(a))'.qs),
                                       proj = ref_compile(base, 'a, b'.qs, 'scale(b, dot(a, b) / dot(b, b))'.qs),
                                       orth = ref_compile(base, 'a, b'.qs, 'minus(a, scale(b, dot(a, b) / dot(b, b)))'.qs)]

                      -where [ref_compile(functions, formals, body) = $.compile('(function (_formals) {return _e})'.qs /~replace/
                                                                                {_formals: formals,
                                                                                 _e:       body |~replace| functions %v*[new $.ref(x)] -seq})],

// Matrix functions.
// Most of these are standard textbook functions, though there some of them are peculiar to this data representation. In particular, all matrix coordinates are unrolled; this means that some
// weird optimizations can happen. Vector functions were essentially flat; there was very little repetitive access to sub-arrays. This isn't true of matrices, however. Consider a simple
// coordinate-wise addition function over 2x2 matrices:

// | plus(a, b) = [[a[0][0] + b[0][0], a[0][1] + b[0][1]],
//                 [a[1][0] + b[1][0], a[1][1] + b[1][1]]]

// Each top-level sub-array in a and b is accessed twice (that is, a[0], a[1], b[0], and b[1]), and unless the Javascript runtime is clever enough to prove their invariance, these loads will
// happen twice. Rather than explicitly loading the sub-arrays twice, better is to perform common subexpression elimination and allocate local variables to cache the lookups:

// | plus = function (a, b) {
//     var a0 = a[0], a1 = a[1], b0 = b[0], b1 = b[1];
//     return [[a0[0] + b0[0], a0[1] + b0[1]],
//             [a1[0] + b1[0], a1[1] + b1[1]]];
//   };

// Matrix functions are structured roughly the same way as vector functions from a compilation perspective. It's a bit more complicated here because there are two levels of reduction instead of
// one. Some things are also complexified by conditions on the matrix size; for instance, the determinant only exists for square matrices. (Fortunately, this library only provides functions for
// square matrices.)

// There are some compromises made for performance. In particular, matrix and vector functions are untyped, so it doesn't make sense to compute a cross product in the usual vector way. (That is,
// set the top row of the matrix to contain vectors instead of scalars.) Instead, the cross product is a special form of the determinant; this preserves the untyped representation. In addition to
// things like this, various safety rules are ignored; for instance, there is no size-checking despite the fact that every function operates only on matrices of specific dimensions.

         base_m_functions_for(n)          = capture [plus  = componentwise(n, 'a, b'.qs, 'a[i][j] + b[i][j]'.qs),  scale     = componentwise(n, 'a, b'.qs, 'a[i][j] * b'.qs),
                                                     minus = componentwise(n, 'a, b'.qs, 'a[i][j] - b[i][j]'.qs),  transpose = componentwise(n, 'a'.qs, 'a[j][i]'.qs),

                                                     times = r3(n, 'a, b'.qs, '[x]'.qs, 'x, y'.qs, '[x]'.qs, 'x, y'.qs, 'x'.qs, 'x + y'.qs, 'a[i][k] * b[k][j]'.qs),

                                                     determinant = r3(n, 'a'.qs, 

                                                                                               ]})(caterwaul);

// Generated by SDoc 
