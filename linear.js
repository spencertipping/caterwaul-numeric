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

// Implementation.
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

// Note that cross products aren't handled by the vector library; these are considered to be matrix functions, and the determinant formula is computed specifically for a given matrix size rather
// than being explicitly generalized.

  $.linear = capture [vector(n,    prefix) = {} /base_functions /-$.merge/ high_level_v_functions_for(base_functions) /-rename/ prefix -where [base_functions = base_v_functions_for(n)],
                      matrix(n, m, prefix) = raise [new Error('not yet implemented')]],

  where [rename(o, prefix)                = o %k*['#{prefix}#{x}'] -seq,

         base_v_functions_for(n)          = capture [plus  = r(n, 'a, b'.qs, '[x]'.qs, 'x, y'.qs, 'a[i] + b[i]'.qs),  times = r(n, 'a, b'.qs, '[x]'.qs, 'x, y'.qs, 'a[i] * b[i]'.qs),
                                                     minus = r(n, 'a, b'.qs, '[x]'.qs, 'x, y'.qs, 'a[i] - b[i]'.qs),  scale = r(n, 'a, b'.qs, '[x]'.qs, 'x, y'.qs, 'a[i] * b'.qs),

                                                     dot   = r(n, 'a, b'.qs, 'x'.qs, 'x + y'.qs, 'a[i] * b[i]'.qs),   norm  = r(n, 'a'.qs, 'Math.sqrt(x)'.qs, 'x + y'.qs, 'a[i] * a[i]'.qs)]

                                    -where [r(n, formals, wrap, fold, each) = '(function (_formals) {return _e})'.qs.replace(
                                                                              {_formals: formals,
                                                                               _e:       wrap.replace({x: n[n] *[each.replace({i: '#{x}'})] /[fold.replace({x: x0, y: x})] -seq})})],

         high_level_v_functions_for(base) = capture [unit = inline(base, 'a'.qs,    'scale(a, 1 / norm(a))'.qs),
                                                     proj = inline(base, 'a, b'.qs, 'scale(b, dot(a, b) / dot(b, b))'.qs),
                                                     orth = inline(base, 'a, b'.qs, 'minus(a, scale(b, dot(a, b) / dot(b, b)))'.qs)]

                                    -where [inline(base, formals, body) = '(function (_formals) {_vars; return _e})'.qs.replace(
                                                                          {_formals: formals, _vars: intermediate_formals, _e: beta_expanded_body})

                                                                  -where [intermediate_formals = [],    // <- constructed statefully
                                                                          destructured_bases   = base %v*destructure -seq,
                                                                          beta_expanded_body   = ]

                                                                  -where [destructure(f) = 'function (_formals) {return _body}'.qs.match(f /!$.parse)]]]})(caterwaul);

// Generated by SDoc 
