 caterwaul.module( 'vector' , (function (qs,qs1,qs2,qs3,qs4,qs5,qs6,qs7,qs8,qs9,qsa,qsb,qsc,qsd,qse,qsf,qsg,qsh,qsi,qsj,qsk,qsl,qsm,qsn,qso,qsp,qsq,qsr,qss,qst,qsu,qsv,qsw,qsx,qsy,qsz,qs10,qs11,qs12,qs13,qs14,qs15,qs16,qs17,qs18,qs19,qs1a,qs1b,qs1c,qs1d,qs1e,qs1f,qs1g,qs1h,qs1i) {var result= ( function ($) { (function () {var generator =function (base, composite) { ; return function (n, prefix) { ; return(function () {var compiled_base = base(n) ; return rename($.merge( {} ,compiled_base, composite(compiled_base)) , prefix)}) .call(this)}} , rename =function (o, prefix) { ; return(function (xs) {var x, x0, xi, xl, xr;var xr = new xs.constructor() ; for (var x in xs) if (Object.prototype.hasOwnProperty.call(xs, x)) xr[ ( '' + (prefix || "") + '' + (x) + '')] = xs[x] ; return xr}) .call(this, o)} , compile_function =function (xs, e) { ; return(function () {var tree = (qs) .replace( {_formals: xs, _e: e}) ; return(function (it) {return it.tree =tree, it}) .call(this, ($.compile( tree)))}) .call(this)} , base_v =function (n) { ; return(function () {var r =function (n, formals, wrap, fold, each) { ; return(function () {var body = ( wrap) .replace( {x: (function (xs) {var x, x0, xi, xl, xr;for (var x0 = xs[0] , xi = 1, xl = xs.length; xi < xl; ++xi) x = xs[xi] , x0 = ( (fold) .replace( {x: x0, y: x})) ; return x0}) .call(this, (function (xs) {var x, x0, xi, xl, xr;for (var xr = new xs.constructor() , xi = 0, xl = xs.length; xi < xl; ++xi) x = xs[xi] , xr.push( ( (each) .replace( {i: ( '' + (x) + '')}))) ; return xr}) .call(this, (function (i, u, s) {if ( (u - i) * s <= 0) return [] ; for (var r = [] , d = u - i; d >= 0 ? i < u: i > u; i += s) r.push(i) ; return r}) ( (0) , (n) , (1))))}) ; return compile_function( formals, body)}) .call(this)} ; return{plus: r(n,qs1,qs2,qs3,qs4) , times: r(n,qs5,qs6,qs7,qs8) , minus: r(n,qs9,qsa,qsb,qsc) , scale: r(n,qsd,qse,qsf,qsg) , dot: r(n,qsh,qsi,qsj,qsk) , norm: r(n,qsl,qsm,qsn,qso) , min: r(n,qsp,qsq,qsr,qss) , macv: r(n,qst,qsu,qsv,qsw) , max: r(n,qsx,qsy,qsz,qs10) , macs: r(n,qs11,qs12,qs13,qs14) , mixv: r(n,qs15,qs16,qs17,qs18) , mixs: r(n,qs19,qs1a,qs1b,qs1c)}}) .call(this)} , composite_v =function (base) { ; return(function () {var ref_compile =function (functions, formals, body) { ; return(function () {var new_body = ( body) .replace( (function (xs) {var x, x0, xi, xl, xr;var xr = new xs.constructor() ; for (var k in xs) if (Object.prototype.hasOwnProperty.call(xs, k)) x = xs[k] , xr[k] = (new $.expression_ref(x.tree)) ; return xr}) .call(this, functions)) ; return compile_function( formals, new_body)}) .call(this)} ; return{unit: ref_compile(base,qs1d,qs1e) , proj: ref_compile(base,qs1f,qs1g) , orth: ref_compile(base,qs1h,qs1i)}}) .call(this)} ; return $.vector = generator(base_v, composite_v)}) .call(this)}) ;result.caterwaul_expression_ref_table = {qs: ( " new caterwaul.syntax( \"(\" ,new caterwaul.syntax( \"function\" ,new caterwaul.syntax( \"(\" ,new caterwaul.syntax( \"_formals\")) .prefix( \" \") ,new caterwaul.syntax( \"{\" ,new caterwaul.syntax( \"return\" ,new caterwaul.syntax( \"_e\") .prefix( \" \"))) .prefix( \" \")))") ,qs1: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"b\") .prefix( \" \"))") ,qs2: ( " new caterwaul.syntax( \"[\" ,new caterwaul.syntax( \"x\"))") ,qs3: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"x\") ,new caterwaul.syntax( \"y\") .prefix( \" \"))") ,qs4: ( " new caterwaul.syntax( \"+\" ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"i\")) ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"b\") .prefix( \" \") ,new caterwaul.syntax( \"i\"))) .prefix( \" \")") ,qs5: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"b\") .prefix( \" \"))") ,qs6: ( " new caterwaul.syntax( \"[\" ,new caterwaul.syntax( \"x\"))") ,qs7: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"x\") ,new caterwaul.syntax( \"y\") .prefix( \" \"))") ,qs8: ( " new caterwaul.syntax( \"*\" ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"i\")) ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"b\") .prefix( \" \") ,new caterwaul.syntax( \"i\"))) .prefix( \" \")") ,qs9: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"b\") .prefix( \" \"))") ,qsa: ( " new caterwaul.syntax( \"[\" ,new caterwaul.syntax( \"x\"))") ,qsb: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"x\") ,new caterwaul.syntax( \"y\") .prefix( \" \"))") ,qsc: ( " new caterwaul.syntax( \"-\" ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"i\")) ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"b\") .prefix( \" \") ,new caterwaul.syntax( \"i\"))) .prefix( \" \")") ,qsd: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"b\") .prefix( \" \"))") ,qse: ( " new caterwaul.syntax( \"[\" ,new caterwaul.syntax( \"x\"))") ,qsf: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"x\") ,new caterwaul.syntax( \"y\") .prefix( \" \"))") ,qsg: ( " new caterwaul.syntax( \"*\" ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"i\")) ,new caterwaul.syntax( \"b\") .prefix( \" \")) .prefix( \" \")") ,qsh: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"b\") .prefix( \" \"))") ,qsi: ( " new caterwaul.syntax( \"x\")") ,qsj: ( " new caterwaul.syntax( \"+\" ,new caterwaul.syntax( \"x\") ,new caterwaul.syntax( \"y\") .prefix( \" \")) .prefix( \" \")") ,qsk: ( " new caterwaul.syntax( \"*\" ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"i\")) ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"b\") .prefix( \" \") ,new caterwaul.syntax( \"i\"))) .prefix( \" \")") ,qsl: ( " new caterwaul.syntax( \"a\")") ,qsm: ( " new caterwaul.syntax( \"()\" ,new caterwaul.syntax( \".\" ,new caterwaul.syntax( \"Math\") ,new caterwaul.syntax( \"sqrt\")) ,new caterwaul.syntax( \"x\"))") ,qsn: ( " new caterwaul.syntax( \"+\" ,new caterwaul.syntax( \"x\") ,new caterwaul.syntax( \"y\") .prefix( \" \")) .prefix( \" \")") ,qso: ( " new caterwaul.syntax( \"*\" ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"i\")) ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"a\") .prefix( \" \") ,new caterwaul.syntax( \"i\"))) .prefix( \" \")") ,qsp: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"b\") .prefix( \" \"))") ,qsq: ( " new caterwaul.syntax( \"[\" ,new caterwaul.syntax( \"x\"))") ,qsr: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"x\") ,new caterwaul.syntax( \"y\") .prefix( \" \"))") ,qss: ( " new caterwaul.syntax( \"()\" ,new caterwaul.syntax( \".\" ,new caterwaul.syntax( \"Math\") ,new caterwaul.syntax( \"min\")) ,new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"i\")) ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"b\") .prefix( \" \") ,new caterwaul.syntax( \"i\"))))") ,qst: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"b\") .prefix( \" \")) ,new caterwaul.syntax( \"c\") .prefix( \" \"))") ,qsu: ( " new caterwaul.syntax( \"[\" ,new caterwaul.syntax( \"x\"))") ,qsv: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"x\") ,new caterwaul.syntax( \"y\") .prefix( \" \"))") ,qsw: ( " new caterwaul.syntax( \"+\" ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"i\")) ,new caterwaul.syntax( \"*\" ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"b\") .prefix( \" \") ,new caterwaul.syntax( \"i\")) ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"c\") .prefix( \" \") ,new caterwaul.syntax( \"i\"))) .prefix( \" \")) .prefix( \" \")") ,qsx: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"b\") .prefix( \" \"))") ,qsy: ( " new caterwaul.syntax( \"[\" ,new caterwaul.syntax( \"x\"))") ,qsz: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"x\") ,new caterwaul.syntax( \"y\") .prefix( \" \"))") ,qs10: ( " new caterwaul.syntax( \"()\" ,new caterwaul.syntax( \".\" ,new caterwaul.syntax( \"Math\") ,new caterwaul.syntax( \"max\")) ,new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"i\")) ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"b\") .prefix( \" \") ,new caterwaul.syntax( \"i\"))))") ,qs11: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"b\") .prefix( \" \")) ,new caterwaul.syntax( \"c\") .prefix( \" \"))") ,qs12: ( " new caterwaul.syntax( \"[\" ,new caterwaul.syntax( \"x\"))") ,qs13: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"x\") ,new caterwaul.syntax( \"y\") .prefix( \" \"))") ,qs14: ( " new caterwaul.syntax( \"+\" ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"i\")) ,new caterwaul.syntax( \"*\" ,new caterwaul.syntax( \"b\") .prefix( \" \") ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"c\") .prefix( \" \") ,new caterwaul.syntax( \"i\"))) .prefix( \" \")) .prefix( \" \")") ,qs15: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"b\") .prefix( \" \")) ,new caterwaul.syntax( \"c\") .prefix( \" \"))") ,qs16: ( " new caterwaul.syntax( \"[\" ,new caterwaul.syntax( \"x\"))") ,qs17: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"x\") ,new caterwaul.syntax( \"y\") .prefix( \" \"))") ,qs18: ( " new caterwaul.syntax( \"+\" ,new caterwaul.syntax( \"*\" ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"i\")) ,new caterwaul.syntax( \"(\" ,new caterwaul.syntax( \"-\" ,new caterwaul.syntax( \"1\") ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"b\") .prefix( \" \") ,new caterwaul.syntax( \"i\"))) .prefix( \" \")) .prefix( \" \")) .prefix( \" \") ,new caterwaul.syntax( \"*\" ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"c\") .prefix( \" \") ,new caterwaul.syntax( \"i\")) ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"b\") .prefix( \" \") ,new caterwaul.syntax( \"i\"))) .prefix( \" \")) .prefix( \" \")") ,qs19: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"b\") .prefix( \" \")) ,new caterwaul.syntax( \"c\") .prefix( \" \"))") ,qs1a: ( " new caterwaul.syntax( \"[\" ,new caterwaul.syntax( \"x\"))") ,qs1b: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"x\") ,new caterwaul.syntax( \"y\") .prefix( \" \"))") ,qs1c: ( " new caterwaul.syntax( \"+\" ,new caterwaul.syntax( \"*\" ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"i\")) ,new caterwaul.syntax( \"(\" ,new caterwaul.syntax( \"-\" ,new caterwaul.syntax( \"1\") ,new caterwaul.syntax( \"b\") .prefix( \" \")) .prefix( \" \")) .prefix( \" \")) .prefix( \" \") ,new caterwaul.syntax( \"*\" ,new caterwaul.syntax( \"[]\" ,new caterwaul.syntax( \"c\") .prefix( \" \") ,new caterwaul.syntax( \"i\")) ,new caterwaul.syntax( \"b\") .prefix( \" \")) .prefix( \" \")) .prefix( \" \")") ,qs1d: ( " new caterwaul.syntax( \"a\")") ,qs1e: ( " new caterwaul.syntax( \"()\" ,new caterwaul.syntax( \"scale\") ,new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"/\" ,new caterwaul.syntax( \"1.0\") .prefix( \" \") ,new caterwaul.syntax( \"()\" ,new caterwaul.syntax( \"norm\") .prefix( \" \") ,new caterwaul.syntax( \"a\"))) .prefix( \" \")))") ,qs1f: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"b\") .prefix( \" \"))") ,qs1g: ( " new caterwaul.syntax( \"()\" ,new caterwaul.syntax( \"scale\") ,new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"b\") ,new caterwaul.syntax( \"/\" ,new caterwaul.syntax( \"()\" ,new caterwaul.syntax( \"dot\") .prefix( \" \") ,new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"b\") .prefix( \" \"))) ,new caterwaul.syntax( \"()\" ,new caterwaul.syntax( \"dot\") .prefix( \" \") ,new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"b\") ,new caterwaul.syntax( \"b\") .prefix( \" \")))) .prefix( \" \")))") ,qs1h: ( " new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"b\") .prefix( \" \"))") ,qs1i: ( " new caterwaul.syntax( \"()\" ,new caterwaul.syntax( \"minus\") ,new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"()\" ,new caterwaul.syntax( \"scale\") .prefix( \" \") ,new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"b\") ,new caterwaul.syntax( \"/\" ,new caterwaul.syntax( \"()\" ,new caterwaul.syntax( \"dot\") .prefix( \" \") ,new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"a\") ,new caterwaul.syntax( \"b\") .prefix( \" \"))) ,new caterwaul.syntax( \"()\" ,new caterwaul.syntax( \"dot\") .prefix( \" \") ,new caterwaul.syntax( \",\" ,new caterwaul.syntax( \"b\") ,new caterwaul.syntax( \"b\") .prefix( \" \")))) .prefix( \" \")))))")} ;return(result)}) .call(this,new caterwaul.syntax( "(" ,new caterwaul.syntax( "function" ,new caterwaul.syntax( "(" ,new caterwaul.syntax( "_formals")) .prefix( " ") ,new caterwaul.syntax( "{" ,new caterwaul.syntax( "return" ,new caterwaul.syntax( "_e") .prefix( " "))) .prefix( " "))) ,new caterwaul.syntax( "," ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "b") .prefix( " ")) ,new caterwaul.syntax( "[" ,new caterwaul.syntax( "x")) ,new caterwaul.syntax( "," ,new caterwaul.syntax( "x") ,new caterwaul.syntax( "y") .prefix( " ")) ,new caterwaul.syntax( "+" ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "i")) ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "b") .prefix( " ") ,new caterwaul.syntax( "i"))) .prefix( " ") ,new caterwaul.syntax( "," ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "b") .prefix( " ")) ,new caterwaul.syntax( "[" ,new caterwaul.syntax( "x")) ,new caterwaul.syntax( "," ,new caterwaul.syntax( "x") ,new caterwaul.syntax( "y") .prefix( " ")) ,new caterwaul.syntax( "*" ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "i")) ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "b") .prefix( " ") ,new caterwaul.syntax( "i"))) .prefix( " ") ,new caterwaul.syntax( "," ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "b") .prefix( " ")) ,new caterwaul.syntax( "[" ,new caterwaul.syntax( "x")) ,new caterwaul.syntax( "," ,new caterwaul.syntax( "x") ,new caterwaul.syntax( "y") .prefix( " ")) ,new caterwaul.syntax( "-" ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "i")) ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "b") .prefix( " ") ,new caterwaul.syntax( "i"))) .prefix( " ") ,new caterwaul.syntax( "," ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "b") .prefix( " ")) ,new caterwaul.syntax( "[" ,new caterwaul.syntax( "x")) ,new caterwaul.syntax( "," ,new caterwaul.syntax( "x") ,new caterwaul.syntax( "y") .prefix( " ")) ,new caterwaul.syntax( "*" ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "i")) ,new caterwaul.syntax( "b") .prefix( " ")) .prefix( " ") ,new caterwaul.syntax( "," ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "b") .prefix( " ")) ,new caterwaul.syntax( "x") ,new caterwaul.syntax( "+" ,new caterwaul.syntax( "x") ,new caterwaul.syntax( "y") .prefix( " ")) .prefix( " ") ,new caterwaul.syntax( "*" ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "i")) ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "b") .prefix( " ") ,new caterwaul.syntax( "i"))) .prefix( " ") ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "()" ,new caterwaul.syntax( "." ,new caterwaul.syntax( "Math") ,new caterwaul.syntax( "sqrt")) ,new caterwaul.syntax( "x")) ,new caterwaul.syntax( "+" ,new caterwaul.syntax( "x") ,new caterwaul.syntax( "y") .prefix( " ")) .prefix( " ") ,new caterwaul.syntax( "*" ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "i")) ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "a") .prefix( " ") ,new caterwaul.syntax( "i"))) .prefix( " ") ,new caterwaul.syntax( "," ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "b") .prefix( " ")) ,new caterwaul.syntax( "[" ,new caterwaul.syntax( "x")) ,new caterwaul.syntax( "," ,new caterwaul.syntax( "x") ,new caterwaul.syntax( "y") .prefix( " ")) ,new caterwaul.syntax( "()" ,new caterwaul.syntax( "." ,new caterwaul.syntax( "Math") ,new caterwaul.syntax( "min")) ,new caterwaul.syntax( "," ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "i")) ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "b") .prefix( " ") ,new caterwaul.syntax( "i")))) ,new caterwaul.syntax( "," ,new caterwaul.syntax( "," ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "b") .prefix( " ")) ,new caterwaul.syntax( "c") .prefix( " ")) ,new caterwaul.syntax( "[" ,new caterwaul.syntax( "x")) ,new caterwaul.syntax( "," ,new caterwaul.syntax( "x") ,new caterwaul.syntax( "y") .prefix( " ")) ,new caterwaul.syntax( "+" ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "i")) ,new caterwaul.syntax( "*" ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "b") .prefix( " ") ,new caterwaul.syntax( "i")) ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "c") .prefix( " ") ,new caterwaul.syntax( "i"))) .prefix( " ")) .prefix( " ") ,new caterwaul.syntax( "," ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "b") .prefix( " ")) ,new caterwaul.syntax( "[" ,new caterwaul.syntax( "x")) ,new caterwaul.syntax( "," ,new caterwaul.syntax( "x") ,new caterwaul.syntax( "y") .prefix( " ")) ,new caterwaul.syntax( "()" ,new caterwaul.syntax( "." ,new caterwaul.syntax( "Math") ,new caterwaul.syntax( "max")) ,new caterwaul.syntax( "," ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "i")) ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "b") .prefix( " ") ,new caterwaul.syntax( "i")))) ,new caterwaul.syntax( "," ,new caterwaul.syntax( "," ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "b") .prefix( " ")) ,new caterwaul.syntax( "c") .prefix( " ")) ,new caterwaul.syntax( "[" ,new caterwaul.syntax( "x")) ,new caterwaul.syntax( "," ,new caterwaul.syntax( "x") ,new caterwaul.syntax( "y") .prefix( " ")) ,new caterwaul.syntax( "+" ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "i")) ,new caterwaul.syntax( "*" ,new caterwaul.syntax( "b") .prefix( " ") ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "c") .prefix( " ") ,new caterwaul.syntax( "i"))) .prefix( " ")) .prefix( " ") ,new caterwaul.syntax( "," ,new caterwaul.syntax( "," ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "b") .prefix( " ")) ,new caterwaul.syntax( "c") .prefix( " ")) ,new caterwaul.syntax( "[" ,new caterwaul.syntax( "x")) ,new caterwaul.syntax( "," ,new caterwaul.syntax( "x") ,new caterwaul.syntax( "y") .prefix( " ")) ,new caterwaul.syntax( "+" ,new caterwaul.syntax( "*" ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "i")) ,new caterwaul.syntax( "(" ,new caterwaul.syntax( "-" ,new caterwaul.syntax( "1") ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "b") .prefix( " ") ,new caterwaul.syntax( "i"))) .prefix( " ")) .prefix( " ")) .prefix( " ") ,new caterwaul.syntax( "*" ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "c") .prefix( " ") ,new caterwaul.syntax( "i")) ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "b") .prefix( " ") ,new caterwaul.syntax( "i"))) .prefix( " ")) .prefix( " ") ,new caterwaul.syntax( "," ,new caterwaul.syntax( "," ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "b") .prefix( " ")) ,new caterwaul.syntax( "c") .prefix( " ")) ,new caterwaul.syntax( "[" ,new caterwaul.syntax( "x")) ,new caterwaul.syntax( "," ,new caterwaul.syntax( "x") ,new caterwaul.syntax( "y") .prefix( " ")) ,new caterwaul.syntax( "+" ,new caterwaul.syntax( "*" ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "i")) ,new caterwaul.syntax( "(" ,new caterwaul.syntax( "-" ,new caterwaul.syntax( "1") ,new caterwaul.syntax( "b") .prefix( " ")) .prefix( " ")) .prefix( " ")) .prefix( " ") ,new caterwaul.syntax( "*" ,new caterwaul.syntax( "[]" ,new caterwaul.syntax( "c") .prefix( " ") ,new caterwaul.syntax( "i")) ,new caterwaul.syntax( "b") .prefix( " ")) .prefix( " ")) .prefix( " ") ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "()" ,new caterwaul.syntax( "scale") ,new caterwaul.syntax( "," ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "/" ,new caterwaul.syntax( "1.0") .prefix( " ") ,new caterwaul.syntax( "()" ,new caterwaul.syntax( "norm") .prefix( " ") ,new caterwaul.syntax( "a"))) .prefix( " "))) ,new caterwaul.syntax( "," ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "b") .prefix( " ")) ,new caterwaul.syntax( "()" ,new caterwaul.syntax( "scale") ,new caterwaul.syntax( "," ,new caterwaul.syntax( "b") ,new caterwaul.syntax( "/" ,new caterwaul.syntax( "()" ,new caterwaul.syntax( "dot") .prefix( " ") ,new caterwaul.syntax( "," ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "b") .prefix( " "))) ,new caterwaul.syntax( "()" ,new caterwaul.syntax( "dot") .prefix( " ") ,new caterwaul.syntax( "," ,new caterwaul.syntax( "b") ,new caterwaul.syntax( "b") .prefix( " ")))) .prefix( " "))) ,new caterwaul.syntax( "," ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "b") .prefix( " ")) ,new caterwaul.syntax( "()" ,new caterwaul.syntax( "minus") ,new caterwaul.syntax( "," ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "()" ,new caterwaul.syntax( "scale") .prefix( " ") ,new caterwaul.syntax( "," ,new caterwaul.syntax( "b") ,new caterwaul.syntax( "/" ,new caterwaul.syntax( "()" ,new caterwaul.syntax( "dot") .prefix( " ") ,new caterwaul.syntax( "," ,new caterwaul.syntax( "a") ,new caterwaul.syntax( "b") .prefix( " "))) ,new caterwaul.syntax( "()" ,new caterwaul.syntax( "dot") .prefix( " ") ,new caterwaul.syntax( "," ,new caterwaul.syntax( "b") ,new caterwaul.syntax( "b") .prefix( " ")))) .prefix( " "))))))) ;