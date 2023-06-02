(*  OCaml short names handy for numcomp

    © Philippe Strauss, 2009 - 2012  *)


let pi = 4. *. atan 1.

let pr = Printf.printf
let spr = Printf.sprintf

let foi = float_of_int
let iof = int_of_float

let soi = string_of_int
let ios = int_of_string

let floori x = iof (floor x)
let ceili x = iof (ceil x)

let arlen = Array.length
let balen = Bigarray.Array1.dim

let hfind = Hashtbl.find

let half_down x =
    if (x mod 2) = 0 then x / 2 else (x-1) / 2

let half_up x =
    if (x mod 2) = 0 then x / 2 else (x+1) / 2

(*  TESTED 2013-04-16  w. bessel filters s polynomial  *)
let pow n exp =
	assert (n >= 0 && exp >= 0) ;
    let rec loop i = if i = 0 then 1 else n*loop (i-1) in
    loop exp

(*  TESTED 2013-04-16  w. bessel filters s polynomial  *)
let rec fact n =
	assert (n>=0) ;
    if n < 2 then 1 else n * fact (n-1)

exception NegativeNotApplicable

let log10_ex x = if x > 0. then log10 x else raise NegativeNotApplicable

let polar_from_complex z =
    (Complex.norm z, Complex.arg z)

let complex_from_polar = Complex.polar

(*  TESTED 2013-04-14  w. chebychev filters s poles  *)
let asinh x = log (x +. sqrt (x**2. +. 1.))

(*  TESTED 2013-04-14  w. chebychev filters s poles  *)
let acosh x =
    assert (x >= 1.) ;
    log (x +. sqrt (x**2. -. 1.))

(* let spr_farray floats =
    BatPrintf.sprintf2 "%a" (BatArray.print ~sep:"; " BatFloat.print) floats

let prf floats = Printf.printf "%s\n" (spr_farray floats) *)

let spr_complex z =
	let radius, arc = polar_from_complex z in
	let degr = (arc /. pi) *. 180. in
	Printf.sprintf "re = %+e, im = %+e; radius = %e, arc = %+e = %+f°" Complex.(z.re) Complex.(z.im) radius arc degr

let spr_conj z =
	let radius, arc = polar_from_complex z in
	let degr = (arc /. pi) *. 180. in
	Printf.sprintf "re = %+e, im = +/- % e; radius= %e, arc = % e = %+f°" Complex.(z.re) Complex.(z.im) radius arc degr

let spr_angle arc =
	let degr = (arc /. pi) *. 180. in
	Printf.sprintf "<| %+e [rad] = %+f°" arc degr

