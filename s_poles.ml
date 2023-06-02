(*  one-dimensional discrete-time signal processing, with audio in mind.

    IIR filters design/implementation helpers

    © Philippe Strauss, 2011 - 2013
    License: see accompanying LICENSE.txt or LICENSE.md file.

    References :

        [OppSch]    : Discrete-Time Signal Processing, Oppenheim & Schafer, Prentice Hall, 3rd ed.
        [ProMan]    : Digital Signal Processing, Proakis & Manolakis, Macmillan, 2nd ed.
        [OppWill]   : Signals & Systems, Oppenheim & Willsky, Prentice Hall, Intl 2nd ed.
        [ElFilt]    : Electronic Filter Design Handbook, 4th Ed.
        [IngPro]    : Digital Signal Processing using MATLAB, Ingle & Proakis  *)


open Oneliners
(* open Zroots *)


(*  About IIR ordering :

VV : Unless you can check through all possible arrangements, common practice
is put sections in the order of decreasing Q, keep together nearest
poles and zeroes, alternate HPF and LPF type of sections. The result
would be generally not too far from optimal.

RBJ or DW : The alternate HPF/LPF doesn't happen while there are remaining zeros not
on either z=1 or z=-1.  then they're BPFs.  but i agree, start with the
least stable pair of poles and match the zeros that are closest to them.   *)


(*  Common analog / Laplace transform domain IIR filters characteristics functions  *)


(*  BUTTERWORTH poles & magnitude  *)

(*  see [OppWill] p. 704, eq. 9.147 - TESTED *)
let butterworth_pole_arc order m =
    assert(2*m < order) ;
    pi *. (2. *. (foi m) +. 1.) /. (2. *. (foi order)) +. (pi /. 2.)

(*  polar to ortho w. omega_cutoff scaling
    see [OppWill] p. 704, eq. 9.148 - TESTED *)
let butterworth_pole ?kwc:(kw=1.0) order m =
    let open Complex in
    Complex.mul {re=kw; im=0.0} (Complex.exp {re=0.0; im=(butterworth_pole_arc order m)})

(* max flat magnitude function *)
let denom_butterworth_mag f fc order =
    sqrt (1. +. (f /. fc) ** (2. *. (foi order)))

let butterworth_mag_lp f fc order =
    1. /. (denom_butterworth_mag f fc order)

let butterworth_mag_hp f fc order =
    (* numerator: from the top of my head, CHECK ! *)
    ((f /. fc) ** (foi order)) /. denom_butterworth_mag f fc order


(*  LINKWITZ-RILEY magnitude *)

let lr_mag_lp f fc order =
    assert(order mod 2 = 0);
    (butterworth_mag_lp f fc (order/2)) ** 2.

let lr_mag_hp f fc order =
    assert(order mod 2 = 0);
    (butterworth_mag_hp f fc (order/2)) ** 2.


(*  CHEBYCHEV type I & II poles/zeroes  *)

let cheby_theta_m order m =
    assert (m>0 && m<=order) ;
    (pi/.2.) *. ((2.*.(foi m)-.1.)/.(foi order))

let cheby_pole_raw_epsilon epsilon order m =
    let open Complex in
    let a = (1. /. (foi order)) *. (asinh (1./.epsilon)) in
    let kr = sinh a and ki = cosh a in
    let bpole = butterworth_pole order m in
    {re = bpole.re *. kr; im = bpole.im *. ki}

let cheby1_pole ripple_db order m =
    let cheby1_epsilon ripple =
        (sqrt (10. ** (ripple /. 10.) -. 1.)) in
    let epsilon = cheby1_epsilon ripple_db in
    cheby_pole_raw_epsilon epsilon order m

(*  odd order needs one more zero at infinity
    corresponding to a div by zero at theta_m = +/- pi
    here beyond assert fence *)
let cheby2_zero order m0 =
    let open Complex in
    let m = m0+1 in
    assert(m <= half_down order) ;
    let th = cheby_theta_m order m in
    Complex.inv {re = 0. ; im = -.cos th}

let cheby2_pole floor_db order m =
    let cheby2_epsilon floor =
        10. ** (floor /. 20.) in
    let epsilon = cheby2_epsilon floor_db in
    Complex.inv (cheby_pole_raw_epsilon epsilon order m)


(*  BESSEL s polynomials and poles  *)

(*  TESTED 2013-04-16 using order 5, against values on display on wikipedia bessel filter page  *)
let bessel_ak n k = fact (2*n - k) / ((pow 2 (n-k)) * (fact k) * fact (n-k))

(*  TESTED 2013-04-18  *)
let bessel_poly order =
    Array.init (order+1) (fun i -> foi (bessel_ak order i))

(*  FIX : suppress poles with neg im part (implicit conjugation)  *)
(* let bessel_poles order =
    Zroots.zroots (bessel_poly order) *)
