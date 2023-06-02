(*  one-dimensional discrete-time signal processing, with audio in mind.

    IIR filters design/implementation helpers : Bilinear transform

    © Philippe Strauss, 2011 - 2013
    License: see accompanying LICENSE.txt or LICENSE.md file.

    References :

        [OppSch]    : Discrete-Time Signal Processing, Oppenheim & Schafer, Prentice Hall, 3rd ed.
        [ProMan]    : Digital Signal Processing, Proakis & Manolakis, Macmillan, 2nd ed.
        [OppWill]   : Signals & Systems, Oppenheim & Willsky, Prentice Hall, Intl 2nd ed.
        [ElFilt]    : Electronic Filter Design Handbook, 4th Ed.
        [IngPro]    : Digital Signal Processing using MATLAB, Ingle & Proakis  *)


open Complex
open Oneliners


(*  there are two ways in the literature to define bilinear warping :
    2. *. tan (w/.2.) or tan (w/.2.)
    using the second imply valuating 2/Ts to one *)
let omega_warp fc fs =
    tan (fc *. pi /. fs)

(*  Discrete time to continuous time warping of omega *)
let omega_dt2ct fc fs =
    2. *. fs *. omega_warp fc fs
 

(*  Lowpass all poles bilinear biquads z coefficients from s plane conjugate pole coeffs
    wc:  Tangent warped, discrete time domain pulsation in normalized radian per second (Fs/.2. = pi)
         See [ElFilt] p. 11 eq. 2-2 about the basics, wc is nammed FSF in the text.

    All-pole lowpass to (Wc) lowpass S domain freq. xform embedded,
    plus bilinear s = f(z^(-1)) variable substitution (hand/head crafted :-)
    TESTED 2013-04-15 in t_bilin_lr_impulse.ml and t_bilin_butter_impulse.ml  *)
let allpoles_lp wc pz =
    let pn2 = Complex.norm2 pz in
    let wc2pn2 = wc ** 2. *. pn2 in
    let s1 = 2. *. wc *. pz.re in
    let d2 = 1. +. s1 +. wc2pn2
    and d1 = 2. *. (wc2pn2 -. 1.)
    and d0 = 1. -. s1 +. wc2pn2 in
    let n2 = wc2pn2
    and n1 = 2. *. wc2pn2
    and n0 = wc2pn2 in
    ([|n0/.d0; n1/.d0; n2/.d0|], [|1.; d1/.d0; d2/.d0|])

(*  TESTED 2013-04-15 in t_bilin_lr_impulse.ml and t_bilin_butter_impulse.ml *)
let allpoles_hp wc pz =
    let pn2 = Complex.norm2 pz in
    let wc2 = wc ** 2. in
    let s1 = 2. *. wc *. pz.re in
    let d2 = pn2 +. s1 +. wc2
    and d1 = 2. *. (wc2 -. pn2)
    and d0 = pn2 -. s1 +. wc2 in
    let n2 = pn2
    and n1 = (-.2.) *. pn2
    and n0 = pn2 in
    ([|n0/.d0; n1/.d0; n2/.d0|], [|1.; d1/.d0; d2/.d0|])

(*  b0/a0 overall gain normalization for when using only polar conjugate lti_xfer functions
    = biquad numerator0/denominator0. TESTED 2013-04-15 against chebychev1_7 test
    pz is S domain unitary pulsation pole  *)
let norm2_lp wc pz =
    let pn2 = Complex.norm2 pz in
    let wc2pn2 = wc ** 2. *. pn2 in
    let res = wc2pn2 /. (1. -. 2. *. wc *. pz.re +. wc2pn2) in
    assert (res > 0.) ;
    res

let norm2_hp wc pz =
    let pn2 = Complex.norm2 pz in
    let res = pn2 /. (pn2 -. 2. *. wc *. pz.re +. wc ** 2.) in
    assert (res > 0.) ;
    res    

(*  Single real S pole gain normalization for lowpass
    TESTED 2013-04-15 against chebychev1_7 test  *)
let norm1_lp wc pr =
    let res = (-1. *. wc *. pr) /. (1.-.wc*.pr) in
    assert (res > 0.) ;
    res


(*  bilinear transform for a single pole conjugate pair kept in polar form  
    TESTED in t_bilin_polar.ml 2013-04-12  *)
let pz_bilin w_dt s_pz =
    let sigma = s_pz.re *. w_dt in
    let omega = s_pz.im *. w_dt in
    (Complex.div {re = 1.+.sigma ; im = omega} {re = 1.-.sigma ; im = (-.omega)})

(*  TESTED in t_bilin_polar.ml 2013-04-12  *)
let pz_bilin_polar w_dt s_pz =
    let z = pz_bilin s_pz w_dt in
    polar_from_complex z
