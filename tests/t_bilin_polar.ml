(*  vector signal processing

    IIR tests

    © Philippe Strauss, 2012  *)


open Oneliners
include T_common
open Vec_sp
open Complex
open Lti_xfer
open Bilinear
open S_poles


let fc = 2000.
let fs = 44100
let ffs = foi fs


let () =

    let wc_dt = omega_warp fc ffs in
    let wc_ct = omega_dt2ct fc ffs in
    let kwc = wc_ct /. (fc *. 2. *. pi) in
    pr "\nwc_dt=%f [normalized rad/s rel. Fs/2]; wc_ct=%f [rad/s]; warping factor : wc_ct/(2*fc*pi)=%f\n" wc_dt wc_ct kwc ;

	let p0 = butterworth_pole 4 0 in
	let p1 = butterworth_pole 4 1 in
	pr "\nS domain normalized to unity pulsation Butterworth 4th order poles: p0 : %s, p1 : %s\n" (spr_conj p0) (spr_conj p1);

	let r0, arc0 = pz_bilin_polar p0 wc_dt in
	let r1, arc1 = pz_bilin_polar p1 wc_dt in
	pr  "\npolar z domain conjugate poles:\n\nr0 = %f, arc0 = %f, r1 = %f, arc1 = %f\n" r0 arc0 r1 arc1 ;

	let g0 = norm2_lp wc_dt p0 in
	let g1 = norm2_lp wc_dt p1 in
	pr "\ng0=%f; g1=%f\n" g0 g1 ;

	let z00, z01, z02 = poly2_from_conjugate r0 arc0 in
	let z10, z11, z12 = poly2_from_conjugate r1 arc1 in
	pr "\npoly2_from_conjugate :\n\npoly0 = (%f, %f, %f)\npoly1 = (%f, %f, %f)\n" z00 z01 z02 z10 z11 z12 ;

	let b_num0, a_denom0 = allpoles_lp wc_dt p0 in
	let b_num1, a_denom1 = allpoles_lp wc_dt p1 in
	pr "\nb/a polynomials :\n\nb0 = %s; a0 = %s;\nb1 = %s; a1 = %s\n\n" (spr_farray b_num0) (spr_farray a_denom0) (spr_farray b_num1) (spr_farray a_denom1) ;

	begin match conjugate_from_poly2 (z00, z01, z02) with
	| Conjugate pole_c -> (let r, arc = polar_from_complex pole_c in Printf.printf "conjugate_from_poly2 poly0 = %s\n" (spr_conj pole_c))
	| Real pole_r -> (Printf.printf "conjugate_from_poly2 poly0: single real root at %f\n" pole_r) end ;

	let thb0 = butterworth_pole_arc 7 0
	and thb1 = butterworth_pole_arc 7 1
	and thb2 = butterworth_pole_arc 7 2
	and thb3 = butterworth_pole_arc 7 3 in
	pr "\nButterworth order 7 angles: 0:%s, 1:%s 2:%s, 3:%s\n" (spr_angle thb0) (spr_angle thb1) (spr_angle thb2) (spr_angle thb3) ;

	let ch2z0 = cheby2_zero 7 0
	and ch2z1 = cheby2_zero 7 1
	and ch2z2 = cheby2_zero 7 2
	(* and ch2z3 = cheby2_zero 7 3 *) in
	pr "\nChebychev type II order 7 zeroes : 1:%s, 2:%s, 3:%s\n" (spr_conj ch2z0) (spr_conj ch2z1) (spr_conj ch2z2) (* spr_conj ch2z3 *)
