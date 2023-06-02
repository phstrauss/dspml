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


(* found in ADI app note *)
let linphase005deg_n8 = [|
	{re = -.0.8195 ; im = 0.3711} ;
	{re = -.0.7930 ; im = 1.1054} ;
	{re = -.0.7213 ; im = 1.8134} ;
	{re = -.0.5341 ; im = 2.4761} ;
|]


let () =

    let wc_dt = omega_warp fc ffs in
    let wc_ct = omega_dt2ct fc ffs in
    let kwc = wc_ct /. (fc *. 2. *. pi) in
    Printf.printf "\nwc_dt=%f [normalized rad/s rel. Fs/2]; wc_ct=%f [rad/s]; warping factor : wc_ct/(2*fc*pi)=%f\n" wc_dt wc_ct kwc ;

	(* ***** *)

	let lph8_polar = Array.map (fun c -> pz_polar_bilin c wc_dt) linphase005deg_n8 in

	let lph8_poly2 = Array.map (fun c -> lp_bilin_biquad wc_dt c) linphase005deg_n8 in ()
