(*  vector signal processing

    IIR tests

    © Philippe Strauss, 2012  *)


open Oneliners
include T_common
open Vec_sp
open Bilinear
open S_poles
open Filter_process
open Complex


let fcl = 600.
let fch = 400.
let fs = 44100
let ffs = foi fs
let alen = 8192


let () =

    let wcl_dt = omega_warp fcl ffs in
    let wcl_ct = omega_dt2ct fcl ffs in
    let kwcl = wcl_ct /. (fcl *. 2. *. pi) in
    Printf.printf "wcl_dt=%f; wcl_ct=%f [rad/s]; wcl_ct/(2*fc*pi)=%f\n" wcl_dt wcl_ct kwcl ;

    let wch_dt = omega_warp fch ffs in
    let wch_ct = omega_dt2ct fch ffs in
    let kwch = wch_ct /. (fch *. 2. *. pi) in
    Printf.printf "wch_dt=%f; wch_ct=%f [rad/s]; wch_ct/(2*fc*pi)=%f\n" wch_dt wch_ct kwch ;

	let p0 = butterworth_pole 12 0 in
	let p1 = butterworth_pole 12 1 in
	let p2 = butterworth_pole 12 2 in
	let p3 = butterworth_pole 12 3 in
	let p4 = butterworth_pole 12 4 in
	let p5 = butterworth_pole 12 5 in

    let clp0 = allpoles_lp wcl_dt p0 in
    let clp1 = allpoles_lp wcl_dt p1 in
    let clp2 = allpoles_lp wcl_dt p2 in
    let clp3 = allpoles_lp wcl_dt p3 in
    let clp4 = allpoles_lp wcl_dt p4 in
    let clp5 = allpoles_lp wcl_dt p5 in

    let chp0 = allpoles_hp wch_dt p0 in
    let chp1 = allpoles_hp wch_dt p1 in
    let chp2 = allpoles_hp wch_dt p2 in
    let chp3 = allpoles_hp wch_dt p3 in
    let chp4 = allpoles_hp wch_dt p4 in
    let chp5 = allpoles_hp wch_dt p5 in

    let z0_lp = Array.make 3 0. in
    let z1_lp = Array.make 3 0. in
    let z2_lp = Array.make 3 0. in
    let z3_lp = Array.make 3 0. in
    let z4_lp = Array.make 3 0. in
    let z5_lp = Array.make 3 0. in

    let z0_hp = Array.make 3 0. in
    let z1_hp = Array.make 3 0. in
    let z2_hp = Array.make 3 0. in
    let z3_hp = Array.make 3 0. in
    let z4_hp = Array.make 3 0. in
    let z5_hp = Array.make 3 0. in

    let samp = Array.make alen 0. in
    samp.(0) <- 1.0 ;

    iir2_df2 (fst clp0) (invsign (snd clp0)) z0_lp samp alen ;
    iir2_df2 (fst clp1) (invsign (snd clp1)) z1_lp samp alen ;
    iir2_df2 (fst clp2) (invsign (snd clp2)) z2_lp samp alen ;
    iir2_df2 (fst clp3) (invsign (snd clp3)) z3_lp samp alen ;
    iir2_df2 (fst clp4) (invsign (snd clp4)) z4_lp samp alen ;
    iir2_df2 (fst clp5) (invsign (snd clp5)) z5_lp samp alen ;

    iir2_df2 (fst chp0) (invsign (snd chp0)) z0_hp samp alen ;
    iir2_df2 (fst chp1) (invsign (snd chp1)) z1_hp samp alen ;
    iir2_df2 (fst chp2) (invsign (snd chp2)) z2_hp samp alen ;
    iir2_df2 (fst chp3) (invsign (snd chp3)) z3_hp samp alen ;
    iir2_df2 (fst chp4) (invsign (snd chp4)) z4_hp samp alen ;
    iir2_df2 (fst chp5) (invsign (snd chp5)) z5_hp samp alen ;

    plot samp fs Time ;
    plot samp fs Frequency ;



