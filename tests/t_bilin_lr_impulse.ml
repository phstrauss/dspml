(*  vector signal processing

    IIR tests

    © Philippe Strauss, 2011, 2012 *)


open Oneliners
include T_common
open Bigarray
open Complex
open Bilinear
open S_poles
open Filter_process


let fc = 500.
let fs = 44100
let ffs = foi fs
let alen = 4096


let () =

    let wc_dt = omega_warp fc ffs in
    let wc_ct = omega_dt2ct fc ffs in
    let kwc = wc_ct /. (fc *. 2. *. pi) in
    Printf.printf "wc_dt=%f; wc_ct=%f [rad/s]; wc_ct/(2*fc*pi)=%f\n" wc_dt wc_ct kwc ;
 
    let p2 = butterworth_pole 2 0 in
    let coeffs_lp = allpoles_lp wc_dt p2 in

    let samp_lp = Array.make alen 0. in
    let z0_lp = Array.make 3 0. in
    let z1_lp = Array.make 3 0. in
    samp_lp.(0) <- 1.0 ;

    (* see [IngPro] p. 48 about the next line *)
    (* iir2_df2 [|1.; 0.; 0.|] [|-1.; 1.; -0.9|] z0 tsamples alen ; *)

    (* matlab butter wc = 0.1 * pi *)
    (* iir2_df2 [|0.0675; 0.1349; 0.0675|] (invsign [|1.0; -1.1430; 0.4128|]) z0 tsamples alen ; *)

    iir2_df1 (fst coeffs_lp) (invsign (snd coeffs_lp)) z0_lp z1_lp samp_lp alen ;
    iir2_df1 (fst coeffs_lp) (invsign (snd coeffs_lp)) z0_lp z1_lp samp_lp alen ;

    plot samp_lp fs Time ;
    plot samp_lp fs Frequency ;

    let coeffs_hp = allpoles_hp wc_dt p2 in

    let samp_hp = Array.make alen 0. in
    let z0_hp = Array.make 3 0. in
    let z1_hp = Array.make 3 0. in
    samp_hp.(0) <- 1.0 ;

    iir2_df1 (fst coeffs_hp) (invsign (snd coeffs_hp)) z0_hp z1_hp samp_hp alen ;
    iir2_df1 (fst coeffs_hp) (invsign (snd coeffs_hp)) z0_hp z1_hp samp_hp alen ;

    plot samp_hp fs Time ;
    plot samp_hp fs Frequency
