(*  vector signal processing

    IIR tests

    © Philippe Strauss, 2012  *)


open Oneliners
open Vec_sp
open Lti_xfer
open Complex
open S_poles
open Bilinear
module P = Gnuplot.Array


let () =

    let wn = Array.init 1400 (fun x -> foi (x-700) /. 100.) in

    let cminusone = {re=(-1.); im=0.} in
    let numerator = [cminusone; cminusone; cminusone; cminusone; cminusone; cminusone; cminusone] in
    (* let numerator = [cminusone;] in *)

    let cheby1_pole7 = cheby1_pole 1.0 7 in
    let py0 = cheby1_pole7 0 in
    let py1 = cheby1_pole7 1 in
    let py2 = cheby1_pole7 2 in
    let py3 = cheby1_pole7 3 in
    pr "\nS domain normalized 7th order Chebychev I poles location: re0=%f, im0=%f, re1=%f, im1=%f, re2=%f, im2=%f, re3=%f, im3=%f\n" py0.re py0.im py1.re py1.im py2.re py2.im py3.re py3.im;
    let wc_dt = omega_warp 4000. 44100. in
    let pz_bilin_ch1_7 = pz_bilin wc_dt in
    let pyw0 = pz_bilin_ch1_7 py0 in
    let pyw1 = pz_bilin_ch1_7 py1 in
    let pyw2 = pz_bilin_ch1_7 py2 in
    let pyw3 = pz_bilin_ch1_7 py3  in
    let denominator = [pyw0; (Complex.conj pyw0); pyw1; (Complex.conj pyw1); pyw2; (Complex.conj pyw2); pyw3] in
    let norm = norm2_lp wc_dt in
    let g0 = (norm py0)
    and g1 = (norm py1)
    and g2 = (norm py2)
    and g3 = (norm1_lp wc_dt py3.re) in
    let g = g0*.g1*.g2*.g3 in
    pr "wc_dt=%f; b0/a0 gain normalization: %f %f %f %f, overall=%e\n" wc_dt g0 g1 g2 g3 g;

    let lti_cheby1_7_xfer = xfer_prod1 numerator denominator ~g0:g in
    let lti_cheby1_7_phase = parg_prod1 numerator denominator in
    let lti_cheby1_7_grd = grd_prod1 numerator denominator 0. in
    let lti_cheby1_7_grd' = grd'_prod1 numerator denominator 0. in

    let hjw = Array.map (fun w -> lti_cheby1_7_xfer w) wn in
    let hmag = cpower_db hjw in
    let ph = Array.map (fun w -> lti_cheby1_7_phase w) wn in
    let grd = Array.map (fun w -> lti_cheby1_7_grd w) wn in
    let grd' = Array.map (fun w -> lti_cheby1_7_grd' w) wn in

    let g = P.init ~xsize:800. ~ysize:600. Gnuplot.X in
    P.env g (0.) (pi) (-80.) (10.) ;
    P.pen g 1 ;
    P.xy g wn hmag ;
    P.close g ;

    let g = P.init ~xsize:800. ~ysize:600. Gnuplot.X in
    P.env g (0.) (pi) (-.4.*.pi) (0.) ;
    P.pen g 2 ;
    P.xy g wn ph ;
    P.close g ;

    let g = P.init ~xsize:800. ~ysize:600. Gnuplot.X in
    P.env g (0.) (pi) (0.) (60.) ;
    P.pen g 3 ;
    P.xy g wn grd ;
    P.close g ;

    let g = P.init ~xsize:800. ~ysize:600. Gnuplot.X in
    P.env g (0.) (pi) (-1500.) (1500.) ;
    P.pen g 4 ;
    P.xy g wn grd' ;
    P.close g

