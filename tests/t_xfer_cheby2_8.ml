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

    (* let cminusone = {re=(-1.); im=0.} in *)
    (* let numerator = [cminusone; cminusone; cminusone; cminusone; cminusone; cminusone; cminusone] in *)
    (* let numerator = [cminusone;] in *)

    let cheby2_pole8 = cheby2_pole (-.40.) 7 in
    let py0 = cheby2_pole8 0
    and py1 = cheby2_pole8 1
    and py2 = cheby2_pole8 2
    and py3 = cheby2_pole8 3 in
    let wc_dt = omega_warp 4000. 44100. in
    let pz_bilin_wc = pz_bilin wc_dt in
    let pyw0 = pz_bilin_wc py0
    and pyw1 = pz_bilin_wc py1
    and pyw2 = pz_bilin_wc py2
    and pyw3 = pz_bilin_wc py3 in
    let denominator = [pyw0; (Complex.conj pyw0); pyw1; (Complex.conj pyw1); pyw2; (Complex.conj pyw2); pyw3; (*Complex.conj pyw3*)] in

    let cheby2_zero_8 = cheby2_zero 7 in
    let zy0 = cheby2_zero_8 0
    and zy1 = cheby2_zero_8 1
    and zy2 = cheby2_zero_8 2
    (* and zy3 = cheby2_zero_8 3 *) in
    let zyw0 = pz_bilin_wc zy0
    and zyw1 = pz_bilin_wc zy1
    and zyw2 = pz_bilin_wc zy2
    (* and zyw3 = pz_bilin_wc zy3 *) in
    let numerator = [zyw0; (Complex.conj zyw0); zyw1; (Complex.conj zyw1); zyw2; (Complex.conj zyw2); (* zyw3; (Complex.conj zyw3) *)] in

    let norm = norm2_lp wc_dt in
    let gd0 = (norm py0)
    and gd1 = (norm py1)
    and gd2 = (norm py2)
    (* and gd3 = (norm py3) in *)
    and gd3 = norm1_lp wc_dt py3.re in
    let gn0 = (norm zy0)
    and gn1 = (norm zy1)
    and gn2 = (norm zy2)
    (* and gn3 = (norm zy3) *) in
    let g = (gd0*.gd1*.gd2*.gd3) /. (gn0*.gn1*.gn2(* *.gn3 *)) in
    pr "wc_dt=%f; b0/a0 gain normalization: %e\n" wc_dt g;

    let lti_cheby2_8_xfer = xfer_prod1 numerator denominator ~g0:g in
    let lti_cheby2_8_phase = parg_prod1 numerator denominator in
    let lti_cheby2_8_grd = grd_prod1 numerator denominator 0. in
    let lti_cheby2_8_grd' = grd'_prod1 numerator denominator 0. in

    let hjw = Array.map (fun w -> lti_cheby2_8_xfer w) wn in
    let hmag = cpower_db hjw in
    let ph = Array.map (fun w -> lti_cheby2_8_phase w) wn in
    let grd = Array.map (fun w -> lti_cheby2_8_grd w) wn in
    let grd' = Array.map (fun w -> lti_cheby2_8_grd' w) wn in

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
