(*  vector signal processing

    IIR tests

    © Philippe Strauss, 2012  *)


open Oneliners
include T_common
open Vec_sp
open Lti_xfer
open Complex


let () =

    (*  ---------------------------------
        OppSch 2nd ed., p. 457, eq. 7.37
        LP Butterworth example of order 6
        ---------------------------------s

        Numerator:
        ----------
        >> roots([1,1])
        ans = -1

        (with g0 = 0.0007378)

        Denominator:
        ------------
        >> roots([1, -1.2686, 0.7051])
        ans =
           0.6343 + 0.5502i
           0.6343 - 0.5502i

        >> roots([1, -1.0106, 0.3583])
        ans =
           0.5053 + 0.3209i
           0.5053 - 0.3209i

        >> roots([1, -0.9044, 0.2155])
        ans =
           0.4522 + 0.1050i
           0.4522 - 0.1050i  *)

    let wn = Array.init 1400 (fun x -> foi (x-700) /. 100.) in

    let g0 = 0.0007378 in
    let cminusone = {re=(-1.); im=0.} in
    let numerator = [cminusone; cminusone; cminusone; cminusone; cminusone; cminusone] in
    let denominator = [ {re=0.6343; im=0.5502}; {re=0.6343; im=(-0.5502)};
                        {re=0.5053; im=0.3209}; {re=0.5053; im=(-0.3209)};
                        {re=0.4522; im=0.1050}; {re=0.4522; im=(-0.1050)} ] in

    let lti_butter_lp6_xfer = xfer_prod1 numerator denominator ~g0:g0 in
    let lti_butter_lp6_phase = parg_prod1 numerator denominator in
    let lti_butter_lp6_grd = grd_prod1 numerator denominator 0. in
    let lti_butter_lp6_grd' = grd'_prod1 numerator denominator 0. in

    let hjw = Array.map (fun w -> lti_butter_lp6_xfer w) wn in
    let hmag = cpower_db hjw in
    let ph = Array.map (fun w -> lti_butter_lp6_phase w) wn in
    let grd = Array.map (fun w -> lti_butter_lp6_grd w) wn in
    let grd' = Array.map (fun w -> lti_butter_lp6_grd' w) wn in

    let g = P.init ~xsize:800. ~ysize:600. Gnuplot.X in
    P.env g (0.) (pi) (-80.) (10.) ;
    P.pen g 1 ;
    P.xy g wn hmag ;
    P.close g ;

    let g = P.init ~xsize:800. ~ysize:600. Gnuplot.X in
    P.env g (0.) (pi) (-.3.*.pi) (0.) ;
    P.pen g 2 ;
    P.xy g wn ph ;
    P.close g ;

    let g = P.init ~xsize:800. ~ysize:600. Gnuplot.X in
    P.env g (0.) (pi) (0.) (15.) ;
    P.pen g 3 ;
    P.xy g wn grd ;
    P.close g ;

    let g = P.init ~xsize:800. ~ysize:600. Gnuplot.X in
    P.env g (0.) (pi) (-30.) (30.) ;
    P.pen g 4 ;
    P.xy g wn grd' ;
    P.close g

