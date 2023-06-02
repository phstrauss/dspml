(*  vector signal processing

    IIR tests

    © Philippe Strauss, 2012  *)


open Oneliners
include T_common
open Vec_sp
open Lti_xfer
open Complex


let () =

    let pz0 = Complex.mul {re=0.9; im=0.} (Complex.exp {re=0.; im=pi/.2.}) in
	let lti_1st = xfer_prod1 [pz0] [Complex.zero] in

	let wn = Array.init 1400 (fun x -> foi (x-700) /. 100.) in

	let hw1 = Array.map (fun w -> lti_1st w) wn in
	let hmag1 = cpower_db hw1 in

    let hw2 = Array.map (fun w -> xfer_factor1 0.9 (pi/.2.) w) wn in
    let hmag2 = cpower_db hw2 in

    let mag3 = Array.map (fun w -> mag_factor1 0.9 (pi/.2.) w) wn in

    let g = P.init ~xsize:800. ~ysize:600. Gnuplot.X in
    P.env g (* ~xlog:true *) (-7.) (7.) (-40.) (40.) ;
    P.pen g 1 ;
    P.xy g wn hmag1 ;
    P.pen g 2 ;
    P.xy g wn hmag2 ;
    P.pen g 3 ;
    P.xy g wn (log_db20 mag3) ;
    P.close g ;

    (*  MATLAB :

        see [IngPro] p. 48 :

        b=[1,0,0]; a=[1, -1, 0.9];

        >> roots([1,0,0])
        ans =
             0
             0

        >> roots([1, -1, 0.9])
        ans =
           0.5000 + 0.8062i
           0.5000 - 0.8062i  *)

    let lti_2nd_poly = xfer_poly [Complex.one; Complex.zero; Complex.zero] [Complex.one; {re=(-1.); im=0.}; {re=0.9; im=0.}] in

    let num = [Complex.zero; Complex.zero] in
    let pole = {re=0.5; im=0.8062} in
    let denom = [pole; (Complex.conj pole)] in
    let lti_2nd_prod = xfer_prod1 num denom in

    let hw4 = Array.map (fun w -> lti_2nd_poly w) wn in
    let hmag4 = cpower_db hw4 in

    let hw5 = Array.map (fun w -> lti_2nd_prod w) wn in
    let hmag5 = cpower_db hw5 in

    let r, arc = polar_from_complex pole in
    let lti_2nd_xfer = xfer_conjugate2 r arc in
    let hw6 = Array.map (fun w -> Complex.inv (lti_2nd_xfer w)) wn in
    let hmag6 = cpower_db hw6 in

    let g = P.init ~xsize:800. ~ysize:600. Gnuplot.X in
    P.env g (* ~xlog:true *) (-7.) (7.) (-40.) (40.) ;
    P.pen g 1 ;
    P.xy g wn hmag4 ;
    P.pen g 2 ;
    P.xy g wn hmag5 ;
    P.pen g 3 ;
    P.xy g wn hmag6 ;
    P.close g ;
    (* seems same output as matlab freqz *)


