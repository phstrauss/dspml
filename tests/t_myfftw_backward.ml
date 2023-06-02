open Myfftw
open Oneliners
open Complex
open Vec_sp
module P = Gnuplot.Array


let len1 = 256
let len1' = foi len1

let len2 = 255
let len2' = foi len2

let pr = Printf.printf

let print_cvalues carr =
    let len_v = arlen carr in
    for i = 0 to len_v - 1 do
        pr "sp.(%d) = {re=%f; im=%f}, mag=%f\n" i carr.(i).re carr.(i).im (Complex.norm carr.(i))
    done

let plot pen signal =
    let len = Array.length signal in
    let len' = foi len in
    let g = P.init ~xsize:800. ~ysize:600. Gnuplot.X in
    P.env g ~xlog:false 0. len' (-0.1) 0.1 ;
    P.pen g pen ;
    let x = Array.init len (fun i -> foi i) in
    P.xy g x signal ;
    P.close g


let () =

    (* ***** Complex-Complex EVEN, freq = 0, DC ***** *)

    pr "Testing complex IFFT of EVEN length %d.\nPress a key to begin...\n%!" len1 ;
    let _ = input_char stdin in

    let ftc_b = new fft_complex_t len1 FFT.Backward in

    let sp_dc = Array.make len1 Complex.zero in
    sp_dc.(0) <- {re = 1.0; im = 0.0} ;
    let sig_dc = Array.make len1 Complex.zero in
    ftc_b #fillin_c2c sp_dc ;
    ftc_b #exec sig_dc ;

    energy2c (energyc sp_dc) sig_dc ;
    plot 1 (Array.init len1 (fun i -> sig_dc.(i).re)) ;

    (* ***** Complex-Complex even, freq = 1 ***** *)

    let sp_1 = Array.make len1 Complex.zero in
    sp_1.(1) <- {re = 0.0; im = -0.5} ;
    sp_1.(len1-1) <- {re = 0.0; im = 0.5} ;
    let sig_1 = Array.make len1 Complex.zero in
    ftc_b #fillin_c2c sp_1 ;
    ftc_b #exec sig_1 ;

    energy2c (energyc sp_1) sig_1 ;
    plot 2 (Array.init len1 (fun i -> sig_1.(i).re)) ;

    (* ***** Complex-Complex even, freq = 5 ***** *)

    let sp_5 = Array.make len1 Complex.zero in
    sp_5.(5) <- {re = 0.0; im = -0.5} ;
    sp_5.(len1-5) <- {re = 0.0; im = 0.5} ;
    let sig_5 = Array.make len1 Complex.zero in
    ftc_b #fillin_c2c sp_5 ;
    ftc_b #exec sig_5 ;

    energy2c (energyc sp_5) sig_5 ;
    plot 3 (Array.init len1 (fun i -> sig_5.(i).re)) ;

    (* ***** Complex-Complex even, freq = 5 + 10 ***** *)

    let sp_5_10 = Array.make len1 Complex.zero in
    sp_5_10.(5) <- {re = 0.0; im = -0.25} ;
    sp_5_10.(len1-5) <- {re = 0.0; im = 0.25} ;
    sp_5_10.(10) <- {re = 0.0; im = -0.25} ;
    sp_5_10.(len1-10) <- {re = 0.0; im = 0.25} ;
    let sig_5_10 = Array.make len1 Complex.zero in
    ftc_b #fillin_c2c sp_5_10 ;
    ftc_b #exec sig_5_10 ;

    energy2c (energyc sp_5_10) sig_5_10 ;
    plot 4 (Array.init len1 (fun i -> sig_5_10.(i).re)) ;

    (* ***** Complex-Complex even, freq = nyquist ***** *)

    let sp_ny = Array.make len1 Complex.zero in
    sp_ny.(len1/2) <- {re = 1.0; im = 0.} ;
    let sig_ny = Array.make len1 Complex.zero in
    ftc_b #fillin_c2c sp_ny ;
    ftc_b #exec sig_ny ;

    energy2c (energyc sp_ny) sig_ny ;
    plot 5 (Array.init len1 (fun i -> sig_ny.(i).re)) ;

    (* ***** Complex-Complex ODD, freq = 0, DC ***** *)

    pr "Testing complex IFFT of ODD length %d.\nPress a key to begin...\n%!" len2 ;
    let _ = input_char stdin in

    let ftc_b = new fft_complex_t len2 FFT.Backward in

    let sp_dc = Array.make len2 Complex.zero in
    sp_dc.(0) <- {re = 1.0; im = 0.0} ;
    let sig_dc = Array.make len2 Complex.zero in
    ftc_b #fillin_c2c sp_dc ;
    ftc_b #exec sig_dc ;

    energy2c (energyc sp_dc) sig_dc ;
    plot 1 (Array.init len2 (fun i -> sig_dc.(i).re)) ;

    (* ***** Complex-Complex odd, freq = 1 ***** *)

    let sp_1 = Array.make len2 Complex.zero in
    sp_1.(1) <- {re = 0.0; im = -0.5} ;
    sp_1.(len2-1) <- {re = 0.0; im = 0.5} ;
    let sig_1 = Array.make len2 Complex.zero in
    ftc_b #fillin_c2c sp_1 ;
    ftc_b #exec sig_1 ;

    energy2c (energyc sp_1) sig_1 ;
    plot 2 (Array.init len2 (fun i -> sig_1.(i).re)) ;

    (* ***** Complex-Complex odd, freq = 5 ***** *)

    let sp_5 = Array.make len2 Complex.zero in
    sp_5.(5) <- {re = 0.0; im = -0.5} ;
    sp_5.(len2-5) <- {re = 0.0; im = 0.5} ;
    let sig_5 = Array.make len2 Complex.zero in
    ftc_b #fillin_c2c sp_5 ;
    ftc_b #exec sig_5 ;

    energy2c (energyc sp_5) sig_5 ;
    plot 3 (Array.init len2 (fun i -> sig_5.(i).re)) ;

    (* ***** Complex-Complex odd, freq = 5 + 10 ***** *)

    let sp_5_10 = Array.make len2 Complex.zero in
    sp_5_10.(5) <- {re = 0.0; im = -0.25} ;
    sp_5_10.(len2-5) <- {re = 0.0; im = 0.25} ;
    sp_5_10.(10) <- {re = 0.0; im = -0.25} ;
    sp_5_10.(len2-10) <- {re = 0.0; im = 0.25} ;
    let sig_5_10 = Array.make len2 Complex.zero in
    ftc_b #fillin_c2c sp_5_10 ;
    ftc_b #exec sig_5_10 ;

    energy2c (energyc sp_5_10) sig_5_10 ;
    plot 4 (Array.init len2 (fun i -> sig_5_10.(i).re)) ;

    (* ***** Complex-Complex ODD, freq = NEAREST OF nyquist ***** *)

    let sp_ny = Array.make len2 Complex.zero in
(*  sp_ny.(half len2) <- {re = 0.5; im = 0.} ;
    sp_ny.((half len2)+1) <- {re = 0.5; im = 0.} ; *)
    sp_ny.(half_down len2) <- {re = 0.5; im = 0.} ;
    sp_ny.((half_down len2)+1) <- {re = 0.5; im = 0.} ;
    let sig_ny = Array.make len2 Complex.zero in
    ftc_b #fillin_c2c sp_ny ;
    ftc_b #exec sig_ny ;

    energy2c (energyc sp_ny) sig_ny ;
    plot 5 (Array.init len2 (fun i -> sig_ny.(i).re)) ;

    (* ***** Half-complex EVEN, freq = 0, DC ***** *)

    pr "Testing half-complex IFFT of EVEN length %d.\nPress a key to begin...\n%!" len1 ;
    let _ = input_char stdin in

    let ftr_b = new fft_real_t len1 in
    let hlen1 = hc_len_complex len1 in

    let sp_dc = Array.make hlen1 Complex.zero in
    sp_dc.(0) <- {re = 1.; im = 0.} ;
    let sig_dc = Array.make len1 0. in
    ftr_b #exec_back sp_dc sig_dc;

    energy2 (energyc sp_dc) sig_dc ;
    plot 1 sig_dc ;

    (* ***** Half-complex even, freq = 1 ***** *)

    let sp_1 = Array.make hlen1 Complex.zero in
    sp_1.(1) <- {re = 0.0; im = -1.0} ;
    let sig_1 = Array.make len1 0. in
    ftr_b #exec_back sp_1 sig_1;

    energy2 (energyc sp_1) sig_1 ;
    plot 2 sig_1 ;

    (* ***** Half-complex even, freq = 5 ***** *)

    let sp_5 = Array.make hlen1 Complex.zero in
    sp_5.(5) <- {re = 0.; im = -1.} ;
    let sig_5 = Array.make len1 0. in
    ftr_b #exec_back sp_5 sig_5;

    energy2 (energyc sp_5) sig_5 ;
    plot 3 sig_5 ;

    (* ***** Half-complex even, freq = 5 + 10 ***** *)

    let sp_5_10 = Array.make hlen1 Complex.zero in
    sp_5_10.(5) <- {re = 0.0; im = -0.5} ;
    sp_5_10.(10) <- {re = 0.0; im = -0.5} ;
    let sig_5_10 = Array.make len1 0. in
    ftr_b #exec_back sp_5_10 sig_5_10;

    energy2 (energyc sp_5_10) sig_5_10 ;
    plot 4 sig_5_10 ;

    (* ***** Half-complex even, freq = nyquist ***** *)

    let sp_ny = Array.make (hlen1) Complex.zero in
    sp_ny.(hlen1-1) <- {re = 1.0; im = 0.} ;
    let sig_ny = Array.make len1 0. in
    ftr_b #exec_back sp_ny sig_ny;

    energy2 (energyc sp_ny) sig_ny ;
    plot 5 sig_ny ;

    (* ***** Half-complex ODD, freq = 0, DC ***** *)

    pr "Testing half-complex IFFT of ODD length %d.\nPress a key to begin...\n%!" len2 ;
    let _ = input_char stdin in

    let ftr_b = new fft_real_t len2 in
    let hlen2 = hc_len_complex len2 in

    let sp_dc = Array.make hlen2 Complex.zero in
    sp_dc.(0) <- {re = 1.; im = 0.} ;
    let sig_dc = Array.make len2 0. in
    ftr_b #exec_back sp_dc sig_dc;

    energy2 (energyc sp_dc) sig_dc ;
    plot 1 sig_dc ;

    (* ***** Half-complex odd, freq = 1 ***** *)

    let sp_1 = Array.make hlen2 Complex.zero in
    sp_1.(1) <- {re = 0.0; im = -1.0} ;
    let sig_1 = Array.make len2 0. in
    ftr_b #exec_back sp_1 sig_1;

    energy2 (energyc sp_1) sig_1 ;
    plot 2 sig_1 ;

    (* ***** Half-complex odd, freq = 5 ***** *)

    let sp_5 = Array.make hlen2 Complex.zero in
    sp_5.(5) <- {re = 0.; im = -1.} ;
    let sig_5 = Array.make len2 0. in
    ftr_b #exec_back sp_5 sig_5;

    energy2 (energyc sp_5) sig_5 ;
    plot 3 sig_5 ;

    (* ***** Half-complex odd, freq = 5 + 10 ***** *)

    let sp_5_10 = Array.make hlen2 Complex.zero in
    sp_5_10.(5) <- {re = 0.0; im = -0.5} ;
    sp_5_10.(10) <- {re = 0.0; im = -0.5} ;
    let sig_5_10 = Array.make len2 0. in
    ftr_b #exec_back sp_5_10 sig_5_10;

    energy2 (energyc sp_5_10) sig_5_10 ;
    plot 4 sig_5_10 ;

    (* ***** Half-complex ODD, NEAREST freq of nyquist ***** *)

    let sp_ny = Array.make (hlen2) Complex.zero in
    sp_ny.(hlen2-1) <- {re = 1.0; im = 0.} ;
    let sig_ny = Array.make len2 0. in
    ftr_b #exec_back sp_ny sig_ny;

    energy2 (energyc sp_ny) sig_ny ;
    plot 5 sig_ny ;

