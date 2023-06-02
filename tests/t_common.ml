open Oneliners
open Vec_sp
module P = Gnuplot.Array
open Myfftw


type domain_t = Time | Frequency


let plot tsamples fs domain =
    let len = Array.length tsamples in
    let ffs = foi fs in
    let g = P.init ~xsize:800. ~ysize:600. Gnuplot.X in
    match domain with
    | Time -> begin
        let xt = tgrid fs len in
        P.env g ~xlog:false (0.) (30e-3) (-1.) (1.) ;
        P.pen g 1 ;
        P.xy g xt tsamples ;
        P.close g
    end
    | Frequency -> begin
        let ftr = new fft_real_t len in
        let len_c = hc_len_complex len in
        let xf = fgrid (fs/2) len_c in
        let sp1_tr = Array.make len_c Complex.zero in
        ftr #exec_fwd tsamples sp1_tr ~scale:(1./.(sqrt 2.));
        let sp1_norm = cmagnitude sp1_tr in
        let sp1_db = log_db20 sp1_norm in
        P.env g ~xlog:true (10.) (ffs/.2.) (-80.) (10.) ;
        P.pen g 3 ;
        P.xy g (nodc xf) (nodc sp1_db) ;
        P.close g
    end
