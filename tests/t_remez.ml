open Vec_sp
open Common
module P = Gnuplot.Array


let () =

    let r12 = 2. ** (1./.12.) in
    let r24 = 2. ** (1./.24.) in

    let notes_l = List.rev (Smooth.chroma_notes 440.5 11025. 12.) in
    let notes = Array.of_list (BatList.take 14 notes_l) in
    Printf.printf "chromatic notes freq.: %s\n\n" (sprint_farray notes) ;

    (* test specification in frequencies bands, conventional API *)

    let tr = 0.25 /. r12 in

    let edges = [| 0.; tr; 0.25; 0.5 |] in
    let desired = [|1.; 0.|] in
    let weights = [| 1. ; 10. |] in
    let iter, h = Remez.remez_from_bands 371 2 edges desired weights Remez.Standard 0.0001 in
    Printf.printf "%d iterations, filter coefficients: %s\n\n" iter (sprint_farray h) ;

    let hp = Zeropad_trim.zeropad_min_pow2 h 12 in
    let len2 = arlen hp in
    let spectrum = Array.make (half len2) (Complex.zero) in
    let fft = new Myfftw.fft_real_t len2 in
    fft #exec_fwd hp spectrum ;
    let fsp = fgrid 1 len2 in

    let g = P.init ~xsize:800. ~ysize:400. Gnuplot.X in
    P.env g (* ~xlog:true *) (0.) (0.5) (-140.) (10.) ;

    P.pen g 1 ;
    P.xy g fsp (cpower_db spectrum) ;

    (* new, flexible remez_from_grid API *)

    let ncoeffs = 401 in
    let ncoeffs_r = (ncoeffs - 1) / 2 in
    let fdensity = 8 in
    let nf = ncoeffs_r * fdensity in

    let f = Array.init nf (fun i -> foi i *. 0.5 /. foi nf) in

    let desired = Array.make nf 0. in
    Smooth.one_bell_log 0.25 (* la3 is 0.01 for 44.1 k sampling :*) r12 f desired ;
    let weight = Array.map (fun x -> 10. -. 9. *. x) desired in

    let iter, h = Remez.remez_from_grid ncoeffs f desired weight nf Remez.Standard 0.0001 in
    Printf.printf "%d iterations, filter coefficients: %s\n\n" iter (sprint_farray h) ;

    let hp = Zeropad_trim.zeropad_min_pow2 h 13 in
    let len2 = arlen hp in
    let spectrum = Array.make (half len2) (Complex.zero) in
    let fft = new Myfftw.fft_real_t len2 in
    fft #exec_fwd hp spectrum ;
    let fsp = fgrid 1 len2 in

    P.pen g 3 ;
    P.xy g fsp (cpower_db spectrum) ;

    P.close g