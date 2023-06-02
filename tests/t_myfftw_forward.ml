open Myfftw
open Oneliners
open Complex
open Vec_sp


(* val fold_left : ('a -> 'b -> 'a) -> 'a -> 'b array -> 'a
Array.fold_left f x a computes f (... (f (f x a.(0)) a.(1)) ...)
where n is the length of the array a.
val fold_right : ('b -> 'a -> 'a) -> 'b array -> 'a -> 'a
Array.fold_right f a x computes f a.(0) (f a.(1) ( ...
where n is the length of the array a. *)

let len1 = 4096
let len1' = foi len1

let len2 = 4095
let len2' = foi len2

let k1 = 200.
let k2 = 231.

let disp_epsilon = 1e-3

let pr = Printf.printf

let print_values carr =
    let len_v = arlen carr in
    for i = 0 to len_v - 1 do
        if Complex.norm carr.(i) > disp_epsilon then
            pr "sp.(%d) = {re=%f; im=%f}, mag=%f\n" i carr.(i).re carr.(i).im (Complex.norm carr.(i))
    done


let () =

    let sig_dc = Array.init len1 (fun i -> 1.0) in
    let e_dc = energy sig_dc in
    let sig_1 = Array.init len1 (fun i -> let i' = foi i in sin (2.*.pi*.i'/.len1')) in
    let e_1 = energy sig_1 in
    let sig_k1 = Array.init len1 (fun i -> let i' = foi i in sin (2.*.pi*.i'*.k1/.len1')) in
    let e_k1 = energy sig_k1 in
    let sig_k12 = Array.init len1 (fun i -> let i' = foi i in
           0.5 *. sin (2.*.pi*.i'*.k1/.len1')
        +. 0.5 *. sin (2.*.pi*.i'*.k2/.len1')) in
    let e_k12 = energy sig_k12 in
    let sig_ny = Array.init len1 (fun i -> let i' = foi i in cos (pi*.i')) in
    let e_ny = energy sig_ny in

    let ftc_f = new fft_complex_t len1 FFT.Forward in

    let sp_dc = Array.make len1 Complex.zero in
    let sp_1 = Array.make len1 Complex.zero in
    let sp_k1 = Array.make len1 Complex.zero in
    let sp_k12 = Array.make len1 Complex.zero in
    let sp_ny = Array.make len1 Complex.zero in

    pr "Testing complex FFT of length %d with real input:\n" len1 ;
    pr "sig_dc : input array all 1.0\n" ;
    ftc_f #fillin_r2c sig_dc ; ftc_f #exec sp_dc ;
    print_values sp_dc ;
    energy2c e_dc sp_dc ;
    pr "sig_1 : sine of rel. freq of 1.0\n" ;
    ftc_f #fillin_r2c sig_1 ; ftc_f #exec sp_1 ;
    print_values sp_1 ;
    energy2c e_1 sp_1 ;
    pr "sig_k1 : sine of rel. freq of %f\n" k1 ;
    ftc_f #fillin_r2c sig_k1 ; ftc_f #exec sp_k1 ;
    print_values sp_k1 ;
    energy2c e_k1 sp_k1 ;
    pr "sig_k12 : sine of rel. freq of %f and %f, HALF AMPLITUDE\n" k1 k2 ;
    ftc_f #fillin_r2c sig_k12 ; ftc_f #exec sp_k12 ;
    print_values sp_k12 ;
    energy2c e_k12 sp_k12 ;
    pr "sig_ny : cosine at nyquist freq\n" ;
    ftc_f #fillin_r2c sig_ny ; ftc_f #exec sp_ny ;
    print_values sp_ny ;
    energy2c e_ny sp_ny ;

    pr "\n" ;

    let ftc_r = new fft_real_t len1 in

    let len_c = hc_len_complex len1 in

    let sp_dc = Array.make len_c Complex.zero in
    let sp_1 = Array.make len_c Complex.zero in
    let sp_k1 = Array.make len_c Complex.zero in
    let sp_k12 = Array.make len_c Complex.zero in
    let sp_ny = Array.make len_c Complex.zero in

    pr "Now testing halfcomplex FFT still of length %d:\n" len1 ;
    pr "sig_dc : input array all 1.0\n" ;
    ftc_r #exec_fwd sig_dc sp_dc ;
    print_values sp_dc ;
    energy2c e_dc sp_dc ;
    pr "sig_1 : sine of rel. freq of 1.0\n" ;
    ftc_r #exec_fwd sig_1 sp_1 ;
    print_values sp_1 ;
    energy2c e_1 sp_1 ;
    pr "sig_k1 : sine of rel. freq of %f\n" k1 ;
    ftc_r #exec_fwd sig_k1 sp_k1 ;
    print_values sp_k1 ;
    energy2c e_k1 sp_k1 ;
    pr "sig_k12 : sine of rel. freq of %f and %f, HALF AMPLITUDE\n" k1 k2 ;
    ftc_r #exec_fwd sig_k12 sp_k12 ;
    print_values sp_k12 ;
    energy2c e_k12 sp_k12 ;
    pr "sig_ny : cosine at nyquist freq\n";
    ftc_r #exec_fwd sig_ny sp_ny ;
    print_values sp_ny ;
    energy2c e_ny sp_ny ;

    (* fft_complex_t of even length are OK rel. normalization *)

    pr "\n" ;

    let sig_dc = Array.init len2 (fun i -> 1.0) in
    let e_dc = energy sig_dc in
    let sig_1 = Array.init len2 (fun i -> let i' = foi i in sin (2.*.pi*.i'/.len2')) in
    let e_1 = energy sig_1 in
    let sig_k1 = Array.init len2 (fun i -> let i' = foi i in sin (2.*.pi*.i'*.k1/.len2')) in
    let e_k1 = energy sig_k1 in
    let sig_k12 = Array.init len2 (fun i -> let i' = foi i in
           0.5 *. sin (2.*.pi*.i'*.k1/.len2')
        +. 0.5 *. sin (2.*.pi*.i'*.k2/.len2')) in
    let e_k12 = energy sig_k12 in
    let sig_ny = Array.init len2 (fun i -> let i' = foi i in cos (pi*.i'*.(len2'-.1.)/.len2')) in
    let e_ny = energy sig_ny in

    let ftc_f = new fft_complex_t len2 FFT.Forward in

    let sp_dc = Array.make len2 Complex.zero in
    let sp_1 = Array.make len2 Complex.zero in
    let sp_k1 = Array.make len2 Complex.zero in
    let sp_k12 = Array.make len2 Complex.zero in
    let sp_ny = Array.make len2 Complex.zero in

    pr "Testing complex FFT of length %d with real input:\n" len2 ;
    pr "sig_dc : input array all 1.0\n" ;
    ftc_f #fillin_r2c sig_dc ; ftc_f #exec sp_dc ;
    print_values sp_dc ;
    energy2c e_dc sp_dc ;
    pr "sig_1 : sine of rel. freq of 1.0\n" ;
    ftc_f #fillin_r2c sig_1 ; ftc_f #exec sp_1 ;
    print_values sp_1 ;
    energy2c e_1 sp_1 ;
    pr "sig_k1 : sine of rel. freq of %f\n" k1 ;
    ftc_f #fillin_r2c sig_k1 ; ftc_f #exec sp_k1 ;
    print_values sp_k1 ;
    energy2c e_k1 sp_k1 ;
    pr "sig_k12 : sine of rel. freq of %f and %f, HALF AMPLITUDE\n" k1 k2 ;
    ftc_f #fillin_r2c sig_k12 ; ftc_f #exec sp_k12 ;
    print_values sp_k12 ;
    energy2c e_k12 sp_k12 ;
    pr "sig_ny : cosine NEAREST of nyquist freq\n" ;
    ftc_f #fillin_r2c sig_ny ; ftc_f #exec sp_ny ;
    print_values sp_ny ;
    energy2c e_ny sp_ny ;

    (* fft_complex_t also OK rel. normalization for odd length *)

    pr "\n" ;

    let ftc_r = new fft_real_t len2 in

    let len_c = hc_len_complex len2 in

    let sp_dc = Array.make len_c Complex.zero in
    let sp_1 = Array.make len_c Complex.zero in
    let sp_k1 = Array.make len_c Complex.zero in
    let sp_k12 = Array.make len_c Complex.zero in
    let sp_ny = Array.make len_c Complex.zero in

    pr "Now testing halfcomplex FFT still of length %d:\n" len2 ;
    pr "sig_dc : input array all 1.0\n" ;
    ftc_r #exec_fwd sig_dc sp_dc ;
    print_values sp_dc ;
    energy2c e_dc sp_dc ;
    pr "sig_1 : sine of rel. freq of 1.0\n" ;
    ftc_r #exec_fwd sig_1 sp_1 ;
    print_values sp_1 ;
    energy2c e_1 sp_1 ;
    pr "sig_k1 : sine of rel. freq of %f\n" k1 ;
    ftc_r #exec_fwd sig_k1 sp_k1 ;
    print_values sp_k1 ;
    energy2c e_k1 sp_k1 ;
    pr "sig_k12 : sine of rel. freq of %f and %f, HALF AMPLITUDE\n" k1 k2 ;
    ftc_r #exec_fwd sig_k12 sp_k12 ;
    print_values sp_k12 ;
    energy2c e_k12 sp_k12 ;
    pr "sig_ny : cosine NEAREST of nyquist freq\n" ;
    ftc_r #exec_fwd sig_ny sp_ny ;
    print_values sp_ny ;
    energy2c e_ny sp_ny ;
