(*  one-dimensional discrete-time signal processing, with audio in mind.

    Tiny object wrapper around fftw for most used cases
    ocaml arrays needs to be allocated once beforehand for efficiency
    fftw bigarrays and plans have the lifespan of class instantiation

    © Philippe Strauss, 2012, 2015
    License: see accompanying LICENSE.txt or LICENSE.md file.

    References :

        [OppSch]    : Discrete-Time Signal Processing, Oppenheim & Schafer, Prentice Hall, 3rd ed.
        [ProMan]    : Digital Signal Processing, Proakis & Manolakis, Macmillan, 2nd ed.
        [OppWill]   : Signals & Systems, Oppenheim & Willsky, Prentice Hall, Intl 2nd ed.
        [ElFilt]    : Electronic Filter Design Handbook, 4th Ed.
        [IngPro]    : Digital Signal Processing using MATLAB, Ingle & Proakis  *)


open Oneliners
open Bigarray
open Complex
module FFT = Fftw3.D


(* type 'a plan (** FFTW plan. *)
   type c2c     (** [c2c plan] usual discrete Fourier transform,
                    from complex to complex *)
   type r2c     (** [r2c plan] real to complex transform *)
   type c2r     (** [c2r plan] complex to real transform *)
   type r2r     (** [r2r plan] real to real transform *)

   type r2r_kind =
    | R2HC (** real to halfcomplex *)
    | HC2R (** halfcomplex to real *) *)

class fft_complex_t dim direction =
    let x = FFT.Array1.create FFT.complex c_layout dim in
    let y = FFT.Array1.create FFT.complex c_layout dim in (* Fourier vars *)
    let plan = FFT.Array1.dft ~meas:FFT.Estimate direction x y in
    object
        method fillin_r2c din =
            Array.iteri (fun i v -> x.{i} <- { re = v ; im = 0.0 }) din ;
        method fillin_c2c din =
            Array.iteri (fun i v -> x.{i} <- v) din ;
        method exec ?(scale = 1. /. Stdlib.sqrt (foi dim)) dout =
            (* let scale = 1. /. Stdlib.sqrt (foi dim) in *) (* 2015-11-16 after energy bench *)
            FFT.exec plan ;
            Array.iteri (fun i _ -> let yi = y.{i} in
                dout.(i) <- {re = yi.re *. scale; im = yi.im *. scale}) dout
    end

(* fftw3.pdf :

   § 2.5.1 The Halfcomplex-format DFT
   An r2r kind of FFTW_R2HC (r2hc) corresponds to an r2c DFT (see Section 2.3 [One-Dimensional DFTs of
   Real Data], page 6) but with “halfcomplex” format output, and may sometimes be faster and/or more
   convenient than the latter. The inverse hc2r transform is of kind FFTW_HC2R. This consists of the
   non-redundant half of the complex output for a 1d real-input DFT of size n, stored as a sequence
   of n real numbers (double) in the format:

   r0, r1, r2, ..., rn/2,  i(n+1)/2−1, ..., i2, i1

    Here, rk is the real part of the kth output, and ik is the imaginary part.
    (Division by 2 is rounded down.) For a halfcomplex array hc[n], the kth component thus has its
    real part in hc[k] and its imaginary part in hc[n-k], with the exception of k == 0 or n/2
    (the latter only if n is even)—in these two cases, the imaginary part is zero due to symmetries
    of the real-input DFT, and is not stored. Thus, the r2hc transform of n real values is a
    halfcomplex array of length n, and vice versa for hc2r.
*)

(* len of complex array as of function of HC array length *)
let hc_len_complex lenhc = half_up (lenhc + 1)

(* half complex -> complex coefficients layout *)
let hc2c scale halfcomplex complex =
    let len_hc = balen halfcomplex in
    (* let scale = 1. /. Stdlib.sqrt (foi len_hc) in *)
    let len_c = arlen complex in
    assert (len_c = (hc_len_complex len_hc)) ;
    let even = ref (if (len_hc mod 2 = 0) then true else false) in
    let stop = if !even then len_hc/2-1 else (len_hc-1) / 2 in  
    for k = 1 to stop do
        complex.(k) <- { re = halfcomplex.{k} *. scale *. Stdlib.sqrt 2. ;
                         im = halfcomplex.{len_hc-k} *. scale *. Stdlib.sqrt 2. } 
    done ;
    if !even then (
        complex.(len_c-1) <- {re = halfcomplex.{len_hc/2} *. scale ; im = 0.}
    ) ;
    (* hermitian sym. and periodicity prop. *)
    complex.(0) <- { re = halfcomplex.{0} *. scale ; im = 0. }

let c2hc scale complex halfcomplex =
    let len_c = arlen complex in
    let len_hc = balen halfcomplex in
    (* let scale = 1. /. Stdlib.sqrt (foi len_hc) in *)
    assert (len_c = (hc_len_complex len_hc)) ;
    let even = ref (if (len_hc mod 2 = 0) then true else false) in
    let stop = if !even then len_c-2 else len_c-1 in
    for k = 1 to stop do (* FIX : think len_hc 256, len_c = 129 => overlap! seems ok for odd HC length *)
        halfcomplex.{k} <- complex.(k).re *. scale /. Stdlib.sqrt 2.;
        halfcomplex.{len_hc-k} <- complex.(k).im *. scale /. Stdlib.sqrt 2. ;
    done ;
    if !even then (
        halfcomplex.{len_hc/2} <- complex.(len_c-1).re *. scale
    ) ;
    halfcomplex.{0} <- complex.(0).re *. scale


class fft_real_t dim =
    let x = FFT.Array1.create FFT.float c_layout dim
    and y = FFT.Array1.create FFT.float c_layout dim in
    let plan_fwd = FFT.Array1.r2r FFT.R2HC ~meas:FFT.Estimate x y
    and plan_back = FFT.Array1.r2r FFT.HC2R ~meas:FFT.Estimate x y in
    object         
        method exec_fwd ?(scale = 1. /. Stdlib.sqrt (foi dim)) din dout =
            Array.iteri (fun i v -> x.{i} <- v) din ;
            FFT.exec plan_fwd ;
            hc2c scale y dout
        method exec_back ?(scale = 1. /. Stdlib.sqrt (foi dim)) din dout =
            c2hc scale din x ;
            FFT.exec plan_back ;
            Array.iteri (fun i _ -> dout.(i) <- y.{i}) dout
    end
