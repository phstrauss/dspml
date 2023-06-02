(*  one-dimensional discrete-time signal processing, with audio in mind.

    numerical integration, differentiation on numerical data/signals.

    © Philippe Strauss, 2021
    License: see accompanying LICENSE.txt or LICENSE.md file.

    References :

        [OppSch]    : Discrete-Time Signal Processing, Oppenheim & Schafer, Prentice Hall, 3rd ed.
        [ProMan]    : Digital Signal Processing, Proakis & Manolakis, Macmillan, 2nd ed.
        [OppWill]   : Signals & Systems, Oppenheim & Willsky, Prentice Hall, Intl 2nd ed.
        [ElFilt]    : Electronic Filter Design Handbook, 4th Ed.
        [IngPro]    : Digital Signal Processing using MATLAB, Ingle & Proakis  *)


open Oneliners


(* trapezoidal integration *)
let integrate_trapezoidal ?(scale=1.) fs signal =
    let len = arlen signal in
    let accum = ref 0. in
    let h = 1. /. (foi fs) in
    for i = 0 to len-2 do
        let height = (signal.(i) +. signal.(i+1)) /. 2. in
        accum := !accum +. height *. h *. scale ;
    done ;
    !accum

let deriv_naive signal =
    let len = arlen signal in
    let diff = Array.make (len-1) 0. in
    for i = 0 to len-2 do
        diff.(i) <- signal.(i+1) -. signal.(i) ;
    done ;
    diff

type coeffs_t = Coeffs of float array | File of string

let filter_fir coeffs signal =
    let c = match coeffs with
        | Coeffs x -> x
        | File path -> Sp_misc.read_filter_coeffs path in
    let conv = new Spectral.convolve_fft_t (* ~scale_fwd:1. *) (* ~scale_back:1. *) c (arlen signal) in (* WW: still needs scaling done right ! *)
        conv #exec signal

let deriv_fir = filter_fir
let hilbert_xform_fir = filter_fir