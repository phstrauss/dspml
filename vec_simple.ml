(*  one-dimensional discrete-time signal processing, with audio in mind.

    Basics helpers when fiddling with float or complex 1D vectors

    © Philippe Strauss, 2009 - 2012
    License: see accompanying LICENSE.txt or LICENSE.md file.

    References :

        [OppSch]    : Discrete-Time Signal Processing, Oppenheim & Schafer, Prentice Hall, 3rd ed.
        [ProMan]    : Digital Signal Processing, Proakis & Manolakis, Macmillan, 2nd ed.
        [OppWill]   : Signals & Systems, Oppenheim & Willsky, Prentice Hall, Intl 2nd ed.
        [ElFilt]    : Electronic Filter Design Handbook, 4th Ed.
        [IngPro]    : Digital Signal Processing using MATLAB, Ingle & Proakis  *)


open Oneliners
open Complex
open Vec_search
module Log = Logs


let creal cvec =
    Array.map (fun x -> x.re) cvec

let cimag cvec =
    Array.map (fun x -> x.im) cvec

let cmagnitude cvec =
    Array.map (fun x -> Complex.norm x) cvec

(* square each elements of a vector *)
let rpower y =
    Array.map (fun x -> x ** 2.) y

let cpower cvec =
    Array.map (fun x -> Complex.norm2 x) cvec

let cargument cvec =
    Array.map (fun x -> Complex.arg x) cvec

let scale k vec =
    Array.map (fun x -> x *. k) vec

(* common vector operations for signals - the 10 following mostly for interfacing graphical representation *)

let eps_db_pow = 1e-14
let eps_db_field = 1e-7

(* decibel of a power kind of value, as opposed to a field value *)
let log_db10 y =
    Array.map (fun x -> if x > eps_db_pow then 10. *. log10 x else -(140.0)) y

(* decibel of a field value *)
let log_db20 y =
    Array.map (fun x -> if x > eps_db_field then 20. *. log10 x else -(140.0)) y

let cpower_db y =
    Array.map (fun x ->
        let powval = Complex.norm2 x in
        if powval > eps_db_pow then 10. *. log10 powval else -(140.)
    ) y

let rfield_db y =
    Array.map (fun x ->
        let powval = x ** 2. in
        if powval > eps_db_pow then 10. *. (log10 powval) else -(140.)
    ) y

(*  originally in myfftw.ml : ENERGY PRESERVING NORMALIZATION DONE I HOPE RIGHT THIS TIME
    2015-11-16 (1/√N - physicists)  *)

let cenergy x =
    Array.fold_left (fun accum y -> Complex.norm2 y +. accum) 0. x

let energy x =
    Array.fold_left (fun accum y -> y**2. +. accum) 0. x

let cenergy2 e1 x =
    let e2 = cenergy x in
    pr "cenergy2 : e1=%f, e2=%f, e2/e1=%f, e1/e2=%f\n%!" e1 e2 (e2/.e1) (e1/.e2) 

let energy2 e1 x =
    let e2 = energy x in
    pr "energy2 : e1=%f, e2=%f, e2/e1=%f, e1/e2=%f\n%!" e1 e2 (e2/.e1) (e1/.e2)

(* remove element 0, the continuous current component - for log transformation on the x-axis *)
let nodc y =
    let len = arlen y in
    Array.sub y 1 (len-1)

(* generate a time x-axis *)
let tgrid sfreq len =
    Array.init len (fun i -> foi i /. foi sfreq)

(* generate a freq x-axis *)
let fgrid sfreq len =
    Array.init len (fun i -> foi i *. foi sfreq /. foi len)

let fgrid2 sfreq len =
    Array.init len (fun i -> foi i *. foi sfreq /. (foi len *. 2.))

let fgrid_half sfreq len =
    Array.sub (fgrid sfreq len) 0 (half_down len)

let fgrid_log sfreq len =
    Array.map (fun x -> log10_ex x) (nodc (fgrid sfreq len))

let fgrid_half_log sfreq len =
    Array.map (fun x -> log10_ex x) (nodc (fgrid_half sfreq len))

let idx2time idx sfreq =
    foi idx /. foi sfreq


let energy_time_curve compact =
    let (tmax, minmax) = Minmax.find_min_or_max compact in
    let amp_scaled = Array.map (fun x -> x /. minmax) compact in
    let power = Array.map (fun x -> 10. *. log10 (x ** 2.)) amp_scaled in
    (power, tmax, minmax)


(* mostly for putting the positive time part of a hilbert filter at the
   beginning of an array, the negative time part at the end *)
let split_nyquist hn out_len =
    let out = Array.make out_len 0. in
    let len_h = arlen hn in let len_h2 = half_down len_h in
    assert (out_len >= len_h) ;
    (* split the hilbert coeffs in the upper half, aliased negative frequencies *)
    for i = 0 to (len_h2-1) do
        out.(i) <- hn.(i+len_h2) ;
        out.(out_len-i-1) <- hn.(len_h2-i-1)
    done ;
    out

(* used with howmuch = len/2, freq. above nyquist are translated below and vice versa *)
let splice_nyquist arr howmuch =
    let len = arlen arr in
    Array.append (Array.sub arr (len-howmuch) howmuch) (Array.sub arr 0 howmuch)

(* do a upper even symetry of an array - output double the size of it's input *)
let even_from_single half_spectra scale =
    let len = arlen half_spectra in
    let adj = if scale then 0.5 else 1. in
    let mirror i = if i < len then adj *. half_spectra.(i) else adj *. half_spectra.(2*len - i - 1) in
    Array.init (2 * len) mirror

let odd_from_single half_spectra scale =
    let len = arlen half_spectra in
    let adj = if scale then 0.5 else 1. in
    let mirror i = if i < len then adj *. half_spectra.(i) else (-.adj) *. half_spectra.(2*len - i - 1) in
    Array.init (2 * len) mirror

let single_from_even full_spectra =
    let len = arlen full_spectra in
    let ret = if (len mod 2 = 0) then
        Array.init (len/2) (fun i -> full_spectra.(i) +. full_spectra.(len-1-i))
    else
        Array.init (half_down len + 1) (fun i -> full_spectra.(i) +. full_spectra.(len-1-i)) in
        (* CHECK: is the middle elt in odd len case to be divided by two ? *)
    ret


let vec_field_avg mag fgrid fl fh =
    let idx_l = index_of_freq fgrid fl
    and idx_h = index_of_freq fgrid fh
    and sum = ref 0. in
    for i = idx_l to idx_h do
        (* we want avg power of spectrum normalized to 1.0 *)
        sum := !sum +. mag.(i) ;
    done ;
    (!sum /. foi (idx_h - idx_l))

(* used for normalize just below - BEWARE do not pass a squared or power value *)
let vec_rms_avg mag fgrid fl fh =
    let idx_l = index_of_freq fgrid fl
    and idx_h = index_of_freq fgrid fh
    and sum = ref 0. in
    for i = idx_l to idx_h do
        (* we want avg power of spectrum normalized to 1.0 *)
        sum := !sum +. mag.(i) ** 2. ;
    done ;
    Stdlib.sqrt (!sum /. foi (idx_h - idx_l))

let spectrum_normalize mag fgrid fl fh =
    let avg = vec_rms_avg mag fgrid fl fh
    and len_m = arlen mag in
    [%log debug "normalize_spectrum: sqrt of power avg: %f" avg] ;
    for i = 0 to (len_m-1) do
        mag.(i) <- mag.(i) /. avg ;
    done ;
    mag
