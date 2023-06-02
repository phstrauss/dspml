(*  one-dimensional discrete-time signal processing, with audio in mind.

    Windowing functions around cosinus fade-in-out slopes

    © Philippe Strauss, 2009 - 2012
    License: see accompanying LICENSE.txt or LICENSE.md file.

    References :

        [OppSch]    : Discrete-Time Signal Processing, Oppenheim & Schafer, Prentice Hall, 3rd ed.
        [ProMan]    : Digital Signal Processing, Proakis & Manolakis, Macmillan, 2nd ed.
        [OppWill]   : Signals & Systems, Oppenheim & Willsky, Prentice Hall, Intl 2nd ed.
        [ElFilt]    : Electronic Filter Design Handbook, 4th Ed.
        [IngPro]    : Digital Signal Processing using MATLAB, Ingle & Proakis  *)


open Oneliners


(*  Raised cosinus based windows and cross-fades  *)

type fade_direction = In | Out

(* raised cosine sample coeff, not array *)
let cos_fade direction l n =
    match direction with
    | Out -> 0.5 *. (1. +. cos (foi n *. pi /. foi l))
    | In ->  0.5 *. (1. -. cos (foi n *. pi /. foi l))

(* 2 arrays, cross raised cosine fade *)
let x_cos_fade l =
    let out = Array.make_matrix 2 l 0. in
    for i = 0 to (l-1) do
        out.(0).(i) <- cos_fade In l i ;
        out.(1).(i) <- cos_fade Out l i
    done ;
    out

(* 2 arrays, 2 fade in/out pairs *)
let x_cos_fade_window l l0 l1 =
    (* assert(l0 < l1 < l) ; *)
    let out = Array.make_matrix 2 l 0. in
    let ltr = l1 - l0 in
    for i = 0 to (l-1) do
        let point = if i < l0 then (0., 1.)
              else if i >= l0 && i < l1 then (cos_fade In ltr (i-l0), cos_fade Out ltr (i-l0))
              else if i >= l1 then (1., 0.)
              else (0., 0.) in
        out.(0).(i) <- fst point ;
        out.(1).(i) <- snd point
    done ;
    out

(* tukey window *)
let tukey l in0 in1 out1 out0 =
    let out = Array.make l 0. in
    for i = in0 to out0 do
        out.(i) <-  if      i < in1                 then cos_fade In  (in1-in0) (i-in0)
                    else if i >= in1 && i < out1    then 1.
                    else if i >= out1               then cos_fade Out (out0-out1) (i-out1)
                    else 0.
    done ;
    out

(* same as above, arg based on time, not samples *)
let tukey_time l sr tin0 tin1 tout1 tout0 =
    let srf  = foi sr in
    let in0  = iof (tin0  *. srf)
    and in1  = iof (tin1  *. srf)
    and out1 = iof (tout1 *. srf)
    and out0 = iof (tout0 *. srf) in
    tukey l in0 in1 out1 out0

(* Hann coefficient *)
let hann_coeff l n =
    0.5 *. (1. -. cos (foi n *. pi *. 2. /. foi (l-1)))

(* Hann window *)
let hann_window l =
    let finit = hann_coeff l in
    Array.init l finit

(* Hamming coefficient *)
let hamming_coeff l n =
    0.54 -. 0.46 *. cos (foi n *. pi *. 2. /. foi (l-1))

(* Hamming window *)
let hamming_window l =
    let finit = hamming_coeff l in
    Array.init l finit
