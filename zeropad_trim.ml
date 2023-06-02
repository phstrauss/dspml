(*  one-dimensional discrete-time signal processing, with audio in mind.

    Zero padding and trimming helpers
    License: see accompanying LICENSE.txt or LICENSE.md file.

    © Philippe Strauss, 2009 - 2012, 2020

    References :

        [OppSch]    : Discrete-Time Signal Processing, Oppenheim & Schafer, Prentice Hall, 3rd ed.
        [ProMan]    : Digital Signal Processing, Proakis & Manolakis, Macmillan, 2nd ed.
        [OppWill]   : Signals & Systems, Oppenheim & Willsky, Prentice Hall, Intl 2nd ed.
        [ElFilt]    : Electronic Filter Design Handbook, 4th Ed.
        [IngPro]    : Digital Signal Processing using MATLAB, Ingle & Proakis  *)


open Oneliners
open Vec_simple
open Vec_search
open Minmax
module Log = Logs


let upper_pow2 n =
    let iter = ref 0 in
    let rec div2 n =
        if n > 1 then begin
            incr iter ;
            div2 (n/2)
        end else () in
    div2 n ;
    let exp = !iter + 1 in
    (exp, iof (2. ** (foi exp)))

let zeropad_upper_pow2 samples =
    let len = arlen samples in
    let _, next_pow2 = upper_pow2 len in
    let zeroes = Array.make (next_pow2 - len) 0.0 in
    Array.append samples zeroes

let zeropad_min_pow2 samples exp2 =
    let len1 = arlen samples in
    let lenreq = iof (2. ** (foi exp2)) in
    if lenreq > len1 then begin
        let zeroes = Array.make (lenreq - len1) 0.0 in
        Array.append samples zeroes
    end else samples


let trim_before_fixed compact samples_before =
    let len = arlen compact in
    let idmax, _ = find_min_or_max compact in
    let start =
        if idmax >= samples_before
        then idmax - samples_before
        else 0 in
    Array.sub compact start (len-start)

(* exception Threshold_unfit_for_data *)

(* WW : uninspired naming
   EE : threshold not fit for data *)
let triggers compact th_before th_after =
    let compact_db = rfield_db compact in
    let start, _ = try
        above_from_begin compact_db th_before with AboveFrom_Not_found -> (0, 0.) in
    let stop, _ = try
        above_from_end compact_db th_after with AboveFrom_Not_found -> (arlen compact, 0.) in
    let idxmax, _ = find_min_or_max compact in
    [%log debug "triggers: start=%d, max=%d, stop=%d, len=%d" start idxmax stop (stop-start+1)] ;
    (start, idxmax, stop)

let halfmax compact =
    let (_, minmax) = find_min_or_max compact in
    let half = abs_float minmax /. 2. in
    let idx = ref 0 in
    while (abs_float compact.(!idx)) < half do
        incr idx
    done ;
    [%log debug "halfmax: trigg=%d" !idx] ;
    !idx

let triggers_halfmax compact th_before th_after =
    let start, _, stop = triggers compact th_before th_after in
    let trigg = halfmax compact in
    (start, trigg, stop)

let trim_bothends impulse th_before th_after =
    let start, _, stop = triggers impulse th_before th_after in
    Array.sub impulse start (stop-start+1)
