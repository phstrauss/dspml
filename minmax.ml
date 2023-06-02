(*  one-dimensional discrete-time signal processing, with audio in mind.

    looks for extrema(s) on numerical data.

    © Philippe Strauss, 2009 - 2012, 2015
    License: see accompanying LICENSE.txt or LICENSE.md file.

    References :

        [OppSch]    : Discrete-Time Signal Processing, Oppenheim & Schafer, Prentice Hall, 3rd ed.
        [ProMan]    : Digital Signal Processing, Proakis & Manolakis, Macmillan, 2nd ed.
        [OppWill]   : Signals & Systems, Oppenheim & Willsky, Prentice Hall, Intl 2nd ed.
        [ElFilt]    : Electronic Filter Design Handbook, 4th Ed.
        [IngPro]    : Digital Signal Processing using MATLAB, Ingle & Proakis  *)


open Oneliners


(* find local minimas and maximas on a relatively smooth curve, gives them in
   x ascending order in a list. CAUTION: min and max are interleaved in a single
   output list.
   N.B.: four next funcs, four perf. bug for the price of one *)
(* let local_min_max1 arr =
    let diff data =
        let len = arlen data in
        let d' = Array.make len 0. in
        for i = 0 to (len-2) do
            d'.(i+1) <- data.(i+1) -. data.(i)
        done ; d' in
    let arr' = diff arr in
    let len = arlen arr in
    let minmax = ref [] in
    for i = 0 to (len-2) do
        if (arr'.(i) *. arr'.(i+1)) < 0. then begin
            [%log debug "local_min_max1: local minmax found at index %d; value %f" i arr.(i)];
            minmax :=  ( (i , arr.(i)) :: !minmax )
        end ;
    done ;
    List.rev !minmax *)

(*  Find maximas and minimas using a simple naive derivative,
    put them into two lists (one for maxs, one for mins) for later processing  *)
let local_min_max process_index signal =
    let n = arlen signal in
    let lmax = ref [] and lmin = ref [] in
    let diff = ref 0. and prevdiff = ref 0. in
    for i=0 to n-2 do
        diff := signal.(i+1) -. signal.(i) ;
        if (!prevdiff *. !diff) < 0. then begin
            if      !diff <  0. then lmax := ((process_index i, signal.(i)) :: !lmax)
            else if !diff >= 0. then lmin := ((process_index i, signal.(i)) :: !lmin)
        end ;
        prevdiff := !diff
    done ;
    (List.rev !lmax, List.rev !lmin)

(* find location and value of global minima and maxima *)
let abs_minmax arr =
    let len = arlen arr in
    (* FIXME : if len > 0, raise exception *)
    let max = ref arr.(0) in
    let min = ref arr.(0) in
    let idx_max = ref 0 in
    let idx_min = ref 0 in
    for i = 0 to len-1 do
        if arr.(i) > !max then begin
            max := arr.(i) ;
            idx_max := i
        end ;
        if arr.(i) < !min then begin
            min := arr.(i) ;
            idx_min := i
        end ;
    done ;
    (!idx_min, !min, !idx_max, !max)

let global_extrema_both = abs_minmax

(* find highest magnitude data point *)
let find_min_or_max arr =
    let len = arlen arr in
    let abs_arr = Array.map (fun x -> abs_float x) arr in
    let max = ref 0.
    and tmax = ref 0 in
    for i = 0 to len - 1 do
        if abs_arr.(i) > !max then begin
            max := abs_arr.(i) ;
            tmax := i
        end
    done ;
    (!tmax, arr.(!tmax))

let global_magnitude_max = find_min_or_max
