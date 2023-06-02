(*  one-dimensional discrete-time signal processing, with audio in mind.

    About finding a value inside a vector

    © Philippe Strauss, 2009 - 2012, 2015
    License: see accompanying LICENSE.txt or LICENSE.md file.

    References :

        [OppSch]    : Discrete-Time Signal Processing, Oppenheim & Schafer, Prentice Hall, 3rd ed.
        [ProMan]    : Digital Signal Processing, Proakis & Manolakis, Macmillan, 2nd ed.
        [OppWill]   : Signals & Systems, Oppenheim & Willsky, Prentice Hall, Intl 2nd ed.
        [ElFilt]    : Electronic Filter Design Handbook, 4th Ed.
        [IngPro]    : Digital Signal Processing using MATLAB, Ingle & Proakis  *)


open Oneliners
module Log = Logs


(* find the index of a value in an ascending ordered vector of floats *)
let index_of_float_exhaustive grid value =
    (* FIXME raise exception on value out of grid *)
    let len = arlen grid in
    (* if Log.threshold_level = Log.Trace then begin
        for i = 0 to len-1 do
            LOG "index_of_float_exhaustive grid.(%d) = %f" i grid.(i) LEVEL TRACE ;
        done ;
    end ; *)
    let idx = ref (len-1) in
    while grid.(!idx) > value do
        idx := !idx - 1 ;
    done ;
    [%log debug "index_of_float_exhaustive: found index %d for value %f" !idx value];
    !idx

exception IndexOfFloatBoundary

(* find the index of a value in an uniformly sampled ascending ordered vector of floats *)
let index_of_float_lingrid grid value =
    (* FIXME: handle value out of grid *)
    let len_g = arlen grid in
    let index =
        if len_g > 2 then
            let step = grid.(1) -. grid.(0) in
            floori ((value -. grid.(0)) /. step)
        else
            index_of_float_exhaustive grid value
    in
    let retindex = if index < len_g - 1 then begin
        let eps_m = value -. grid.(index)
        and eps_p = grid.(index+1) -. value in
        if eps_m < eps_p then index else index+1
    end else index in
    [%log debug "index_of_float_lingrid: value=%f; retidx=%d; grid.(retindex)=%f" value retindex grid.(retindex)];
    retindex

let index_of_freq fgrid f =
    index_of_float_lingrid fgrid f

(*  binsearch result sample :
    
    let farr = [|0.12; 1.24; 1.43; 2.57; 3.09; 4.28; 4.29; 4.93; 7.09; 7.12; 9.12111; 9.85|] in
    let pprint_binsearch value arr =
        match (Sig_proc.binsearch arr value) with
        | Sig_proc.Out_of_interval -> Printf.printf "value out of interval\n" ;
        | Sig_proc.Exact idx -> Printf.printf "exact match: value %f at index %d\n" value idx ;
        | Sig_proc.Bounded (lo, hi) -> Printf.printf "value %f between lo index %d and hi index %d\n" value lo hi in
    pprint_binsearch 5. farr ;
    pprint_binsearch 7.1 farr ;
    pprint_binsearch 4.93 farr ;
    
    gives:
    
    value 5.000000 between lo index 7 and hi index 8
    value 7.100000 between lo index 8 and hi index 9
    exact match: value 4.930000 at index 7 *)

type binsearch_result = Exact of int | Bounded of int * int | Out_of_interval
			   
let binsearch sorted value =
    let len_s = arlen sorted in
    let idx_min = ref 0 in
    let idx_max = ref (len_s - 1) in
    let idx_mid = ref 0 in
    let search = ref true in
    let result = ref Out_of_interval in
    assert(sorted.(len_s -1) > sorted.(0)) ;
    while !search do
        idx_mid := (!idx_min + (!idx_max - !idx_min) / 2) ;
        if value > sorted.(!idx_mid) then begin
            idx_min := !idx_mid + 1
        end else begin
            idx_max := !idx_mid - 1 ;
        end ;
        if sorted.(!idx_mid) = value then begin
            search := false ;
            result := Exact !idx_mid ;
        end else if (!idx_min > !idx_max) then begin
            search := false ;
            result := Bounded (!idx_max, !idx_min) ;
        end ;
    done ;
    !result

exception AboveFrom_Not_found

let above_from_begin (compact: float array) threshold =
    let i = ref 0 in
    let len = arlen compact in
    while compact.(!i) < threshold do
        incr i ;
        if !i >= (len-1) then begin
	    raise AboveFrom_Not_found
	end
    done ;
    (!i, compact.(!i))

let above_from_end (compact: float array) threshold =
    let len = arlen compact in
    let i = ref (len-1) in
    while compact.(!i) < threshold do
        decr i ;
        if !i <= 0 then begin
	    raise AboveFrom_Not_found
	end
    done ;
    (!i, compact.(!i))

