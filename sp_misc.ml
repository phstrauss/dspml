(*  one-dimensional discrete-time signal processing, with audio in mind.

    Misc.

    © Philippe Strauss, 2009 - 2012; 2020
    License: see accompanying LICENSE.txt or LICENSE.md file.

    References :

        [OppSch]    : Discrete-Time Signal Processing, Oppenheim & Schafer, Prentice Hall, 3rd ed.
        [ProMan]    : Digital Signal Processing, Proakis & Manolakis, Macmillan, 2nd ed.
        [OppWill]   : Signals & Systems, Oppenheim & Willsky, Prentice Hall, Intl 2nd ed.
        [ElFilt]    : Electronic Filter Design Handbook, 4th Ed.
        [IngPro]    : Digital Signal Processing using MATLAB, Ingle & Proakis  *)


open Oneliners
open Vec_search
module Log = Logs


(* read filters coeffs from a matlab saved text file - a file with one float in sci/engineering format per line *)
let read_filter_coeffs filename =
    (* FIX: exceptions: Scanf.Scan_failure / Failure / End_of_file *)
    let scanb = Scanf.Scanning.from_file filename in
    let get_next_coeff scanb =
        Scanf.bscanf scanb "%f\n" (fun x -> x)
    and ended = ref false
    and h_list = ref [] in 
    while not !ended do
        try ( h_list := !h_list @ [ (get_next_coeff scanb) ] ; () )
        with End_of_file -> ended := true
    done ;
    Array.of_list !h_list

(* simple histogram *)
let histogram data min max nbin =
    let len_d = arlen data in
    let histo = Array.make nbin 0 in
    let histo_index = Array.init nbin (fun i -> (foi i *. (max -. min) /. foi nbin) +. min) in
    for i = 0 to (len_d-1) do
        let value = data.(i) in
        (* FIXME: raise exception when value out of bounds *)
        let idx = index_of_float_lingrid histo_index value in
        histo.(idx) <- histo.(idx) + 1 ;
    done ;
    (histo_index, histo)


(*  cross fade two arrays of same length datas in one transition band *)
let xfade dta1 dta2 idxlo idxhi =
    let len = arlen dta1 in
    assert (len = arlen dta2) ;
    let out = Array.make len 0. in
    let len_tr = idxhi - idxlo in
    let xfade_win = Windows_cosine.x_cos_fade len_tr in
    for i = 0 to (len-1) do
        if i <= idxlo then
            out.(i) <- dta1.(i) ;
        if i > idxlo && i < idxhi then
            out.(i) <- xfade_win.(1).(i-idxlo) *. dta1.(i) +. xfade_win.(0).(i-idxlo) *. dta2.(i) ; 
        if i >= idxhi then
            out.(i) <- dta2.(i) ;
    done ;
    out

let signal_fade_bothends signal xgrid xll xlh xhl xhh endsval =
    let idx_ll = index_of_freq xgrid xll
    and idx_lh = index_of_freq xgrid xlh
    and idx_hl = index_of_freq xgrid xhl
    and idx_hh = index_of_freq xgrid xhh in
    [%log info "signal_fade_bothends: ll idx=%d / %f; lh idx=%d / %f; hl idx=%d / %f; hh idx=%d / %f" idx_ll xgrid.(idx_ll) idx_lh xgrid.(idx_lh) idx_hl xgrid.(idx_hl) idx_hh xgrid.(idx_hh)] ;
    let len_tr_l = idx_lh - idx_ll
    and len_tr_h = idx_hh - idx_hl in
    let xfade_l = Windows_cosine.x_cos_fade len_tr_l
    and xfade_h = Windows_cosine.x_cos_fade len_tr_h in
    let len_m = arlen signal in
    [%log debug "signal_fade_bothends: len_m = %d; arlen xgrid = %d" len_m (arlen xgrid)] ;
    assert (len_m = arlen xgrid) ;
    let out = Array.make len_m 0. in
    for i = 0 to (len_m-1) do
        if i <= idx_ll then
            out.(i) <- endsval ;
        if i > idx_ll && i < idx_lh then
            out.(i) <- signal.(i) *. xfade_l.(0).(i-idx_ll) +. xfade_l.(1).(i-idx_ll) *. endsval ;
        if i >= idx_lh && i <= idx_hl then
            out.(i) <- signal.(i) ;
        if i > idx_hl && i < idx_hh then
            out.(i) <- signal.(i) *. xfade_h.(1).(i-idx_hl) +. xfade_h.(0).(i-idx_hl) *. endsval ;
        if i >= idx_hh then
            out.(i) <- endsval ;
    done ;
    out

(*  cross fade at both ends a real only magnitude spectrum with the value of 1.
    was named x_fade_spectrum  *)
let spectrum_fade_unity mag fgrid fll flh fhl fhh =
    signal_fade_bothends mag fgrid fll flh fhl fhh 1.


(* helpers for FIR filter design w LP simplex, remez, ..., !very specific - TESTED *)

type trigf_t = Cosine | Sine

let trig fct f harm =
    let w = 2. *. pi *. f *. (foi harm) in
    match fct with
    | Cosine -> cos w
    | Sine -> sin w

let init_trig_table fct nf ncol =
    let tbl = Array.make_matrix nf ncol 0. in
    let fstep = 0.5 /. (foi nf) in
    for i = 0 to nf - 1 do
        for j = 0 to ncol - 1 do
            tbl.(i).(j) <- trig fct (fstep *. (foi i)) j
        done ;
    done ;
    tbl
