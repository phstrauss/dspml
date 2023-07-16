(*  vector signal processing, mostly related to straightliner/audiofocus inverse convolution
    Phase/GRD processing

    © Philippe Strauss, 2009 - 2012, 2013, 2016  *)


open Oneliners
module Log = Logs


(*  despite the name can gives good result, on a smooth curve,
    zeropad the ir quiet a bit before use.
    WW: naive derivative means prone to offset by 1/2 sample on the x axis  *)
let naive spectrum fgrid =
    let phase         = Array.map (fun x -> Complex.arg x) spectrum in
    let len           = arlen phase in
    assert (arlen fgrid = len) ;
    let unwrapped     = Array.make len 0. in
    let gdelay        = Array.make len 0. in
    let nwrap_arr     = Array.make len 0 in
    let thresh_hi     = pi in
    let thresh_low    = (-.pi) in
    let df = (fgrid.(len-1) -. fgrid.(0)) /. (foi len) in
    let dw = df *. 2. *. pi in
    let nwrap = ref 0 in
    let grd = ref 0. in
    for i = 0 to (len-2) do
        let dph = phase.(i+1) -. phase.(i) in
        if dph > thresh_hi then (
            nwrap := !nwrap - 1 ;
            grd := -. (dph -. 2. *. pi) /. dw ;
        ) else if dph < thresh_low then (
            nwrap := !nwrap + 1 ;
            grd := -. (dph +. 2. *. pi) /. dw ;
        ) else (
            grd := -. dph /. dw ;
        ) ;
        (* N.B.: trying to equalize somehow x axis offset rel. to grd_eq_zero, and it works! *)
        nwrap_arr.(i) <- !nwrap ;
        unwrapped.(i) <- phase.(i+1) +. (foi !nwrap *. 2. *. pi) ;
        gdelay.(i) <- !grd ;
    done ;
    (* wrapped is length len, the others would be len-1 (more formally) *)
    nwrap_arr.(len-1) <- nwrap_arr.(len-2) ;
    unwrapped.(len-1) <- unwrapped.(len-2) ;
    gdelay.(len-1) <- gdelay.(len-2) ;
    (phase, unwrapped, gdelay)
