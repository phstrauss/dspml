(*  one-dimensional discrete-time signal processing, with audio in mind.

    Functions to soft-clip signals with a parabola or atan segment having first derivative
    equal to a unitary slope at the junction point.

    © Philippe Strauss, 2012
    License: see accompanying LICENSE.txt or LICENSE.md file.

    References :

        [OppSch]    : Discrete-Time Signal Processing, Oppenheim & Schafer, Prentice Hall, 3rd ed.
        [ProMan]    : Digital Signal Processing, Proakis & Manolakis, Macmillan, 2nd ed.
        [OppWill]   : Signals & Systems, Oppenheim & Willsky, Prentice Hall, Intl 2nd ed.
        [ElFilt]    : Electronic Filter Design Handbook, 4th Ed.
        [IngPro]    : Digital Signal Processing using MATLAB, Ingle & Proakis  *)

(*  compute the quadratic ax^2 + bx + c coefficients of the upper soft-clipping part of the curve,
    the positive part. actually, only coeff "a" remain, we'll proceed by axis translation  *)

module Log = Logs


type quad_t = {xth: float; xc: float; yth: float; yc: float; a: float}

let quad_a xthp xcp =
    let xj = xcp -. xthp in
    let a = (-.1. /. (2. *. xj)) in
    (* LOG "quad_a: xj=%f, a=%f" xj a LEVEL TRACE ; *)
    (a, 0.5 *. xj)

exception AxisTranslationBeyoundBounds of float

type headroom_t = Headroom_x of float | Headroom_y of float

let clip_level head =
    match head with
    | Headroom_x l -> l
    | Headroom_y l -> 2. *. l

let quad_params x0th head =
    let x0c = x0th +. clip_level head in
    let asq, yj = quad_a x0th x0c in
    [%log debug "quad_params: a=%f, xth=%f, xc=%f, yth=%f, yc=%f" asq x0th x0c yj (x0th +. yj)] ;
    {xth = x0th; xc = x0c; yth = yj; yc = x0th +. yj; a = asq}


module QuadUpper = struct

    let create = quad_params

    (*  x axis translation for the upper soft-clipped part *)
    let axis_trans_x_hi x0 upst =
        if x0 >= upst.xth && x0 <= upst.xc then x0 -. upst.xc else raise (AxisTranslationBeyoundBounds x0)
    let axis_trans_y_hi y1 upst =
        let y0 = upst.yc +. y1 in
        if y0 > upst.yc then raise (AxisTranslationBeyoundBounds y0) ;
        y0

    let square_clip_high x0 upst =
        if x0 < upst.xth then x0
        else if x0 >= upst.xth && x0 <= upst.xc then begin
            let xq = axis_trans_x_hi x0 upst in
            let yq = upst.a *. xq ** 2. in
            axis_trans_y_hi yq upst
        end else upst.yc

end (* Upper *)


module QuadLower = struct

    let create = quad_params

    let axis_trans_x_lo x0 lost =
        if x0 <= lost.xth && x0 >= lost.xc then x0 -. lost.xc else raise (AxisTranslationBeyoundBounds x0)

    let axis_trans_y_lo y1 lost =
        let y0 = lost.yc +. y1 in
        if y0 < lost.yc then raise (AxisTranslationBeyoundBounds y0) ;
        y0

    let square_clip_low x0 lost =
        if x0 > lost.xth then x0
        else if x0 <= lost.xth && x0 >= lost.xc then begin
            let xq = axis_trans_x_lo x0 lost in
            let yq = lost.a *. xq ** 2. in
            axis_trans_y_lo yq lost
        end else lost.yc

end (* Lower *)


module QuadBoth = struct

    type t = quad_t * quad_t

    let create xthm headm xthp headp =
        (quad_params xthm headm, quad_params xthp headp)

    let square_clip x0 tupst =
        if x0 > 0. then QuadUpper.square_clip_high x0 (snd tupst)
        else if x0 < 0. then QuadLower.square_clip_low x0 (fst tupst)
        else 0.

end (* Both *)


module QuadSymetric = struct

    type t = quad_t * quad_t

    let sym_headroom head =
        match head with
        | Headroom_x f -> (Headroom_x (-.f))
        | Headroom_y f -> (Headroom_y (-.f))

    let create xthp headp =
        (quad_params (-.xthp) (sym_headroom headp), quad_params xthp headp)

    let square_clip x0 tupst =
        if x0 > 0. then QuadUpper.square_clip_high x0 (snd tupst)
        else if x0 < 0. then QuadLower.square_clip_low x0 (fst tupst)
        else 0.

end (* Symetric *)


module Atan = struct

    let soft_clip x linlo linhi =
        if x >= linlo && x <= linhi then x
        else if x > linhi then linhi +. (Stdlib.atan (x-.linhi))
        else if x < linlo then linlo +. (Stdlib.atan (x-.linlo))
        else 0. 

end (* module Atan *)

