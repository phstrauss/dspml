(*  hsv.ml
    Hue-Saturation-Value color computation
    © Philippe Strauss  *)

(* functions to generate an HSV color palette
   source: http://fr.wikipedia.org/wiki/Teinte_Saturation_Valeur *)

(* a single color channel, like green, outputs between 0. and 1. *)
let color_channel theta =
    (* theta: a hue angle in degree *)
    let cdiv = floor(theta /. 360.) in
    let ctheta = theta -. (cdiv *. 360.) in
    if (ctheta >= 0.) && (ctheta < 60.) then ctheta /. 60.
    else if (ctheta >= 60.) && (ctheta < 180.) then 1.
    else if (ctheta >= 180.) && (ctheta < 240.) then 1. -. ((ctheta -. 180.) /. 60.)
    else 0.

(* compute hue only, with max saturation and luminosity *)
let hue2rgb theta =
    let r = color_channel (theta-.240.)
    and g = color_channel theta
    and b = color_channel (theta-.120.) in
    (r, g, b)

(* Hue Saturation Value to RGB conversion *)
let hsv2rgb h s v =
    let r1, g1, b1 = hue2rgb h in
    let ofs = v *. (1. -. s) in
    let sc = v -. ofs in
    (r1 *. sc +. ofs, g1 *. sc +. ofs, b1 *. sc +. ofs)

