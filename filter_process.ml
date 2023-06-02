(*  one-dimensional discrete-time signal processing, with audio in mind.

    IIR filters design/implementation helpers
    IIR processing - Direct Forms

    © Philippe Strauss, 2011, 2012
    License: see accompanying LICENSE.txt or LICENSE.md file.

    References :

        [OppSch]    : Discrete-Time Signal Processing, Oppenheim & Schafer, Prentice Hall, 3rd ed.
        [ProMan]    : Digital Signal Processing, Proakis & Manolakis, Macmillan, 2nd ed.
        [OppWill]   : Signals & Systems, Oppenheim & Willsky, Prentice Hall, Intl 2nd ed.
        [ElFilt]    : Electronic Filter Design Handbook, 4th Ed.
        [IngPro]    : Digital Signal Processing using MATLAB, Ingle & Proakis  *)


(* Direct Form 2 *)

(*  BIG FAT WARNING
    YOU NEED TO FLIP THE SIGN OF A/DENOM COEFFS BEFORE USING *)

let invsign coeffs =
    Array.map (fun x -> -.x) coeffs

let iir1_df2 b a1 z samples l =
    for i = 0 to l-1 do
        z.(1) <- z.(0) ;
        z.(0) <- samples.(i) +. a1 *. z.(1) ;
        samples.(i) <- b.(0) *. z.(0) +. b.(1) *. z.(1) ;
    done

let iir2_df2 b a z samples l =
    for i = 0 to l-1 do
        z.(2) <- z.(1) ; z.(1) <- z.(0) ;
        (* a[0] is assumed to be eq. to 1.0, even if never accessed. *)
        z.(0) <- samples.(i) +. a.(1) *. z.(1) +. a.(2) *. z.(2) ;
        samples.(i) <- b.(0) *. z.(0) +. b.(1) *. z.(1) +. b.(2) *. z.(2) ;
    done

let iir3_df2 b a z samples l =
    for i = 0 to l-1 do
        z.(3) <- z.(2) ; z.(2) <- z.(1) ; z.(1) <- z.(0) ;
        z.(0) <- samples.(i) +. a.(1) *. z.(1) +. a.(2) *. z.(2) +. a.(3) *. z.(3) ;
        samples.(i) <- b.(0) *. z.(0) +. b.(1) *. z.(1) +. b.(2) *. z.(2) +. b.(3) *. z.(3);
    done

let iir4_df2 b a z samples l =
    for i = 0 to l-1 do
        z.(4) <- z.(3) ; z.(3) <- z.(2) ; z.(2) <- z.(1) ; z.(1) <- z.(0) ;
        z.(0) <- samples.(i) +. a.(1) *. z.(1) +. a.(2) *. z.(2) +. a.(3) *. z.(3) +. a.(4) *. z.(4);
        samples.(i) <- b.(0) *. z.(0) +. b.(1) *. z.(1) +. b.(2) *. z.(2) +. b.(3) *. z.(3) +. b.(4) *. z.(4) ;
    done


(* Direct Form 1 *)

let iir1_df1 b a1 z1 z2 samples l =
    for i = 0 to l-1 do
        z1.(1) <- z1.(0) ; z1.(0) <- samples.(i) ;
        z2.(1) <- z2.(0) ;
        samples.(i) <- b.(0) *. z1.(0) +. b.(1) *. z1.(1) +.
            a1 *. z2.(1) ;
        z2.(0) <- samples.(i) ;
    done

let iir2_df1 b a z1 z2 samples l =
    for i = 0 to l-1 do
        z1.(2) <- z1.(1) ; z1.(1) <- z1.(0) ; z1.(0) <- samples.(i) ;
        z2.(2) <- z2.(1) ; z2.(1) <- z2.(0) ;
        samples.(i) <- b.(0) *. z1.(0) +. b.(1) *. z1.(1) +. b.(2) *. z1.(2) +.
            a.(1) *. z2.(1) +. a.(2) *. z2.(2) ;
        z2.(0) <- samples.(i) ;
    done

let iir3_df1 b a z1 z2 samples l =
    for i = 0 to l-1 do
        z1.(3) <- z1.(2) ; z1.(2) <- z1.(1) ; z1.(1) <- z1.(0) ; z1.(0) <- samples.(i) ;
        z2.(3) <- z2.(2) ; z2.(2) <- z2.(1) ; z2.(1) <- z2.(0) ;
        samples.(i) <- b.(0) *. z1.(0) +. b.(1) *. z1.(1) +. b.(2) *. z1.(2) +. b.(3) *. z1.(3) +.
            a.(1) *. z2.(1) +. a.(2) *. z2.(2) +. a.(3) *. z2.(3) ;
        z2.(0) <- samples.(i) ;
    done

let iir4_df1 b a z1 z2 samples l =
    for i = 0 to l-1 do
        z1.(4) <- z1.(3) ; z1.(3) <- z1.(2) ; z1.(2) <- z1.(1) ; z1.(1) <- z1.(0) ; z1.(0) <- samples.(i) ;
        z2.(4) <- z2.(3) ; z2.(3) <- z2.(2) ; z2.(2) <- z2.(1) ; z2.(1) <- z2.(0) ;
        samples.(i) <- b.(0) *. z1.(0) +. b.(1) *. z1.(1) +. b.(2) *. z1.(2) +. b.(3) *. z1.(3) +. b.(4) *. z1.(4) +.
            a.(1) *. z2.(1) +. a.(2) *. z2.(2) +. a.(3) *. z2.(3) +. a.(4) *. z2.(4) ;
        z2.(0) <- samples.(i) ;
    done