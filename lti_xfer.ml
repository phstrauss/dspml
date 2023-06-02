(*  one-dimensional discrete-time signal processing, with audio in mind.

    IIR filters.
    LTI SYSTEMS TRANSFER FUNCTION - polynomials and poles/zeros products
    OppSch 2nd ed. p.258 and followings

    © Philippe Strauss, 2011 - 2013
    License: see accompanying LICENSE.txt or LICENSE.md file.

    References :

        [OppSch]    : Discrete-Time Signal Processing, Oppenheim & Schafer, Prentice Hall, 3rd ed.
        [ProMan]    : Digital Signal Processing, Proakis & Manolakis, Macmillan, 2nd ed.
        [OppWill]   : Signals & Systems, Oppenheim & Willsky, Prentice Hall, Intl 2nd ed.
        [ElFilt]    : Electronic Filter Design Handbook, 4th Ed.
        [IngPro]    : Digital Signal Processing using MATLAB, Ingle & Proakis  *)


open Oneliners
open Complex
module P = Stdlib


(*  return z coeffs of 2nd order poly in descending order (or z^(-x) with x
    in increasing order), rank in tuple corresponding to x order. change of sign
    of theta due to conjugation.
    OppSch 2nd ed. eq. 5.71 p. 265 - 2012-11-8 - TESTED in t_bilin_polar.ml 2013-04-11  *)
let poly2_from_conjugate r theta =
    let rcos = 2. *. r *. (cos theta) in
    let r2 = r ** 2. in
    (* z^0 , z^(-1), z^(-2) *)
    (1., -.rcos, r2)

type quad_root_t = Real of float | Conjugate of Complex.t
exception Quadratic_roots_reals

(*  TESTED in t_bilin_polar.ml 2013-04-13  *)
let conjugate_from_poly2_x poly2 =
    let c, b, a = poly2 in
    let discr = b ** 2. -. 4. *. a *. c in
    if discr > 0. then raise Quadratic_roots_reals ;
    if discr = 0. then Real (-.b /. (2.*.a)) else Conjugate {re = (-.b)/.(2.*.a); im = (P.sqrt(4.*.a*.c -. b**2.) /. (2.*.a))}

(*  conjugate_from_poly2_x : solved over x^0, x^1, x^2, here transform over x = z^(-1);
    Complex.conj : flip imaginary sign before being flipped by inv.
    TESTED in t_bilin_polar.ml 2013-04-13  *)
let conjugate_from_poly2 poly2 =
    match conjugate_from_poly2_x poly2 with
    | Conjugate pole_c -> Conjugate (Complex.inv (Complex.conj pole_c))
    | Real pole_r -> Real (1. /. pole_r)

(*  using mathematica : z power polynomial expand/multiplication  *)

let poly3_from_poly2 poly2a poly1 =
    let a, b, c = poly2a
    and u, v = poly1 in
    let z0 = a*.u
    and z1 = b*.u +. a*.v
    and z2 = c*.u +. b*.v
    and z3 = c*.v in
    (z0, z1, z2, z3)

let poly4_from_poly2 poly2a poly2b =
    let a, b, c = poly2a
    and d, e, f = poly2b in
    let z0 = a*.d
    and z1 = b*.d +. a*.e
    and z2 = c*.d +. b*.e +. a*.f
    and z3 = c*.e +. b*.f
    and z4 = c*.f in
    (z0, z1, z2, z3, z4)

let poly6_from_poly2 poly2a poly2b poly2c =
    let a, b, c = poly2a
    and d, e, f = poly2b
    and g, h, i = poly2c in
    let z0 = a*.d*.g
    and z1 = b*.d*.g +. a*.e*.g +. a*.d *.h
    and z2 = c*.d*.g +. b*.e*.g +. a*.f*.g +. b*.d*.h +. a*.e*.h +. a*.d*.i
    and z3 = c*.e*.g +. b*.f*.g +. c*.d*.h +. b*.e*.h +. a*.f*.h +. b*.d*.i +. a*.e*.i
    and z4 = c*.f*.g +. c*.e*.h +. b*.f*.h +. c*.d*.i +. b*.e*.i +. a*.f*.i
    and z5 = c*.f*.h +. c*.e*.i +. b*.f*.i
    and z6 = c*.f*.i in
    (z0, z1, z2, z3, z4, z5, z6)



(* ***** TRANSFER FUNCTION ***** *)


(*  1st order LTI transfer function, params in polar form
    ProMan 5.2.16 - TESTED against xfer .ml tests

    ***** WARNING : ALL "*_factor1" must be used twice for conjugates pole(/zero) pairs,
                    by flipping theta sign on the second use

                    factor1 =^ 1. -. r *. e^(j*.arc) *. e^(-j*.w)  *****  *)

let xfer_factor1 r theta w =
    {re = 1. -. r *. cos (w -. theta); im = r *. sin (w -. theta)}

(*  2nd order LTI complex conjugate pair transfer function, params in polar form
    OppSch 3rd ed. eq. 5.61, 2nd ed. eq. 5.71.

    TESTED 2013-04-14 in t_xfer.ml

    Unit circle z eval of : H(jw) = 1. -. 2.*.r*.(cos theta)*.e^(-j*.w) +. r ** 2. * e^(-j*.2.*.w)

    Closed form "Expand/Collect" of the above by hand/head,
    plus Mathematica check :

    hjw = 1 - 2*r*Cos[theta]*Exp[-I*w] + r^2 *Exp[-I*2*w] ; ExpToTrig[hjw] ->

    1 - 2 r Cos[theta] Cos[w] + r^2 Cos[2 w] + 2 I r Cos[theta] Sin[w] - I r^2 Sin[2 w]

    BEWARE OP. PRECEDENCE : (cos 2.*.w) != (cos (2.*.w))  *)
let xfer_conjugate2 r theta w =
    (*let exp1 = Complex.exp {re=0.; im=(-.w)}
    and exp2 = Complex.exp {re=0.; im=(-. 2. *. w)} in
    let sub = Complex.sub (Complex.one) (Complex.mul {re = 2. *. r *. (cos theta); im = 0.} exp1) in
    Complex.add sub (Complex.mul {re = r ** 2.; im = 0.} exp2)*)
    let rcos = 2. *. r *. (cos theta) in
    let r2 = r ** 2. in
    let w2 = 2. *. w in
    {re = 1. -. rcos *. (cos w) +. r2 *. (cos w2) ;
     im = rcos *. (sin w) -. r2 *. (sin w2) }

(*  First derivative rel. w of the above, for Tribolet GRD computation :
    2 I r Cos[theta] Cos[w] - 2 I r^2 Cos[2 w] + 2 r Cos[theta] Sin[w] - 2 r^2 Sin[2 w]  *)
let xfer'_conjugate2 r theta w =
    let rcos = 2. *. r *. (cos theta) in
    let r2 = r ** 2. in
    let w2 = 2. *. w in
    {re = rcos *. (sin w) -. 2. *. r2 *. (sin w2) ;
     im = rcos *. (cos w) -. 2. *. r2 *. (cos w2) }

(*  SQUARE OF MAGNITUDE OF 1st order LTI, see OppSch p.292 3rd ed, eq. 5.64, 
    2nd ed. p.254-258 and followings - TESTED 2012-11-7  *)
let mag2_factor1 r theta w =
    1. +. r ** 2. -. 2. *. r *. cos (w -. theta)

(* TESTED 2012-11-7 *)
let mag_factor1 r theta w = P.sqrt (mag2_factor1 r theta w)

(*  overall magnitude by the multiplication of 1st order magnitude factors
    g0 = b0/a0 - untested
    poles and zeros: list of p/z in polar form, [(radius0, arc0), (radius1, arc1), ...] *)
let mag_prod1 ?g0:(g0=1.) poles zeros wn =
    let rec numdenom pz =
        match pz with
        | hd :: tl -> (mag2_factor1 (fst hd) (snd hd) wn) *. (numdenom tl)
        | [] -> 1. in
    g0 *. P.sqrt (numdenom zeros) /. P.sqrt (numdenom poles)

(*  val fold_left : ('a -> 'b -> 'a) -> 'a -> 'b list -> 'a
    List.fold_left f a [b1; ...; bn] is f (... (f (f a b1) b2) ...) bn. *)

(*  poles - zeros product
    TESTED 2012-11-7 against OppSch example  *)
let xfer_prod1 ?g0:(g0=1.) c_num d_denom w =
    let diff xk =
        let cexp = {re=0.; im=(-.w)} in
        Complex.sub Complex.one (Complex.mul xk (Complex.exp cexp)) in
    let ck_ = List.map (fun x -> diff x) c_num in
    let dk_ = List.map (fun x -> diff x) d_denom in
    let poly = Complex.div (List.fold_left Complex.mul Complex.one ck_) (List.fold_left Complex.mul Complex.one dk_) in
    Complex.mul poly {re=g0; im=0.}

(*  b/a POLYNOMIAL - order according to rank in list
    TESTED 2012-11-8 against matlab *)
(* let xfer_poly b_num a_denom w =
    let poly_exp k xk =
        let cexp = {re=0.; im=(-.w *. (foi k))} in
        Complex.mul xk (Complex.exp cexp) in  
    let bke = BatList.mapi (fun k bk -> poly_exp k bk) b_num in
    let ake = BatList.mapi (fun k ak -> poly_exp k ak) a_denom in
    Complex.div (List.fold_left Complex.add Complex.zero bke) (List.fold_left Complex.add Complex.zero ake) *)


(* ***** PHASE ***** *)


(* principal value of arg - 2012-11-7 - see parg_prod *)

let parg_factor1 r theta w =
    atan2 (r *. sin (w -. theta)) (1. -. r *. cos (w -. theta))

(*  UNTESTED 2013-04-13  *)
let parg_conjugate2 r theta w =
    let xfer = xfer_conjugate2 r theta w in
    (atan2 xfer.im xfer.re)

(* beware of accumulation sign, plus or minus, with phase/arg and GRD :
   on numerator: +, on denominator: - *)

(*  Principal value of phase - TESTED 2012-11-09 *)
let parg_prod1 c_num d_denom w =
    let principal pz =
        let r, theta = polar_from_complex pz in
        parg_factor1 r theta w in
    let principal_num = List.fold_left
        (fun acc pz ->  acc +. principal pz) 0. c_num in
    let principal_denom = List.fold_left
        (fun acc pz -> acc -. principal pz) 0. d_denom in
    (principal_num +. principal_denom)



(* ***** GROUP DELAY see OppSch 2nd ed eq. 5.55, p. 255 ***** *)


(* 2012-11-7 see grd_prod *)
let grd_factor1 r theta w =
    (r ** 2. -. r *. cos (w -. theta)) /. (mag2_factor1 r theta w)

(*  Conjugate2 GRD and GRD' gives overlong result using mathematica, giving up,
    and can use Tribolet's rels. instead  *)
let grd_conjugate2 r theta w =
    (*  Alternative Tribolet GRD numerator (mathematica) :
        -2 r (r + r^3 + r Cos[2 theta] - (1 + 3 r^2) Cos[theta] Cos[w] + r Cos[2 w])
        Denominator :
        Abs[1 + E^(-2 I w) r^2 - 2 E^(-I w) r Cos[theta]]^2  *)
    let xfer = xfer_conjugate2 r theta w in
    let xfer' = xfer'_conjugate2 r theta w in
    let mag2 = Complex.norm2 xfer in
    (xfer.re *. xfer'.im -. xfer.im *. xfer'.re) /. mag2

(*  Group delay of whole LTI for one point in frequency
    TESTED 2012-11-09 against OppSch 6th order bilin butterworth example
    Also see OppSch 2nd ed. eq. 5.56  *)
let grd_prod1 c_num d_denom d0 w =
    let grd pz =
        let r, theta = polar_from_complex pz in
        grd_factor1 r theta w in
    let grd_num = List.fold_left
        (fun grd_acc pz ->  grd_acc +. grd pz) d0 c_num in
    let grd_denom = List.fold_left
        (fun grd_acc pz -> grd_acc -. grd pz) 0. d_denom in
    (grd_num +. grd_denom)



(* ***** SECOND DERIVATIVE OF PHASE ***** *)


(* found using mathematica - UNTESTED and double check sign *)
let grd'_factor1 r theta w =
    -.((r ** 3. -. r) *. sin(w -. theta)) /. ((mag2_factor1 r theta w) ** 2.)

(* Derivative of group delay of whole LTI for one point in frequency - 2012-11-08 - UNTESTED *)
let grd'_prod1 c_num d_denom d0 w =
    let grd' pz =
        let r, theta = polar_from_complex pz in
        grd'_factor1 r theta w in
    let grd'_num = List.fold_left
        (fun grd'_acc pz ->  grd'_acc +. grd' pz) d0 c_num in
    let grd'_denom = List.fold_left
        (fun grd'_acc pz -> grd'_acc -. grd' pz) 0. d_denom in
    (grd'_num +. grd'_denom)


(*  Eventual TODO :  OppSch 2nd ed. eq. 5.50 (dbmag), 5.54 (phase of pz LTI), 5.56 (grd of pz LTI, alternate form)  *)
