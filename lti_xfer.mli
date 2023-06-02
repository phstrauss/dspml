val poly2_from_conjugate : float -> float -> float * float * float
type quad_root_t = Real of float | Conjugate of Complex.t
exception Quadratic_roots_reals
val conjugate_from_poly2_x : float * float * float -> quad_root_t
val conjugate_from_poly2 : float * float * float -> quad_root_t
val poly3_from_poly2 :
  float * float * float -> float * float -> float * float * float * float
val poly4_from_poly2 :
  float * float * float ->
  float * float * float -> float * float * float * float * float
val poly6_from_poly2 :
  float * float * float ->
  float * float * float ->
  float * float * float ->
  float * float * float * float * float * float * float
val xfer_factor1 : float -> float -> float -> Complex.t
val xfer_conjugate2 : float -> float -> float -> Complex.t
val xfer'_conjugate2 : float -> float -> float -> Complex.t
val mag2_factor1 : float -> float -> float -> float
val mag_factor1 : float -> float -> float -> float
val mag_prod1 :
  ?g0:float -> (float * float) list -> (float * float) list -> float -> float
val xfer_prod1 :
  ?g0:float -> Complex.t list -> Complex.t list -> float -> Complex.t
(* val xfer_poly : Complex.t list -> Complex.t list -> float -> Complex.t *)
val parg_factor1 : float -> float -> float -> float
val parg_conjugate2 : float -> float -> float -> float
val parg_prod1 : Complex.t list -> Complex.t list -> float -> float
val grd_factor1 : float -> float -> float -> float
val grd_conjugate2 : float -> float -> float -> float
val grd_prod1 : Complex.t list -> Complex.t list -> float -> float -> float
val grd'_factor1 : float -> float -> float -> float
val grd'_prod1 : Complex.t list -> Complex.t list -> float -> float -> float
