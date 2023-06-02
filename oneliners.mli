val pi : float
val pr : ('a, out_channel, unit) format -> 'a
val spr : ('a, unit, string) format -> 'a
val foi : int -> float
val iof : float -> int
val soi : int -> string
val ios : string -> int
val floori : float -> int
val ceili : float -> int
val arlen : 'a array -> int
val balen : ('a, 'b, 'c) Bigarray.Array1.t -> int
val hfind : ('a, 'b) Hashtbl.t -> 'a -> 'b
val half_down : int -> int
val half_up : int -> int
val pow : int -> int -> int
val fact : int -> int
exception NegativeNotApplicable
val log10_ex : float -> float
val polar_from_complex : Complex.t -> float * float
val complex_from_polar : float -> float -> Complex.t
val asinh : float -> float
val acosh : float -> float
(* val spr_farray : BatFloat.t BatArray.t -> string
val prf : BatFloat.t BatArray.t -> unit *)
val spr_complex : Complex.t -> string
val spr_conj : Complex.t -> string
val spr_angle : float -> string
