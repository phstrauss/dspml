val butterworth_pole_arc : int -> int -> float
val butterworth_pole : ?kwc:float -> int -> int -> Complex.t
val denom_butterworth_mag : float -> float -> int -> float
val butterworth_mag_lp : float -> float -> int -> float
val butterworth_mag_hp : float -> float -> int -> float
val lr_mag_lp : float -> float -> int -> float
val lr_mag_hp : float -> float -> int -> float
val cheby_theta_m : int -> int -> float
val cheby_pole_raw_epsilon : float -> int -> int -> Complex.t
val cheby1_pole : float -> int -> int -> Complex.t
val cheby2_zero : int -> int -> Complex.t
val cheby2_pole : float -> int -> int -> Complex.t
val bessel_ak : int -> int -> int
val bessel_poly : int -> float array
