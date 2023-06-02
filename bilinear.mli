val omega_warp : float -> float -> float
val omega_dt2ct : float -> float -> float
val allpoles_lp : float -> Complex.t -> float array * float array
val allpoles_hp : float -> Complex.t -> float array * float array
val norm2_lp : float -> Complex.t -> float
val norm2_hp : float -> Complex.t -> float
val norm1_lp : float -> float -> float
val pz_bilin : float -> Complex.t -> Complex.t
val pz_bilin_polar : Complex.t -> float -> float * float
