val integrate_trapezoidal : ?scale:float -> int -> float array -> float
val deriv_naive : float array -> float array
type coeffs_t = Coeffs of float array | File of string
val filter_fir : coeffs_t -> float array -> Complex.t array * Complex.t array
val deriv_fir : coeffs_t -> float array -> Complex.t array * Complex.t array
val hilbert_xform_fir :
  coeffs_t -> float array -> Complex.t array * Complex.t array
