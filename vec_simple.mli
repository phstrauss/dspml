val creal : Complex.t array -> float array
val cimag : Complex.t array -> float array
val cmagnitude : Complex.t array -> float array
val rpower : float array -> float array
val cpower : Complex.t array -> float array
val cargument : Complex.t array -> float array
val scale : float -> float array -> float array
val eps_db_pow : float
val eps_db_field : float
val log_db10 : float array -> float array
val log_db20 : float array -> float array
val cpower_db : Complex.t array -> float array
val rfield_db : float array -> float array
val cenergy : Complex.t array -> float
val energy : float array -> float
val cenergy2 : float -> Complex.t array -> unit
val energy2 : float -> float array -> unit
val nodc : 'a array -> 'a array
val tgrid : int -> int -> float array
val fgrid : int -> int -> float array
val fgrid2 : int -> int -> float array
val fgrid_half : int -> int -> float array
val fgrid_log : int -> int -> float array
val fgrid_half_log : int -> int -> float array
val idx2time : int -> int -> float
val energy_time_curve : float array -> float array * int * float
val split_nyquist : float array -> int -> float array
val splice_nyquist : 'a array -> int -> 'a array
val even_from_single : float array -> bool -> float array
val odd_from_single : float array -> bool -> float array
val single_from_even : float array -> float array
val vec_field_avg : float array -> float array -> float -> float -> float
val vec_rms_avg : float array -> float array -> float -> float -> float
val spectrum_normalize :
  float array -> float array -> float -> float -> float array
