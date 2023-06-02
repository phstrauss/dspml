val read_filter_coeffs : Scanf.Scanning.file_name -> float array
val histogram :
  float array -> float -> float -> int -> float array * int array
val xfade : float array -> float array -> int -> int -> float array
val signal_fade_bothends :
  float array ->
  float array -> float -> float -> float -> float -> float -> float array
val spectrum_fade_unity :
  float array ->
  float array -> float -> float -> float -> float -> float array
type trigf_t = Cosine | Sine
val trig : trigf_t -> float -> int -> float
val init_trig_table : trigf_t -> int -> int -> float array array
