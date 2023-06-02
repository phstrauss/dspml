type fade_direction = In | Out
val cos_fade : fade_direction -> int -> int -> float
val x_cos_fade : int -> float array array
val x_cos_fade_window : int -> int -> int -> float array array
val tukey : int -> int -> int -> int -> int -> float array
val tukey_time :
  int -> int -> float -> float -> float -> float -> float array
val hann_coeff : int -> int -> float
val hann_window : int -> float array
val hamming_coeff : int -> int -> float
val hamming_window : int -> float array
