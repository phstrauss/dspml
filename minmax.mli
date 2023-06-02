val local_min_max :
  (int -> 'a) -> float array -> ('a * float) list * ('a * float) list
val abs_minmax : 'a array -> int * 'a * int * 'a
val global_extrema_both : 'a array -> int * 'a * int * 'a
val find_min_or_max : float array -> int * float
val global_magnitude_max : float array -> int * float
