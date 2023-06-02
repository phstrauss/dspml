(* exception ConvergenceFailed
exception MallocFailed of string *)
type filter_kind_t = Standard | Differentiator | Hilbert
type symmetry_t = Even | Odd
external remez_from_grid :
  int ->
  float array ->
  float array ->
  float array -> int -> filter_kind_t -> float -> int * float array
  = "ml_remez_from_grid_bytecode" "ml_remez_from_grid"
external remez_from_bands :
  int ->
  int ->
  float array ->
  float array -> float array -> filter_kind_t -> float -> int * float array
  = "ml_remez_from_bands_bytecode" "ml_remez_from_bands"
