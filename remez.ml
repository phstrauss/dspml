(*  Jake Janovetz C implementation of remez exchange (or Parks-McClellan) algorithm for
    symetrical FIR filter design
    Objective Caml binding - C stubs file
    (c) Philippe Strauss, philippe@strauss-acoustics.ch, January 2012  *)


(* exception ConvergenceFailed
exception MallocFailed of string

let _ = Callback.register_exception "exn_Convergence" ConvergenceFailed
let _ = Callback.register_exception "exn_Malloc" (MallocFailed "") *)


type filter_kind_t = Standard | Differentiator | Hilbert
type symmetry_t = Even | Odd


external remez_from_grid: int ->
    float array -> float array -> float array ->
    int -> filter_kind_t -> float -> int * float array = "ml_remez_from_grid_bytecode" "ml_remez_from_grid"

external remez_from_bands: int -> int ->
    float array -> float array -> float array ->
    filter_kind_t -> float -> int * float array = "ml_remez_from_bands_bytecode" "ml_remez_from_bands"