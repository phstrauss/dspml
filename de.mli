(* exception MallocFailed
exception ParameterError *)
type strategy_t =
    Best1Exp
  | Rand1Exp
  | Best1Bin
  | Rand1Bin
  | EitherOr
  | Best2Exp
  | Rand2Exp
  | Best2Bin
  | Rand2Bin
type comparison_t = Minimize | Minimize_Magnitude | Maximize
external de_init_uniform :
  int -> int -> float array -> float array -> comparison_t -> unit
  = "ml_de_init_uniform"
external de_init_solution :
  int -> int -> float array -> float array -> comparison_t -> float array -> float array -> unit
  = "ml_de_init_solution_bytecode" "ml_de_init_solution_native"
external de_setup_std : float -> float -> strategy_t -> unit
  = "ml_de_setup_std"
external de_setup_random_fk : float -> float -> float -> float -> unit
  = "ml_de_setup_random_fk"
external de_setup_pcx :
  float -> float -> float -> strategy_t -> float -> float -> unit
  = "ml_de_setup_pcx_bytecode" "ml_de_setup_pcx_native"
external de_solve : int -> unit = "ml_de_solve"
external de_get_solution : unit -> float array = "ml_de_get_solution"
external de_finalize : unit -> unit = "ml_de_finalize"
val de_register_energyfunc :
  (int -> float -> float array -> bool * float) -> unit
