val index_of_float_exhaustive : float array -> float -> int
exception IndexOfFloatBoundary
val index_of_float_lingrid : float array -> float -> int
val index_of_freq : float array -> float -> int
type binsearch_result = Exact of int | Bounded of int * int | Out_of_interval
val binsearch : 'a array -> 'a -> binsearch_result
exception AboveFrom_Not_found
val above_from_begin : float array -> float -> int * float
val above_from_end : float array -> float -> int * float
