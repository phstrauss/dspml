(*  Differential evolution solver heurisitc
    Lester E. Godwin
    Partially rewrote, enhanced (EitherOr and Probabilist Parent Centric evolution strategies,
    ocaml binding) plus maintained by Philippe Strauss, philippe@strauss-acoustics.ch. Spring 2012.  *)


open Common
open De

(* external show : 'a -> string = "%show" *)
let prf floats = Printf.printf "%s\n" (sprint_farray floats)


let epsilon = 1e-8

(* Rosenbrock's function in 2D *)
let rosenbrock_2D ndim best_energy trial =
    let res = (1. -. trial.(0)) ** 2. +. 100. *. (trial.(1) -. trial.(0) ** 2.) ** 2. in
    if (abs_float(res) <= epsilon) then (true, res) else (false, res)

(* Rastrigin's function as found on wikipedia, minima of 0.0 at 0.0 in all dims *)
let rastrigin ndim best_energy trial =
    let sum = ref 0. in
    for i = 0 to ndim-1 do
        sum := !sum +. trial.(i) ** 2. -. 10. *. cos (2. *. pi *. trial.(i))
    done ;
    let ret = !sum +. (foi ndim) *. 10. in
    if (abs_float ret <= epsilon) then (true, ret) else (false, ret)

let schwefel ndim best_energy trial =
    let sum = ref 0. in
    for i = 0 to ndim-1 do
        sum := !sum +. (-1.0) *. trial.(i) *. sin (sqrt (abs_float trial.(i))) ;
    done ;
    let ret = 418.9828872724 *. (foi ndim) +. !sum in
    (false, ret)

let ackley ndim best_energy trial =
    let sumx2 = ref 0.0
    and sumcos = ref 0.0
    and c = 2.0 *. pi
    and a = 20.0
    and b = 0.2 in
    for i=0 to ndim-1 do
        sumx2 := !sumx2 +. trial.(i) ** 2.0 ;
        sumcos := !sumcos +. cos (c *. trial.(i)) ;
    done ;
    let n = foi ndim in
    sumx2 := !sumx2 /. n ;
    sumcos := !sumcos /. n ;
    let ret = -.a *. exp (-.b *. (sqrt !sumx2)) -. (exp !sumcos) +. a +. (exp 1.0) in
    if (abs_float(ret) <= epsilon) then (true, ret) else (false, ret)

let griewangk ndim best_energy trial =
    let sumx2 = ref 0.0
    and prodcos = ref 0.0 in
    for i=0 to ndim-1 do
        prodcos := !prodcos *. (cos (trial.(i)) /. sqrt (foi (i+1))) ;
        sumx2 := !sumx2 +. trial.(i) ** 2.0 ;
    done ;
    let ret = 1.0 +. (!sumx2 /. 4000.0) -. !prodcos in
    if (abs_float ret) <= (epsilon +. 1.0) then (true, ret) else (false, ret)    


let npop = 350
let scale = 0.7
let xover_prob = 0.6
let f_rand_fact = 0.001
let f_pop_fact = 0.0
let pcx_prob = 0.07
let std_strategy = Best1Exp
let pcx_strategy = Rand1Exp
let max_generations = 100000
let sigma2a = 0.7
let sigma2b = 0.3


let () =

    (* start gently : rosenbrock 2D *)
    let ndim = 2 in
    let min = Array.make ndim (-5.) in
    let max = Array.make ndim 5. in
    de_init_uniform ndim npop min max Minimize;
    de_setup_pcx pcx_prob scale xover_prob pcx_strategy sigma2a sigma2b ;
    de_setup_random_fk f_rand_fact f_pop_fact f_rand_fact f_pop_fact ;
    de_register_energyfunc rosenbrock_2D ;
    de_solve 100 ;
    (* print_endline (show (de_get_solution ())) ; *)
    prf (de_get_solution ()) ;
    de_finalize () ;

    (* ackley *)
    let ndim = 100 in
    let min = Array.make ndim (-32.768) in
    let max = Array.make ndim 32.768 in
    de_init_uniform ndim npop min max Minimize;
    de_setup_pcx 0.01 scale xover_prob pcx_strategy sigma2a sigma2b ;
    de_setup_random_fk f_rand_fact f_pop_fact f_rand_fact f_pop_fact ;
    de_register_energyfunc ackley ;
    de_solve 10000 ;
    (* print_endline (show (de_get_solution ())) ; *)
    prf (de_get_solution ()) ;
    de_finalize () ;

    let ndim = 100 in
    let min = Array.make ndim (-32.768) in
    let max = Array.make ndim 32.768 in
    de_init_uniform ndim npop min max Minimize;
    de_setup_std scale xover_prob std_strategy ;
    de_setup_random_fk f_rand_fact f_pop_fact f_rand_fact f_pop_fact ;
    de_register_energyfunc ackley ;
    de_solve 10000 ;
    (* print_endline (show (de_get_solution ())) ; *)
    prf (de_get_solution ()) ;
    de_finalize () ;

    (* griewangk *)
    let ndim = 100 in
    let min = Array.make ndim (-600.0) in
    let max = Array.make ndim 600.0 in
    de_init_uniform ndim npop min max Minimize;
    de_setup_pcx pcx_prob scale xover_prob pcx_strategy sigma2a sigma2b ;
    de_setup_random_fk f_rand_fact f_pop_fact f_rand_fact f_pop_fact ;
    de_register_energyfunc griewangk ;
    de_solve 10000 ;
    (* print_endline (show (de_get_solution ())) ; *)
    prf (de_get_solution ()) ;
    de_finalize () ;

    let ndim = 100 in
    let min = Array.make ndim (-600.0) in
    let max = Array.make ndim 600.0 in
    de_init_uniform ndim npop min max Minimize;
    de_setup_std scale xover_prob std_strategy ;
    de_setup_random_fk f_rand_fact f_pop_fact f_rand_fact f_pop_fact ;
    de_register_energyfunc griewangk ;
    de_solve 10000 ;
    (* print_endline (show (de_get_solution ())) ; *)
    prf (de_get_solution ()) ;
    de_finalize () ;

    (* rastrigin *)
    let ndim = 100 in
    let min = Array.make ndim (-5.12) in
    let max = Array.make ndim 5.12 in
    de_init_uniform ndim npop min max Minimize;
    de_setup_pcx pcx_prob scale xover_prob pcx_strategy sigma2a sigma2b ;
    de_setup_random_fk f_rand_fact f_pop_fact f_rand_fact f_pop_fact ;
    de_register_energyfunc rastrigin ;
    de_solve max_generations ;
    (* print_endline (show (de_get_solution ())) ; *)
    prf (de_get_solution ()) ;
    de_finalize () ;

    let ndim = 100 in
    let min = Array.make ndim (-5.12) in
    let max = Array.make ndim 5.12 in
    de_init_uniform ndim npop min max Minimize;
    de_setup_std scale xover_prob std_strategy ;
    de_setup_random_fk f_rand_fact f_pop_fact f_rand_fact f_pop_fact ;
    de_register_energyfunc rastrigin ;
    de_solve max_generations ;
    (* print_endline (show (de_get_solution ())) ; *)
    prf (de_get_solution ()) ;
    de_finalize () ;

    (* schwefel *)
    let ndim = 100 in
    let min = Array.make ndim (-500.0) in
    let max = Array.make ndim 500.0 in
    de_init_uniform ndim npop min max Minimize;
    de_setup_std scale xover_prob Rand1Exp ;
    de_setup_random_fk f_rand_fact f_pop_fact f_rand_fact f_pop_fact ;
    de_register_energyfunc schwefel ;
    de_solve 10000 ;
    (* print_endline (show (de_get_solution ())) ; *)
    prf (de_get_solution ()) ;
    de_finalize () ;

