/*  Differential evolution solver heurisitc
    Lester E. Godwin
    Partially rewrote, plain ANSII C port,
    enhanced (EitherOr and Probabilist Parent Centric evolution strategies,
    ocaml binding) plus maintained by Philippe Strauss, philippe@strauss-engineering.ch. Spring 2012.  */


#include <stdio.h>
#include <stdbool.h>
#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/fail.h>
#include <dspc/de.h>
#include <dspc/perr.h>


#define BEST1_EXP           0
#define RAND1_EXP           1
#define BEST1_BIN           2
#define RAND1_BIN           3
#define EITHER_OR           4
#define BEST2_EXP           5
#define RAND2_EXP           6
#define BEST2_BIN           7
#define RAND2_BIN           8


/* globals */
de_t prob_data;
// value ml_trial;


/* void init_ml_trial(int ndim) { */
/*     ml_trial = caml_alloc(ndim * Double_wosize, Double_array_tag); */
/*     caml_register_global_root(&ml_trial); */
/* } */


double ml_energy_function(int n, double bestEnergy, double *trial, bool *bAtSolution) {
    // bug with the following line : easily lost somewhere between caml_callback3
    // and entering ocaml energy function
    // CAMLlocal4(ml_ret, ml_found, ml_energy, ml_trial);
    value ml_ret, ml_found, ml_energy, ml_trial;
    int i;
    static value *closure_energy_f = NULL;

    // make global somehow ?
   	ml_trial = caml_alloc(n * Double_wosize, Double_array_tag);

    for (i = 0; i < n; i++)
            Store_double_field(ml_trial, i, trial[i]);
 
    if (closure_energy_f == NULL)
        closure_energy_f = caml_named_value("EnergyFunction");

    //dprintf("ml_energy_function: before caml_callback3\n");
    ml_ret = caml_callback3(*closure_energy_f, Val_int(n), caml_copy_double(bestEnergy), ml_trial);
    //dprintf("ml_energy_function: after caml_callback3\n");

    ml_found = Field(ml_ret, 0);
    ml_energy = Field(ml_ret, 1);

    *bAtSolution = Bool_val(ml_found);

    return (Double_val(ml_energy));	
}


CAMLprim value ml_de_init_uniform(value ml_ndim, value ml_npop, value ml_min, value ml_max, value ml_comparison) {
    CAMLparam5(ml_ndim, ml_npop, ml_min, ml_max, ml_comparison);
    int ndim, npop, comparison, len_min, len_max, i;
    double *min, *max;

    ndim = Int_val(ml_ndim);
    npop = Int_val(ml_npop);
    comparison = Int_val(ml_comparison);

    //init_ml_trial(ndim);

    len_min = Wosize_val(ml_min) / Double_wosize;
    len_max = Wosize_val(ml_max) / Double_wosize;
    /* assert(len_min == len_max); */

    min = (double *) malloc(sizeof(double) * len_min);
    max = (double *) malloc(sizeof(double) * len_max);

    if (min == NULL || max == NULL)
        caml_failwith("ml_de_init_uniform: malloc failed on min/max c arrays");

    for (i = 0; i < len_min; i++)
        min[i] = Double_field(ml_min, i);
    for (i = 0; i < len_max; i++)
        max[i] = Double_field(ml_max, i);

    de_Init_Uniform(&prob_data, ndim, npop, min, max, comparison);

    free(min); free(max);
    
    CAMLreturn(Val_unit);    
}

// void de_Init_Solution(de_t *st, int dim, int popSize, double *min, double *max, int comparison,
//        double *aprioriSolution, double *popGaussianSigmas)
CAMLprim value ml_de_init_solution_native(value ml_ndim, value ml_npop, value ml_min, value ml_max, value ml_comparison,
    value ml_apriori_sol, value ml_pop_sigmas) {

    CAMLparam5(ml_ndim, ml_npop, ml_min, ml_max, ml_comparison);
    CAMLxparam2(ml_apriori_sol, ml_pop_sigmas);
    int ndim, npop, comparison, len_min, len_max, i, len_aprio, len_sigmas;
    double *min, *max, *apriori_sol, *pop_sigmas;

    ndim = Int_val(ml_ndim);
    npop = Int_val(ml_npop);
    comparison = Int_val(ml_comparison);

    //init_ml_trial(ndim);

    len_min = Wosize_val(ml_min) / Double_wosize;
    len_max = Wosize_val(ml_max) / Double_wosize;

    len_aprio = Wosize_val(ml_apriori_sol) / Double_wosize;
    len_sigmas = Wosize_val(ml_pop_sigmas) / Double_wosize;
    if (len_aprio != npop || len_sigmas != npop)
        caml_invalid_argument("ml_de_init_solution_native: Arrays length inconsistency on a priori solution, sigmas rel. to nPop");

    min = (double *) malloc(sizeof(double) * len_min);
    max = (double *) malloc(sizeof(double) * len_max);
    apriori_sol = (double *) malloc(sizeof(double) * npop);
    pop_sigmas = (double *) malloc(sizeof(double) * npop);

    if (min == NULL || max == NULL)
        caml_failwith("ml_de_init_solution_native: malloc failed on min/max c arrays");

    for (i = 0; i < len_min; i++)
        min[i] = Double_field(ml_min, i);
    for (i = 0; i < len_max; i++)
        max[i] = Double_field(ml_max, i);
    for (i = 0; i < npop; i++)
        apriori_sol[i] = Double_field(ml_apriori_sol, i);
    for (i = 0; i < npop; i++)
        pop_sigmas[i] = Double_field(ml_pop_sigmas, i);

    de_Init_Solution(&prob_data, ndim, npop, min, max, comparison, apriori_sol, pop_sigmas);

    free(min); free(max); free(apriori_sol); free(pop_sigmas);
    
    CAMLreturn(Val_unit);
}

CAMLprim value ml_de_init_solution_bytecode(value * argv, int argn) {
    return ml_de_init_solution_native(argv[0], argv[1], argv[2],
        argv[3], argv[4], argv[5], argv[6]);    
}

CAMLprim value ml_de_setup_std(value ml_scale, value ml_xover_prob, value ml_strategy) {
    CAMLparam3(ml_scale, ml_xover_prob, ml_strategy);
    double scale, xover_prob;
    void (*strategy_func)(de_t *st, int candidate);

    /* dumb compiler */
    strategy_func = Best1Exp;

    switch (Int_val(ml_strategy)) {
        case BEST1_EXP:
            strategy_func = Best1Exp;
            break;
        case RAND1_EXP:
            strategy_func = Rand1Exp;
            break;
        case BEST1_BIN:
            strategy_func = Best1Bin;
            break;
        case RAND1_BIN:
            strategy_func = Rand1Bin;
            break;
        case EITHER_OR:
            strategy_func = EitherOr;
            break;
        case BEST2_EXP:
            strategy_func = Best2Exp;
            break;
        case RAND2_EXP:
            strategy_func = Rand2Exp;
            break;
        case BEST2_BIN:
            strategy_func = Best2Bin;
            break;
        case RAND2_BIN:
            strategy_func = Rand2Bin;
            break;
    }

    scale = Double_val(ml_scale);
    xover_prob = Double_val(ml_xover_prob);

    de_SetupStd(&prob_data, scale, xover_prob, strategy_func, ml_energy_function);
    
    CAMLreturn(Val_unit);    
}


#ifdef RANDOMIZED_F_K

CAMLprim value ml_de_setup_random_fk(value ml_f_rand_fact, value ml_f_pop_fact,
    value ml_k_rand_fact, value ml_k_pop_fact) {
    CAMLparam4(ml_f_rand_fact, ml_f_pop_fact, ml_k_rand_fact, ml_k_pop_fact);
    double f_rand_fact, f_pop_fact, k_rand_fact, k_pop_fact;

    f_rand_fact = Double_val(ml_f_rand_fact);
    f_pop_fact = Double_val(ml_f_pop_fact);
    k_rand_fact = Double_val(ml_k_rand_fact);
    k_pop_fact = Double_val(ml_k_pop_fact);

    de_SetupRandomFK(&prob_data, f_rand_fact, f_pop_fact, k_rand_fact, k_pop_fact);
    
    CAMLreturn(Val_unit);    
}

#else

CAMLprim value ml_de_setup_random_fk(value ml_f_rand_fact, value ml_f_pop_fact,
    value ml_k_rand_fact, value ml_k_pop_fact) {
    CAMLparam4(ml_f_rand_fact, ml_f_pop_fact, ml_k_rand_fact, ml_k_pop_fact);
    
    CAMLreturn(Val_unit);    
}

#endif


CAMLprim value ml_de_setup_pcx_native(value ml_pcx_prob, value ml_scale, value ml_xover_prob,
    value ml_strategy, value ml_sigma2a, value ml_sigma2b) {
    CAMLparam5(ml_pcx_prob, ml_scale, ml_xover_prob, ml_strategy, ml_sigma2a);
    CAMLxparam1(ml_sigma2b);

    double pcx_prob, scale, xover_prob, sigma2a, sigma2b;
    void (*strategy_func)(de_t *st, int candidate);

    /* dumb compiler */
    strategy_func = Best1Exp;

    switch (Int_val(ml_strategy)) {
        case BEST1_EXP:
            strategy_func = Best1Exp;
            break;
        case RAND1_EXP:
            strategy_func = Rand1Exp;
            break;
        case BEST1_BIN:
            strategy_func = Best1Bin;
            break;
        case RAND1_BIN:
            strategy_func = Rand1Bin;
            break;
        case EITHER_OR:
            strategy_func = EitherOr;
            break;
        case BEST2_EXP:
            strategy_func = Best2Exp;
            break;
        case RAND2_EXP:
            strategy_func = Rand2Exp;
            break;
        case BEST2_BIN:
            strategy_func = Best2Bin;
            break;
        case RAND2_BIN:
            strategy_func = Rand2Bin;
            break;
    }

    pcx_prob = Double_val(ml_pcx_prob);
    scale = Double_val(ml_scale);
    xover_prob = Double_val(ml_xover_prob);
    sigma2a = Double_val(ml_sigma2a);
    sigma2b = Double_val(ml_sigma2b);

    de_SetupPCX(&prob_data, pcx_prob, scale, xover_prob,
        strategy_func, ml_energy_function, sigma2a, sigma2b);
    
    CAMLreturn(Val_unit);    
}

CAMLprim value ml_de_setup_pcx_bytecode(value * argv, int argn) {
    return ml_de_setup_pcx_native(argv[0], argv[1], argv[2],
        argv[3], argv[4], argv[5]);    
}

CAMLprim value ml_de_solve(value ml_max_iter) {
    CAMLparam1(ml_max_iter);

    int max_iter = Int_val(ml_max_iter);

    de_Solve(&prob_data, max_iter);

    CAMLreturn(Val_unit);
}

CAMLprim value ml_de_get_solution(value unit) {
    CAMLparam1(unit);
    CAMLlocal1(ml_bestsol);
    int i;

    ml_bestsol = caml_alloc(prob_data.nDim * Double_wosize, Double_array_tag);

    for (i = 0; i < prob_data.nDim; ++i)
        Store_double_field(ml_bestsol, i, prob_data.bestSolution[i]);

    CAMLreturn(ml_bestsol);
}

CAMLprim value ml_de_finalize(value unit) {
    CAMLparam1(unit);

    //caml_remove_global_root(&ml_trial);

    de_Finalize(&prob_data);

    CAMLreturn(Val_unit);   
}
