/*  Jake Janovetz C implementation of remez exchange (or Parks-McClellan) algorithm for
    symetrical FIR filter design
    Objective Caml binding - C stubs file
    (c) Philippe Strauss, philippe@strauss-acoustics.ch, January 2012  */


#include <stdio.h>
#include <caml/mlvalues.h>
#include <caml/alloc.h>
#include <caml/memory.h>
#include <caml/callback.h>
#include <caml/fail.h>

#include <dspc/remez.h>


// both next : pretty neat, but always triggers a sigsegv (11)
void raise_malloc(const char *errmsg) {
    static value *exn_Malloc = NULL;

    if (exn_Malloc == NULL) 
        exn_Malloc = caml_named_value("exn_Malloc");

    caml_raise_with_string(*exn_Malloc, errmsg);
}

void raise_convergence(void) {
    static value *exn_Convergence = NULL;

    if (exn_Convergence == NULL) 
        exn_Convergence = caml_named_value("exn_Convergence");

    caml_raise_constant(*exn_Convergence);
}


/*
int remez_from_grid(
   double h[],
   int numtaps,
   int nextrems,
   double Grid[],
   double D[],
   double W[],
   int gridlen,
   int type,
   int symmetry)
*/

CAMLprim value ml_remez_from_grid(
    value ml_numtaps,
    value ml_Grid, value ml_D, value ml_W,
    value ml_gridlen, value ml_ftype,
    value ml_epsilon)
{
    int i, ret;
    int numtaps, gridlen, ftype;
    double *h, *Grid, *D, *W, epsilon;

    CAMLparam5(ml_numtaps, ml_Grid, ml_D, ml_W, ml_gridlen);
    CAMLxparam2(ml_ftype, ml_epsilon);
    CAMLlocal2(retup, ml_h);

    numtaps = Int_val(ml_numtaps);
    gridlen = Int_val(ml_gridlen);
    ftype = Int_val(ml_ftype);
    epsilon = Double_val(ml_epsilon);

    Grid = (double *) malloc(sizeof(double) * gridlen);
    D = (double *) malloc(sizeof(double) * gridlen);
    W = (double *) malloc(sizeof(double) * gridlen);
    h = (double *) malloc(sizeof(double) * numtaps);
    if(Grid == NULL || D == NULL || W == NULL || h == NULL)
        caml_failwith("ml_remez_from_grid : malloc failure on Grid/D/W/h c arrays");

    ml_h = caml_alloc(numtaps * Double_wosize, Double_array_tag);
    retup = caml_alloc(2, 0);

    for (i = 0; i < gridlen; i++) {
        Grid[i] = Double_field(ml_Grid, i);
        D[i] = Double_field(ml_D, i);
        W[i] = Double_field(ml_W, i);
    }

    ret = remez_from_grid(
        h,
        numtaps,
        Grid, D, W,
        gridlen, ftype,
        epsilon
    );

    for (i = 0; i < numtaps; i++)
        Store_double_field(ml_h, i, h[i]);


    Store_field(retup, 0, Val_int(ret));
    Store_field(retup, 1, ml_h);

    free(Grid);
    free(D);
    free(W);

    if (ret == RET_FAILED)
        caml_failwith("ml_remez_from_grid : Convergence problems");

    CAMLreturn(retup);
}

CAMLprim value ml_remez_from_grid_bytecode(value *argv, int argn) {
    return ml_remez_from_grid(argv[0], argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);
}

/*
int remez_from_bands(
    double h[],
    int numtaps,
    int numband,
    double transitions[],
    double des[],
    double weight[],
    int type)
*/

CAMLprim value ml_remez_from_bands(
    value ml_numtaps, value ml_numbands,
    value ml_transitions, value ml_des, value ml_weight,
    value ml_ftype, value ml_epsilon)
{
    int i, ntr, ndes, nw, ret;
    int numtaps, numbands, ftype;
    double *h, *transitions, *des, *weight, epsilon;

    CAMLparam5(ml_numtaps, ml_numbands, ml_transitions, ml_des, ml_weight);
    CAMLxparam2(ml_ftype, ml_epsilon);
    CAMLlocal2(retup, ml_h);

    numtaps = Int_val(ml_numtaps);
    numbands = Int_val(ml_numbands);
    ftype = Int_val(ml_ftype);
    epsilon = Double_val(ml_epsilon);

    ntr = Wosize_val(ml_transitions) / Double_wosize;
    ndes = Wosize_val(ml_des) / Double_wosize;
    nw = Wosize_val(ml_weight) / Double_wosize;
    if (nw != ndes && ndes*2 != ntr)
        caml_invalid_argument("ml_remez_from_bands : Arrays length inconsistency : desired/transitions/weighting");

    transitions = (double *) malloc(sizeof(double) * ntr);
    des = (double *) malloc(sizeof(double) * ndes);
    weight = (double *) malloc(sizeof(double) * nw);
    h = (double *) malloc(sizeof(double) * numtaps);
    if(transitions == NULL || des == NULL || weight == NULL || h == NULL)
        caml_failwith("ml_remez_from_bands : malloc failed on transitions/des/weight/h c arrays");

    ml_h = caml_alloc(numtaps * Double_wosize, Double_array_tag);
    retup = caml_alloc(2, 0);

    for (i = 0; i < ntr; i++)
        transitions[i] = Double_field(ml_transitions, i);
    for (i = 0; i < ndes; i++) {
        des[i] = Double_field(ml_des, i);
        weight[i] = Double_field(ml_weight, i);
    }

    ret = remez_from_bands(
        h,
        numtaps, numbands,
        transitions, des, weight,
        ftype, epsilon
    );

    for (i = 0; i < numtaps; i++)
        Store_double_field(ml_h, i, h[i]);

    Store_field(retup, 0, Val_int(ret));
    Store_field(retup, 1, ml_h);

    free(transitions);
    free(des);
    free(weight);

    if (ret == RET_FAILED)
        caml_failwith("ml_remez_from_bands : Convergence problems");

    CAMLreturn(retup);
}

CAMLprim value ml_remez_from_bands_bytecode(value *argv, int argn) {
    return ml_remez_from_bands(argv[0], argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);
}
