(*  one-dimensional discrete-time signal processing, with audio in mind.

    7.3.4 Normal Deviates by Transformation (Box-Muller) (Numerical Recipes 3rd ed.)

    © Philippe Strauss, summer 2015
    License: see accompanying LICENSE.txt or LICENSE.md file.

    References :

        [OppSch]    : Discrete-Time Signal Processing, Oppenheim & Schafer, Prentice Hall, 3rd ed.
        [ProMan]    : Digital Signal Processing, Proakis & Manolakis, Macmillan, 2nd ed.
        [OppWill]   : Signals & Systems, Oppenheim & Willsky, Prentice Hall, Intl 2nd ed.
        [ElFilt]    : Electronic Filter Design Handbook, 4th Ed.
        [IngPro]    : Digital Signal Processing using MATLAB, Ingle & Proakis  *)


(*  INLINE double normal_double(struct rng_s *st, double mu, double sigma) {
        double v1, v2, rsq, fac;

        if (st->normal_stored == 0.0) {
            do {
                v1 = 2.0 * uniform_double(st) - 1.0;
                v2 = 2.0 * uniform_double(st) - 1.0;
                rsq = v1 * v1 + v2 * v2;
            } while (rsq >= 1.0 || rsq == 0.0);
            fac = sqrt(-2.0 * log(rsq) / rsq);
            st->normal_stored = v1 * fac;
            return mu + sigma * v2 * fac;
        } else {
            fac = st->normal_stored;
            st->normal_stored = 0.0;
            return mu + sigma * fac;
        }
    } *)


type t = { mutable normal_stored : float }

let create () =
    { normal_stored = 0.; }

let normal_double st mu sigma =
    if st.normal_stored = 0. then (
        let loop = ref true in
        let rsq = ref 0.
        and v1 = ref 0.
        and v2 = ref 0. in
        while !loop do
            v1 := 2. *. (Random.float 1.) -. 1. ;
            v2 := 2. *. (Random.float 1.) -. 1. ;
            rsq := !v1 ** 2. +. !v2 ** 2. ;
            if !rsq < 1. && !rsq <> 0. then loop := false
        done ;
        let fac = sqrt (-2. *. (log !rsq) /. !rsq) in
        st.normal_stored <- !v1 *. fac ;
        mu +. sigma *. !v2 *. fac
    ) else (
        let fac = st.normal_stored in
        st.normal_stored <- 0. ;
        mu +. sigma *. fac
    )