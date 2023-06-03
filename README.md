# Libdspml

OCaml Routines for One-dimension Discrete-time Signal Processing on commodity hardware.

- Analog S plane poles computation (s_poles.ml)
- Bilinear Transform (bilinear.ml)
- LTI transfer functions analysis (lti_xfer.ml)
- IIR & FIR filters simple inner processing loops (filter_process.ml)
- Windows (windows_cosine.ml)
- FFTW wrapper (myfftw.ml)
- Remez exchange FIR filter design binding
- Simple single-threaded Differential Evolution solver binding
- Box-Müller random transform (boxMuller.ml)
- Some helpers about finding minimas and maximas in a vector of float (minmax.ml)
- Nonlinear soft-clipping (second order softening curve, nonlin.ml)
- Simple soundfile helper on top of libsndfile (sfio.ml)
- Some search function of values into a float vector (vec_search.ml)
- Simple common quantitave function on vector of flots (vec_simple.ml)
- Cross correlation (taken from Christophe Troestler, xcorr.ml)
- Zero padding and trimming helpers (zeropad_trim.ml)
- Hue-Saturation processing (hsv.ml)
- Misc. see intderiv.ml and spectral.ml

(c) 2023 Philippe Strauss <catseyechandra at proton dot me>

License : LGPLv3, see [https://www.gnu.org/licenses/lgpl-3.0.html](https://www.gnu.org/licenses/lgpl-3.0.html)