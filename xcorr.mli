module FFT = Fftw3.D
module D = Lacaml.D
module Z = Lacaml.Z
val create :
  ('a, 'b) Bigarray.kind ->
  int -> ('a, 'b, Bigarray.fortran_layout) Bigarray.Array1.t
val log2 : float
val minpow2 : int -> int
val copy0 : ?n:int -> ?ofsx:int -> D.vec -> ?ofsy:int -> D.vec -> unit
type scale = Biased | Unbiased | Coeff
val xcorr :
  ?maxlag:int ->
  ?scale:scale ->
  D.vec ->
  D.vec -> (float, FFT.float_elt, Bigarray.fortran_layout) Bigarray.Array1.t
