module FFT = Fftw3.D
class fft_complex_t :
  int ->
  Fftw3.D.dir ->
  object
    method exec : ?scale:float -> Complex.t array -> unit
    method fillin_c2c : Complex.t array -> unit
    method fillin_r2c : float array -> unit
  end
val hc_len_complex : int -> int
val hc2c :
  float -> (float, 'a, 'b) Bigarray.Array1.t -> Complex.t array -> unit
val c2hc :
  float -> Complex.t array -> (float, 'a, 'b) Bigarray.Array1.t -> unit
class fft_real_t :
  int ->
  object
    method exec_back : ?scale:float -> Complex.t array -> float array -> unit
    method exec_fwd : ?scale:float -> float array -> Complex.t array -> unit
  end
