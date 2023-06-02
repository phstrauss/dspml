module FFT = Fftw3.D
class convolve_fft_t :
  ?scale_fwd:float ->
  ?scale_back:float ->
  float array ->
  int ->
  object method exec : float array -> Complex.t array * Complex.t array end
