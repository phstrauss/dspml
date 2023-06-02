(*  one-dimensional discrete-time signal processing, with audio in mind.

    Frequency domain processing

    © Philippe Strauss, 2009 - 2012, 2015, 2021
    License: see accompanying LICENSE.txt or LICENSE.md file.

    References :

        [OppSch]    : Discrete-Time Signal Processing, Oppenheim & Schafer, Prentice Hall, 3rd ed.
        [ProMan]    : Digital Signal Processing, Proakis & Manolakis, Macmillan, 2nd ed.
        [OppWill]   : Signals & Systems, Oppenheim & Willsky, Prentice Hall, Intl 2nd ed.
        [ElFilt]    : Electronic Filter Design Handbook, 4th Ed.
        [IngPro]    : Digital Signal Processing using MATLAB, Ingle & Proakis  *)


open Oneliners
module FFT = Fftw3.D
open Myfftw
module Log = Logs


class convolve_fft_t ?scale_fwd ?scale_back coeft len_signal =
    let len1 = arlen coeft in
    let len_conv = len1+len_signal-1 in
    let _ = [%log debug "convolve_fft_t: lengths: coeft=%d, signalt=%d, overall=%d" len1 len_signal len_conv] in
    let coeftp = Array.append coeft (Array.make (len_conv-len1) 0.) in
    let ftd = new fft_complex_t len_conv FFT.Forward in
    let ftb = new fft_complex_t len_conv FFT.Backward in
    let coef_f = Array.make len_conv Complex.zero in
    let signalf = Array.make len_conv Complex.zero in
    let conv_time = Array.make len_conv Complex.zero in
    let sig_padding = Array.make (len_conv-len_signal) 0. in
    object(self)
        method exec signalt =
            assert (arlen signalt = len_signal) ;
            let signaltp = Array.append signalt sig_padding in
            begin match scale_fwd with
            | Some sc -> ( ftd #fillin_r2c coeftp ; ftd #exec ~scale:sc coef_f ; ftd #fillin_r2c signaltp ; ftd #exec ~scale:sc signalf )
            | None -> ( ftd #fillin_r2c coeftp ; ftd #exec coef_f ; ftd #fillin_r2c signaltp ; ftd #exec signalf ) end ;   
            let prod = Array.mapi (fun i x -> Complex.mul x signalf.(i)) coef_f in
            begin match scale_back with
            | Some sc -> ( ftb #fillin_c2c prod ; ftb #exec ~scale:sc conv_time )
            | None -> ( ftb #fillin_c2c prod ; ftb #exec conv_time ) end ;
            (prod, conv_time)
    end
