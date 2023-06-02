module Sf :
  sig
    type open_mode_t = Sndfile.open_mode_t = READ | WRITE | RDWR
    type seek_mode_t = Sndfile.seek_mode_t = SEEK_SET | SEEK_CUR | SEEL_END
    type major_format_t =
      Sndfile.major_format_t =
        MAJOR_NONE
      | MAJOR_WAV
      | MAJOR_AIFF
      | MAJOR_AU
      | MAJOR_RAW
      | MAJOR_PAF
      | MAJOR_SVX
      | MAJOR_NIST
      | MAJOR_VOC
      | MAJOR_IRCAM
      | MAJOR_W64
      | MAJOR_MAT4
      | MAJOR_MAT5
      | MAJOR_PVF
      | MAJOR_XI
      | MAJOR_HTK
      | MAJOR_SDS
      | MAJOR_AVR
      | MAJOR_WAVEX
      | MAJOR_SD2
      | MAJOR_FLAC
      | MAJOR_CAF
    type minor_format_t =
      Sndfile.minor_format_t =
        MINOR_NONE
      | MINOR_PCM_S8
      | MINOR_PCM_16
      | MINOR_PCM_24
      | MINOR_PCM_32
      | MINOR_PCM_U8
      | MINOR_FLOAT
      | MINOR_DOUBLE
      | MINOR_ULAW
      | MINOR_ALAW
      | MINOR_IMA_ADPCM
      | MINOR_MS_ADPCM
      | MINOR_GSM610
      | MINOR_VOX_ADPCM
      | MINOR_G721_32
      | MINOR_G723_24
      | MINOR_G723_40
      | MINOR_DWVW_12
      | MINOR_DWVW_16
      | MINOR_DWVW_24
      | MINOR_DWVW_N
      | MINOR_DPCM_8
      | MINOR_DPCM_16
    type endianness_t =
      Sndfile.endianness_t =
        ENDIAN_FILE
      | ENDIAN_LITTLE
      | ENDIAN_BIG
      | ENDIAN_CPU
    type file_format_t = Sndfile.file_format_t
    type error =
      Sndfile.error =
        No_error
      | Unrecognised_format
      | System
      | Malformed_file
      | Unsupported_encoding
      | Internal
    exception Error of (error * string)
    type t = Sndfile.t
    val format : major_format_t -> minor_format_t -> file_format_t
    val format_e :
      major_format_t -> minor_format_t -> endianness_t -> file_format_t
    val openfile :
      ?info:open_mode_t * file_format_t * int * int -> string -> t
    val close : t -> unit
    val read : t -> float array -> int
    val write : t -> float array -> int -> int
    val frames : t -> Int64.t
    val samplerate : t -> int
    val channels : t -> int
    val seek : t -> Int64.t -> seek_mode_t -> Int64.t
    val compare : t -> t -> int
  end
type samp_data_atonce_t = Mono of float array | Multi of float array array
exception FileTooLong of Int64.t
exception NrOfChannels of int
val deinterleave : float array -> int -> float array array
val interleave : float array array -> float array
val read_atonce : string -> samp_data_atonce_t * int
val read_mono_atonce : string -> float array * int
val write_atonce : string -> int -> samp_data_atonce_t -> unit
type sf_info_t = {
  fd : Sndfile.t;
  nchan : int;
  sfreq : int;
  frames : Int64.t;
}
val open_sf : string -> sf_info_t
val read_block_init :
  int ->
  ?demuxlen:int ->
  ?overlap:int ->
  ?pcm16_scale:bool ->
  sf_info_t -> (unit -> Int64.t * int) * float array array
val sumframes : int64 -> int -> int64
val write_block_init :
  string ->
  int ->
  int ->
  Sf.file_format_t ->
  int -> float array array -> (?writelen:int -> unit -> unit) * Sf.t
