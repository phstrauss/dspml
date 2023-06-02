(*  one-dimensional discrete-time signal processing, with audio in mind.

    Sound file I/O using libsndfile.

    © Philippe Strauss, 2009 - 2012; 2020
    License: see accompanying LICENSE.txt or LICENSE.md file.

    References :

        [OppSch]    : Discrete-Time Signal Processing, Oppenheim & Schafer, Prentice Hall, 3rd ed.
        [ProMan]    : Digital Signal Processing, Proakis & Manolakis, Macmillan, 2nd ed.
        [OppWill]   : Signals & Systems, Oppenheim & Willsky, Prentice Hall, Intl 2nd ed.
        [ElFilt]    : Electronic Filter Design Handbook, 4th Ed.
        [IngPro]    : Digital Signal Processing using MATLAB, Ingle & Proakis  *)


open Oneliners
module Sf = Sndfile
module Log = Logs


(* basic interface, read / write in one go for small files *)

type samp_data_atonce_t = Mono of float array | Multi of float array array
exception FileTooLong of Int64.t
exception NrOfChannels of int


let deinterleave interleaved chan =
    let frames_tot = arlen interleaved in
    let frames_chan = frames_tot / chan
    and frames_over = frames_tot mod chan in
    assert (frames_over = 0) ;
    let demuxd = Array.make_matrix chan frames_chan 0.0 in
    for i = 0 to frames_tot-1 do
        demuxd.(i mod chan).(i / chan) <- interleaved.(i)
    done ;
    demuxd

let interleave uninterleaved =
    let chan = arlen uninterleaved in
    let nframes = arlen uninterleaved.(0) in
    if ((Int64.mul (Int64.of_int nframes) (Int64.of_int chan)) > (Int64.of_int Sys.max_array_length)) then
        raise (FileTooLong (Int64.of_int nframes)) ;
    let interleaved = Array.make (nframes * chan) 0.0 in
    for i = 0 to (nframes*chan)-1 do
        interleaved.(i) <- uninterleaved.(i mod chan).(i / chan)
    done ;
    interleaved

let read_atonce filename =
    let sfd = Sf.openfile filename in
    let chan = Sf.channels sfd
    and sfreq = Sf.samplerate sfd
    and nframes = Sf.frames sfd in
    [%log debug "read_wav: opened %s: %d ch; sfreq %d Hz; %Ld frames" filename chan sfreq nframes] ;
    if nframes > (Int64.div (Int64.of_int Sys.max_array_length) (Int64.of_int chan)) then
        raise (FileTooLong nframes) ;
    let iframes = (Int64.to_int nframes) in
    let interleaved = Array.make (iframes*chan) 0.0 in
    let readcount = Sf.read sfd interleaved in
    [%log debug "read_wav: %s: read %d items" filename readcount] ;
    match chan with
    | 1 -> ( Mono interleaved, sfreq )
    | _ -> ( Multi (deinterleave interleaved chan), sfreq)

let read_mono_atonce filename =
    match (read_atonce filename) with
    | (Mono data, sfreq) -> (data, sfreq)
    | (Multi data, _) -> raise (NrOfChannels (arlen data.(0)))

let write_atonce filename sfreq data =
    let fmt = Sf.format Sf.MAJOR_WAV Sf.MINOR_FLOAT in
    match data with
    | Mono samples -> begin
        let sfd = Sf.openfile ~info:(Sf.WRITE, fmt, 1, sfreq) filename in
        [%log debug "write_atonce: opened %s: mono; sfreq %d Hz" filename sfreq] ;
        let count = arlen samples in
        let writecount = Sf.write sfd samples count in
        [%log debug "write_atonce: %s: wrote %d items" filename writecount] ;
        Sf.close sfd
    end
    | Multi samples -> begin
        let chan = arlen samples in
        let interleaved = interleave samples in
        [%log debug "write_atonce: %d channels, %d frames %d interleaved" chan (arlen samples.(0)) (arlen interleaved)] ;
        let sfd = Sf.openfile ~info:(Sf.WRITE, fmt, chan, sfreq) filename in
        [%log debug "write_atonce: opened %s: %d channels; sfreq %d Hz" filename chan sfreq] ;
        let count = arlen interleaved in
        let writecount = Sf.write sfd interleaved count in
        [%log debug  "write_atonce: %s: wrote %d items" filename writecount] ;
        Sf.close sfd
    end ;


(* read linear PCM files block by block - tested with MS RIFF WAV files *)

type sf_info_t = {
    fd: Sndfile.t;
    nchan: int;
    sfreq: int;
    frames: Int64.t;
}

let open_sf filename =
    let sfd = Sf.openfile filename in
    let chan = Sf.channels sfd
    and sfreq = Sf.samplerate sfd
    and nframes = Sf.frames sfd in
    [%log debug "open_sf: opened '%s': %d ch; sfreq %d Hz; %Ld frames" filename chan sfreq nframes] ;
    {fd=sfd; nchan=chan; sfreq=sfreq; frames=nframes}

let read_block_init bsize
                    ?(demuxlen=bsize)
                    ?(overlap=0) (* for overlap-save *)
                    ?(pcm16_scale=false) (* normaly useless *)
                    finfo =
    let inlen = bsize * finfo.nchan in
    let interleaved = Array.make (bsize * finfo.nchan) 0.0 in
    let demuxd = Array.make_matrix finfo.nchan demuxlen 0.0 in
    let times = ref 0 in
    let reader () =
        let seek = if !times > 0 then !times * (* finfo.nchan * *) (bsize - overlap) else 0 in
        let offset = Sf.seek finfo.fd (Int64.of_int seek) Sf.SEEK_SET in
        let readcount = Sf.read finfo.fd interleaved in
        [%log debug "read_block reader: offset %d; read %d items" (Int64.to_int offset) readcount] ;
        incr times ;
        (* deinterleave *)
        if pcm16_scale (* libsndfile takes care of it usually *) then (
            for i = 0 to readcount-1 do
                demuxd.(i mod finfo.nchan).(i / finfo.nchan) <- interleaved.(i) /. 32768.
            done
        ) else (
            for i = 0 to inlen-1 do
                demuxd.(i mod finfo.nchan).(i / finfo.nchan) <- interleaved.(i)
            done
        ) ;
        (offset, readcount / finfo.nchan)
    in (reader, demuxd)

let sumframes offset readframes =
    (Int64.add offset (Int64.of_int readframes))

let write_block_init filename channels sfreq format bsize uninterleaved =
    (* let format = Sf.format Sf.MAJOR_WAV Sf.MINOR_FLOAT in *)
    let sfd = Sf.openfile ~info:(Sf.WRITE, format, channels, sfreq) filename in
    let interleaved = Array.make (bsize * channels) 0.0 in
    let writer ?(writelen=bsize) () =
        let count = writelen * channels in
        for i = 0 to count - 1 do
            interleaved.(i) <- uninterleaved.(i mod channels).(i / channels) 
        done ;
        let wrotecount = Sf.write sfd interleaved count in
        [%log debug "write_block writer: wrote %d items" wrotecount] ;
        assert (wrotecount = count)
    in (writer, sfd)
