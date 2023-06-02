(*  hsv.ml
    Hue-Saturation-Value color computation - tests
    (c) Philippe Strauss, philippe@strauss-acoustics.ch  *)

open OUnit
open Hsvpalette

external show : 'a -> string = "%show"

let eps = 1e-12

let cmp_3f x y =
	match x with x1, x2, x3 ->
		match y with y1, y2, y3 ->
			cmp_float ~epsilon:eps x1 y1 && cmp_float ~epsilon:eps x2 y2 && cmp_float 
				~epsilon:eps x3 y3  

let oassert_3f = assert_equal ~cmp:cmp_3f

let () =
	oassert_3f (1., 0., 0.) (hue2rgb 0.) ; (* red *)
	oassert_3f (0., 1., 0.) (hue2rgb 120.) ; (* green *)
	oassert_3f (0., 0., 1.) (hue2rgb 240.) ; (* blue *)
	oassert_3f (1., 1., 0.) (hue2rgb 60.) ; (* yellow *)
	oassert_3f (0., 1., 1.) (hue2rgb 180.) ; (* cyan *)
	oassert_3f (1., 0., 1.) (hue2rgb 300.) ; (* magenta *) 
