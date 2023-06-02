open Oneliners
open S_poles


let () =
	let order = 5 in
	pr "Bessel roots of order %d :\n" order ;
	let broots = bessel_poles order in
(* 	let w_scale = bessel_wc_scaling order in
	pr "w_scale length: %d\n" (arlen w_scale) ; *)
	for i = 0 to arlen broots - 1 do
		(* pr "w_scaling: %s\n" (spr_complex w_scale.(i)) ; *)
		pr "pole: %s\n" (spr_complex broots.(i))
	done
