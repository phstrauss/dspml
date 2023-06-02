type quad_t = {
  xth : float;
  xc : float;
  yth : float;
  yc : float;
  a : float;
}
val quad_a : float -> float -> float * float
exception AxisTranslationBeyoundBounds of float
type headroom_t = Headroom_x of float | Headroom_y of float
val clip_level : headroom_t -> float
val quad_params : float -> headroom_t -> quad_t
module QuadUpper :
  sig
    val create : float -> headroom_t -> quad_t
    val axis_trans_x_hi : float -> quad_t -> float
    val axis_trans_y_hi : float -> quad_t -> float
    val square_clip_high : float -> quad_t -> float
  end
module QuadLower :
  sig
    val create : float -> headroom_t -> quad_t
    val axis_trans_x_lo : float -> quad_t -> float
    val axis_trans_y_lo : float -> quad_t -> float
    val square_clip_low : float -> quad_t -> float
  end
module QuadBoth :
  sig
    type t = quad_t * quad_t
    val create :
      float -> headroom_t -> float -> headroom_t -> quad_t * quad_t
    val square_clip : float -> quad_t * quad_t -> float
  end
module QuadSymetric :
  sig
    type t = quad_t * quad_t
    val sym_headroom : headroom_t -> headroom_t
    val create : float -> headroom_t -> quad_t * quad_t
    val square_clip : float -> quad_t * quad_t -> float
  end
module Atan : sig val soft_clip : float -> float -> float -> float end
