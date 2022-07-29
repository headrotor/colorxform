# colorxform
### Color math for converting between color spaces for RGB + Amber LEDs

This code was written to control the colors of LED light
fixtures that use more than standard RGB colors. Though right now only
an additional amber is supported, these methods could be easily
extended to LEDs or other sources of any color. This is a standalone
module and does not depend on any libraries, and is written to be
easily translated into other languages (hence unPythonic constructs).

The `HSV_to_RGBA()` and `HSV_to_RGBA_CIE()` functions divides the hue
range into four sections corresponding to red, amber, green, and blue,
and crossfades between them to arrive at an approximately correct
hue. 

The `RGBA_to_HSV_CIE()` function is somewhat more complex and is naturally
an approximation as we are converting from a four-dimensional space to
three. The approach is to convert each light source (R, G, B, and A)
into the [CIE xy](https://en.wikipedia.org/wiki/CIE_1931_color_space#Definition_of_the_CIE_XYZ_color_space)
space, where they can be interpolated to find the color of the
resulting mixture. Once this is calculated (by a simple linear mix),
the `xy` coordinates can be converted to HSV.

In this example, the dotted line shows the colors available by mixing
the amber LED with white (from an equal mix of RGB).  The X is the
point at R =  G = B = A = 0.5.

![CIE colorspace plot showing gamut of RGB LEDS plus amber
LED](https://github.com/headrotor/colorxform/blob/main/amber-trajectory.png?raw=true)


There are two versions of the HSV to RGBA conversion, the first
`HSV_to_RGBA()` is described here. The second approach
`HSV_to_RGBA_CIE()` gives colors that are closer to equal angles in
`xy` space and should thus be a perceptually smoother transition
around the hue circle. The following image shows the locus in `xy` space of 
hues equally spaced in 10 degree increments given by `HSV_to_RGBA()`  on the left and `HSV_to_RGBA_CIE()` on the right.

![Equal hue values in naive HSV_to_RGB and perceptually](https://github.com/headrotor/colorxform/blob/main/images/CIE-gamuts.png?raw=true)

<!-- How to change images based on mode  -->
<!-- ![Fancy logo](./dark.png#gh-dark-mode-only) -->
<!-- ![Fancy logo](./light.png#gh-light-mode-only) -->

A benefit of the latter is that a round-trip conversion `HSV` →
`RGBA`→ `HSV` will result in hues that are much closer to the
original. The `HSV_to_RGBA_CIE()` method has an optional parameter
`warp` that when set to true (the default) will prewarp the hue to
better match the `RGBA_to_HSV_CIE()` results.  `warp` uses a 6th-order
polynomial fit to warp the input hue value; if round-trip hue matching
is of less importance then `warp` mat be set to False to skip this
computation.


## Contents of this repository:

Note on variable names<sup>[Note 1](#note1)</sup>

`colorxform.py` --- Python class with methods to convert between RGB, RGBA (red, green blue, & amber), HSV, CIE XYZ and CIE xyY colorspaces with no library
dependencies.  For more information on the CIE colorspaces, see Wikipedia: 
[https://en.wikipedia.org/wiki/CIE_1931_color_space}

There are three principal methods in this class:

* `HSV_to_RGBA(self, h, s, v)` -- This is the orginal HSV-to-RGBA
  conversion routine described here:
  <http://rotormind.com/blog/2015/Generating-RGBY-from-Hue/>
* `HSV_to_RGBA_CIE(self, h, s, v, warp=True)` This is an improved version with better perceptual linearity in the CIE space.
* `RGBA_to_HSV_CIE(self, R, G, B, A)` This is the inverse function
    that converts a color specified in linear red, green, blue and
    amber value, for example as PWM brightnesses for each LED of the
    respective colors, into hue, value (brightness) and saturation
    values. This calls three additional functions:
  1. `RGBA_to_xy(self, R, G, B, A)` -- Converts RGBA values to CIE _xy_
   values.
  2. `xy_to_hue(self, xy_vec, red_offset=None)`-- Converts a CIE _xy_ value to a hue angle by finding the angle from pure red in _xy_ space.
  3. `prewarp(self, x, coeff=None)` -- A hue warping function to help
   correct nonlinearity in the `RGBA_to_HSV_CIE()` function so that the
   round trip `RGBA_to_HSV_CIE(HSV_to_RGB())` conversion results in closer
   to identical hue values. This is called automatically when the
   default `warp=True` parameter is used in `HSV_to_RGBA_CIE()`; there
   is little reason not to use it.


`plot_trajectories.py` --- code to plot LED colors and gamut in CIE
space. This uses the `color_science` python library
<https://www.colour-science.org/> and generates the plot above

`linearized_RGBA_to_hue.ipynb` --- A Juypyter notebook to explore and test the colorxform methods. In particular this calculates the regression coefficients for the hue warping function used in `prewarp()`.

<a name="note1">Note 1</a>: color-space nomenclature is confusing!
Colorspace names are in uppercase except `xyY` to disambiguate `X`,
`Y` from `x`, `y` from CIE `xyY` space. Four-color space is `RGBA`
with `A` for Amber.  (Though `Y` is often used for for amber/Yellow to
disambiguate from alpha channel "A", we use `A` here to disambiguate
from `Y` in `xyY` and `XYZ` color spaces.)



