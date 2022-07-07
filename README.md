# colorxform
Color math for converting between color spaces for RGB + Amber LEDs

`plot_trajectories.py` --- code to plot LED colors and gamut in CIE space. This uses the `color_science` python library [] and generates the plot below.

`colorxform.py` --- code to  convert between RGB, RGBA (red, green blue, & amber), HSV, CIE XYZ and CIE xyY colorspaces with no library dependencies.

This code was written in an attempt to get good color specification for LED light fixtures using more than standard RGB colors. Though right now only amber is supported, these methods could be easily extended to LEDs or other sources of any color. This is a standalone module and does not depend on any libraries, and is written to be easily translated into other languages (hence unPythonic constructs). 

The `HSV_to_RGBA()` function divides the hue range into four sections corresponding to red, amber, green, and blue, and crossfades between them to arrive at an approximately correct hue. More details can be found at this blog post: []

The `RGBA_to_HSV()` function is somewhat more complex and is naturally an approximation as we are converting from a four-dimensional space to three. The approach is to convert each light source (R, G, B, and A) into the CIE XYZ space, where they can be interpolated to find the color of the resulting mixture. Once this is determined (by a simple linear mix), the XYZ coordinates can be converted back to sRGB (albeit with a loss of gamut) and then to HSV.

In this example, the dotted line shows the colors available by mixing the amber LED with white (from an equal mix of RGB).  The X is the point at R = G = B = A = 0.5. 


![CIE colorspace plot showing gamut of RGB LEDS plus amber LED](https://github.com/headrotor//blob/main/amber-trajectory.png?raw=true)
