"""Generate plots of color trajectories from """


# coulour is color_science
# https://colour.readthedocs.io/en/develop/#

import numpy as np
from pprint import pprint
import matplotlib.pyplot as plt

import colour
from colour.plotting import (
    plot_RGB_chromaticities_in_chromaticity_diagram_CIE1931,
    colour_style,
    override_style,
    artist
)
# https://colour.readthedocs.io/en/develop/generated/colour.plotting.plot_RGB_colourspaces_in_chromaticity_diagram_CIE1931.html


#from colour.plotting.common import KwargsArtist


colour_style()
# make figures a little smaller for blog posts
override_style(**{ "figure.figsize": (6., 6.)})
plt.rcParams['savefig.facecolor'] = "0.8"
plt.tight_layout()

# turn off interactive
#plt.ioff()
plt.ion()

if True:

    # xyY coordinates for 590nm orange LED
    # https://www.luxalight.eu/en/cie-convertor

    org_xyY = np.array([0.575151311, 0.424232235, 1.0])

    #print(org_xyY) 

    org_XYZ = colour.xyY_to_XYZ(org_xyY) 
    org_xyY = colour.XYZ_to_xyY(org_XYZ) 
    
    print("orange XYZ:", org_XYZ)

    colourspace = colour.RGB_COLOURSPACES['sRGB']



    chromatic_adaptation_transform = 'Bradford'

  # these are sRGB, D65 transforms from
    # http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
    matrix_XYZ_to_RGB = np.array(
        [[3.24062548, -1.53720797, -0.49862860],
         [-0.96893071, 1.87575606, 0.04151752],
         [0.05571012, -0.20402105, 1.05699594]]
    )


    matrix_RGB_to_XYZ = np.array(
        [[0.41240000, 0.35760000, 0.18050000],
         [0.21260000, 0.71520000, 0.07220000],
         [0.01930000, 0.11920000, 0.95050000]]
    )


    rgb_comp = np.array([1, 1, 1])
    comp_XYZ = colour.RGB_to_XYZ(rgb_comp,
                                 colourspace.whitepoint,
                                 colourspace.whitepoint,
                                 matrix_RGB_to_XYZ,
                                 chromatic_adaptation_transform)  

    red_XYZ = colour.RGB_to_XYZ([1., 0., 0.],
                                 colourspace.whitepoint,
                                 colourspace.whitepoint,
                                 matrix_RGB_to_XYZ,
                                 chromatic_adaptation_transform)  
    grn_XYZ = colour.RGB_to_XYZ([0., 1., 0.],
                                 colourspace.whitepoint,
                                 colourspace.whitepoint,
                                 matrix_RGB_to_XYZ,
                                 chromatic_adaptation_transform)  
    blu_XYZ = colour.RGB_to_XYZ([0., 0., 1.],
                                 colourspace.whitepoint,
                                 colourspace.whitepoint,
                                 matrix_RGB_to_XYZ,
                                 chromatic_adaptation_transform)  

    wht_XYZ = colour.RGB_to_XYZ([1., 1., 1.],
                                 colourspace.whitepoint,
                                 colourspace.whitepoint,
                                 matrix_RGB_to_XYZ,
                                 chromatic_adaptation_transform)  

    
    mix_XYZ = (0.5*org_XYZ + 0.5*comp_XYZ)


    mix_RGB = colour.XYZ_to_RGB(mix_XYZ,
                                colourspace.whitepoint,
                                colourspace.whitepoint,
                                matrix_XYZ_to_RGB,
                                None)  

    org_RGB = colour.XYZ_to_RGB(org_XYZ,
                                colourspace.whitepoint,
                                colourspace.whitepoint,
                                matrix_XYZ_to_RGB,
                                chromatic_adaptation_transform)  

    #tot_HSV = colour.RGB_to_HSV(tot_RGB)


    # print("HSV is:")
    # print(tot_HSV)


    # print("tot_RGB")
    # print(tot_RGB)


    plot_rgbs = [org_RGB, rgb_comp]

    figure = plt.figure()
    
    axes = figure.add_subplot()
    #axes = figure.add_subplot(212)


    colourspace = colour.RGB_COLOURSPACES['sRGB']


    if True:

        RGB = []
        
        for c_XYZ in [red_XYZ, grn_XYZ, blu_XYZ,wht_XYZ]:
            RGB.append(colour.XYZ_to_RGB(
                c_XYZ,
                colourspace.whitepoint,
                colourspace.whitepoint,
                colourspace.matrix_XYZ_to_RGB)
            )
        
        plot_RGB_chromaticities_in_chromaticity_diagram_CIE1931(
            #RGB = RGB,
            RGB = [org_RGB, mix_RGB],
            colourspace = "sRGB",
            #        None,
            axes = axes,
            #        colourspaces=["ACEScg"],
            colourspaces=[],
            scatter_kwargs={'c': 'k', 'marker': 'x'}

        )


        red_xyY = colour.XYZ_to_xyY(red_XYZ) 
        grn_xyY = colour.XYZ_to_xyY(grn_XYZ) 
        blu_xyY = colour.XYZ_to_xyY(blu_XYZ) 
        org_xyY = colour.XYZ_to_xyY(org_XYZ) 
        wht_xyY = colour.XYZ_to_xyY(wht_XYZ) 


        # overplot gamut in black for visibility
        axes.plot([red_xyY[0], grn_xyY[0]],
                  [red_xyY[1], grn_xyY[1]], 'k' , linestyle='dotted')
        axes.plot([blu_xyY[0], grn_xyY[0]],
                  [blu_xyY[1], grn_xyY[1]], 'k' )
        axes.plot([blu_xyY[0], red_xyY[0]],
                  [blu_xyY[1], red_xyY[1]], 'k' )

        axes.plot([red_xyY[0], org_xyY[0]],
                  [red_xyY[1], org_xyY[1]], 'k' )

        axes.plot([org_xyY[0], grn_xyY[0]],
                  [org_xyY[1], grn_xyY[1]], 'k' )


        axes.plot([org_xyY[0], wht_xyY[0]],
                  [org_xyY[1], wht_xyY[1]], 'k', linestyle='dashed' )


        
        
        axes.text(org_xyY[0], org_xyY[1], "orange LED")
        figure.canvas.draw()
        figure.show()
        axes.redraw_in_frame()
        figure.savefig("myplot.png")
        
    
    print("done")

exit()
