import sys
import math


class ColorXform(object):
    def __init__(self):

        ''' compute matrix transform to convert between colorspaces.
        Note the output may be negative or out of range! '''

        # these are sRGB, D65 transforms from
        # http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
        # self.M_XYZ_to_RGB = [[3.24062548, -1.53720797, -0.49862860],
        #                      [-0.96893071, 1.87575606, 0.04151752],
        #                      [0.05571012, -0.20402105, 1.05699594]]

        # self.M_RGB_to_XYZ = [[0.41240000, 0.35760000, 0.18050000],
        #                      [0.21260000, 0.71520000, 0.07220000],
        #                      [0.01930000, 0.11920000, 0.95050000]]

        # these are 'Best RGB', D50 transforms from
        # http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
        self.space = "Best RGB"
        self.M_RGB_to_XYZ = [[0.6326696, 0.2045558, 0.1269946],
                             [0.2284569, 0.7373523, 0.0341908],
                             [0.0000000, 0.0095142, 0.8156958]]
                                

        self.M_XYZ_to_RGB = [[1.7552599, -0.4836786, -0.2530000],
                             [-0.5441336, 1.5068789, 0.0215528],
                             [0.0063467, -0.0175761, 1.2256959]]
                                
        
        # xyY coordinates for 590nm orange LED
        # https://www.luxalight.eu/en/cie-convertor
        self.orange_xyY = [0.575151311, 0.424232235, 1.0]
        # precompute XYZ color coordinates of orange LED
        self.orange_XYZ = self.xyY_to_XYZ(self.orange_xyY)
        #print(f"orange_XYZ: {self.orange_XYZ}")
        
        
    def matrix_mult(self, M, inp):
        # matrix multiply input triple by M
        # obv easier to do wityh a klibrary like numpy but doing
        # it long form here for clarity and library independence
        a =  M[0][0] * inp[0] + M[0][1] * inp[1] + M[0][2] * inp[2] 
        b =  M[1][0] * inp[0] + M[1][1] * inp[1] + M[1][2] * inp[2] 
        c =  M[2][0] * inp[0] + M[2][1] * inp[1] + M[2][2] * inp[2] 
        return(a, b, c)

    
    def hsv_to_rgba(self, h, s, v):
        """Improved implementation of hsv colorspace to red, green, blue and
        amber channels red and amber squashed into 1/3 the hue range
        (instead of 1/2 as in naive) hsv inputs and rgby outputs all
        floats between 0 and 1
        this code is described as in http://rotormind.com/blog/2015/Generating-RGBY-from-Hue/
        """
        # offset h so 0 is pure red, needed to keep code pretty
        h = h - 1.0/12.0
        if h < 0:
            h += 1.0

        if s == 0.0:
            return [v, v, v, v]

        i = int(h*6.0) # what hue range are we in?

                                # v is top flat
        f = (h*6.0) - i         # slope for 1/6 hue range

        b = v*(1.0 - s)         # bottom flat
        d = v*(1.0 - s*f)       # downslope  
        u = v*(1.0 - s*(1.0-f)) # upslope

        i2 = int(h*12.0)        # what hue subrange are we in?
        f2 = (h*12.0) - i2      # slope for 1/12 hue range
        d2 = v*(1.0 - s*f2)       # steep downslope  
        u2 = v*(1.0 - s*(1.0-f2)) # steep upslope

        i2 = i2 % 12

        if i2 == 0:
            return [d2, b, b, v]  # max a, r down steep
        if i2 == 1: 
            return [b, u2, b, v]  # max a, g up steep
        if i2 == 2 or i2 == 3:
            return [b, v, b, d]   # max g, a down slow
        if i2 == 4 or i2 == 5:
            return [b, v, u, b]   # max g, b up slow
        if i2 == 6 or i2 == 7:
            return [b, d, v, b]   # max b, g down slow
        if i2 == 8 or i2 == 9:
            return [u, b, v, b]   # max b, r up slow
        if i2 == 10:
            return [v, b, d2, b]  # max r, b down steep 
        if i2 == 11:
            return [v, b, b, u2] # max r, a up steep

    def XYZ_to_xyY(self, X, Y, Z):
        # from http://www.brucelindbloom.com/index.html?Eqn_xyY_to_XYZ.html
        if (X + Y + Z) == 0.0:
            print("warning, XYZ out of range")
            return(0.333, 0.333, 1)

        x= X/(X+Y+Z)
        y= Y/(X+Y+Z)
        return(x, y, Y)

    def xyY_to_XYZ(self, xyY):
        # from http://www.brucelindbloom.com/index.html?Eqn_xyY_to_XYZ.html

        x = xyY[0]
        y = xyY[1]
        Y = xyY[2]
        
        if y == 0:
            return(0., 0., 0.)

        return( x*Y/y, Y, (1 - x - y)*Y/y)


    def hue_angle(self, xy_vec):
        import math
        ''' calculate hue angle in degrees from CIE space vector'''
        # numpy version
        #angle = (360/(2*math.pi))*math.arctan2(xy_vec[0], xy_vec[1])
        angle = (1/(2*math.pi))*math.atan2(xy_vec[0], xy_vec[1])
        #sat = np.linalg.norm(xy_vec)
        return((angle + 1.)%1.)

    
    def RGBA_to_XYZ(self, R, G, B, A):
        red_XYZ= self.rgb_to_xyz(R,  0.,  0.)
        grn_XYZ= self.rgb_to_xyz(0.,  G,  0.)
        blu_XYZ= self.rgb_to_xyz(0., 0.,  B )

        org590_xy = [0.575151311, 0.424232235]
        org_XYZ = self.xyY_to_XYZ([org590_xy[0], org590_xy[1], A])


        X = red_XYZ[0] + grn_XYZ[0] + blu_XYZ[0] + org_XYZ[0]
        Y = red_XYZ[1] + grn_XYZ[1] + blu_XYZ[1] + org_XYZ[1]
        Z = red_XYZ[2] + grn_XYZ[2] + blu_XYZ[2] + org_XYZ[2]

        return(X, Y, Z)

    def XYZ_to_hue(self, X, Y, Z):
        import math
        ''' convert color/luminance specification in CIE XYZ coordinates to hue angle
            This version uses the colour_science library https://colour.readthedocs.io/'''
        # these angles are in xy space


        #D50 colorpoint from Best RGB colorspace
        white_xy = [ 0.3457,  0.3585]

        #red< 283.52 org< 254.014 grn< 162.59 blu< 33.67
        #print(f"red< {red_angle} org< {org_angle} grn< {grn_angle} blu< {blu_angle}")

        x, y, Y = self.XYZ_to_xyY(X, Y, Z)
        tmp_vec = [white_xy[0] - x, white_xy[1] - y]
        # calculate angle and offset so red is = zero = 360
        red_zero_angle = 5/36.
        return((2.0 - red_zero_angle - self.hue_angle(tmp_vec))%1.)

    def RGBA_to_HSV(self, R, G, B, A):
        maxc = max((R, G, B, A))
        minc = min((R, G, B, A))

        X, Y, Z = self.RGBA_to_XYZ(R,G,B,A)
        hue = self.XYZ_to_hue(X, Y, Z)

        # yes is is that simple
        val = maxc
        sat = 1. - minc
        return(hue, sat, val)

    def rgba_to_hsv(self, r, g, b, a):
        """ convert red, green, blue, amber color components to
        hue, saturation, value triple.
        This works by first converting to XYZ colorspace, interpolating
        in XYZ colorspace, then converting back into RGB space,
        then coverting to HSV
        all inputs/outputs are floats between 0 and 1"""
        
        pass

    
    def rgb_to_xyz(self, r, g, b, gamma=True):
        gval = 2.2
        if gamma:
            r = r**gval
            g = g**gval
            b = b**gval
#        print(r, g, b)
        X, Y, Z = self.matrix_mult(self.M_XYZ_to_RGB,
                                  [r, g, b])
        return(X, Y, Z)
        
    def xyz_to_rgb(self, x, y, z, gamma=True):
        gval = 1/2.2
        R, G, B = self.matrix_mult(self.M_RGB_to_XYZ,[x, y, z])

        if gamma:
            R = R**gval
            G = G**gval
            B = B**gval

        return(R, G, B)


    
    def rgba_to_rgbOLD(self, r, g, b, a):
        """ convert red, green, blue, amber color components to
        hue, saturation, value triple.
        This works by first converting to XYZ colorspace, interpolating
        in XYZ colorspace, then converting back into RGB space,
        then coverting to HSV"""

        rgb_XYZ = self.matrix_mult(self.M_RGB_to_XYZ, [r, g, b])
        
        #print(f"rgb_XYZ: {rgb_XYZ}")
        #print(f"org_XYZ: {self.orange_XYZ}")
 

        
        X =  rgb_XYZ[0] + a*self.orange_XYZ[0]
        Y =  rgb_XYZ[1] + a*self.orange_XYZ[1]
        # last term of this is ~ .001 so can be skipped for efficiency
        Z =  rgb_XYZ[2] + a*self.orange_XYZ[2]
        

        # now convert total XYZ back to rgb
        rgba_RGB = self.matrix_mult(self.M_XYZ_to_RGB,
                                  [X, Y, Z])


        maxc = max(rgba_RGB)
        #print(f"maxc is {maxc}")
        R = rgba_RGB[0]
        G = rgba_RGB[1]
        B = rgba_RGB[2]
            
        # this may now be greater than range so normalize
        if maxc > 1.:
            R = R/maxc
            G = G/maxc
            B = B/maxc
            

        # clamp negative values from underflow    
        if R < 0:
            R = 0.
        if G < 0:
            G = 0.
        if B < 0:
            B = 0.
        
        return(R, G, B)

    def rgb_to_hsv(self, r, g, b):
        # from https://www.niwa.nu/2013/05/math-behind-colorspace-conversions-rgb-hsl/
        #https://stackoverflow.com/questions/359612/how-to-convert-rgb-color-to-hsv/1626175#1626175
        maxc = max((r, g, b))
        minc = min((r, g, b))

        # calculate value (brightness)
        val = 0.5*(maxc + minc)
        # calculate saturation
        if val > 0.5:
            sat = ( maxc-minc)/(2.0-maxc-minc)
        else:
            if (maxc + minc) > 0:
                sat = (maxc-minc)/(maxc+minc)
            else:
                sat = 0
                
        # if maxc > 0:
        #     ##sat = (maxc - minc)/maxc
        #     sat = 1. - (minc/maxc)
        # else:
        #     sat = 0.

        # val = maxc
           

        # calculate hue
        if abs(maxc - minc) < 0.001:
            # color has close to zero saturation, hue is arbitrary, return 0
            hue = 0.
        
        elif (r > b) and (r > g): # if r is max
            print("red")
            hue = (g - b)/(maxc-minc)
        elif ( g > b):         # if g is max
            print("green")
            hue = 2.0 + (b - r)/(maxc - minc)
        else: # blue must be max
            print("blue")
            hue = 4.0 + (r - g)/(maxc - minc)

        # offset h so 0 is pure red, needed to keep code pretty
        hue = hue + 0.5
            
        # convert to range 0-1
        hue =  ((hue/6.) + 1.0)%1.0
        return (hue, sat, val)



    
    def print_hsv(self, hue, sat, val, extra=""):
        print(f" hsv: {hue:.3f} {sat:.3f} {val:.3f} {extra}")            

    def print_rgba(self, r, g, b, a, extra=""):
        print(f"rgba: {r:.3f} {g:.3f} {b:.3f} {a:.3f} {extra}")        

        
if __name__ == '__main__':

    cx = ColorXform()

    steps = 16

    maxhdiff = 0.0
    for hue_degrees in range(0,360,10):
        hue = float(hue_degrees)/360.
        sat = 0.7
        val = 1.0
        print(f" hsv: {hue:.3f} {sat:.3f} {val:.3f}")        
        r, g, b, a = cx.hsv_to_rgba(hue, sat, val)
        #print(f"rgba: {r:.3f} {g:.3f} {b:.3f} {a:.3f}")        
        R, G, B = cx.rgba_to_rgb(r,g,b,a)
        #print(f"RGB: {R:.3f} {G:.3f} {B:.3f}")
        H, S, V = cx.rgb_to_hsv(R, G, B)
        print(f" HSV: {H:.3f} {S:.3f} {V:.3f}")
        huediff = abs(hue-H)
        if huediff > 0.5:
            huediff = abs(1-huediff)
        
        if huediff > maxhdiff:
            maxhdiff = huediff
        print(f"hue diff: {huediff}")

        print("===============")
        print(f"max hue diff: {maxhdiff}")

        
        
