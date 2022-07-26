import sys
import math


class ColorXform(object):
    def __init__(self):
        ''' ColorXform:

        This library is written for easy translation to other languages
        e.g. c++, JavaScript, so many idioms are not as pythonic as they
        could be. Performance and look could be improved if necessary by using
        NumPy in particular. 

        naming conventions: colorspace variable names are confusing!
        Colorspace names are in uppercase except xyY to 
        disambiguate X, Y from x, y
        4-color space is RGBA with A for Amber.
        (Even though some use Y for Yellow to disambiguate from 
        alpha channel "A", we  use A here to disambiguate from 
        Y in xyY and XYZ color spaces.) 
        '''


        # matrix transform to convert between colorspaces.
        #Note the output may be negative or out of range! '''
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

        # gamma value for RGBA scaling
        self.gamma = 2.2

        # CIE xy coordinates for LEDs assuiming monochromatic
        # https://www.luxalight.eu/en/cie-convertor
        # "neopixel RGB": assumed to be 630nm/530nm/475nm from  datasheets

        self.red630_xy = [0.707917792, 0.292027109]
        self.red640_xy = [0.719032942, 0.280934952]
        self.org590_xy = [0.575151311, 0.424232235]
        self.grn530_xy = [0.154722061, 0.805863545]
        self.blu475_xy = [0.109594324, 0.086842511]
        self.uvv405_xy = [0.173020965, 0.004775050]      

        
        # standard white illuminants in CIE xy space, 5000K, 65000K, sRGB
        self.D65_xy =   [0.34567, 0.35850]
        self.D50_xy =   [0.31271, 0.32902]
        self.sRGB1_xy = [0.4557, 0.4211]
        self.EE_xy =    [0.3333, 0.33333]

        # xy value of equal-energy LED mixture (cound by experiment)
        self.LED_white_xy = [0.386846372,  0.40224135]

        # only white we currently use
        self.white_xy  = self.D65_xy
        self.white_XYZ = self.xyY_to_XYZ(self.white_xy[0], self.white_xy[1])


        # to make xy hue angle comparable with conventional HSV angle,
        # make hue go clockwise, and normalize to red = 0.0
        # these are normalize 0-1, multiply by 360 for degrees

     
        # pre-compute  hue angles in xy space.
        self.r_abs_angle = self.xy_to_hue(self.red630_xy, 0.0)
        r_angle = self.r_abs_angle
        self.r_angle = 0.
        self.a_angle = self.xy_to_hue(self.org590_xy, r_angle)
        self.g_angle = self.xy_to_hue(self.grn530_xy, r_angle)
        self.b_angle = self.xy_to_hue(self.blu475_xy, r_angle)

        # hue halfway between red and amber
        RAH_xy = [0.5*(self.red630_xy[0] + self.org590_xy[0]),
                  0.5*(self.red630_xy[1] + self.org590_xy[1])]
        self.RAH = self.xy_to_hue(RAH_xy, r_angle)
        AGH_xy = [0.5*(self.grn530_xy[0] + self.org590_xy[0]),
                  0.5*(self.grn530_xy[1] + self.org590_xy[1])]
        self.AGH = self.xy_to_hue(AGH_xy, r_angle)
        GBH_xy = [0.5*(self.grn530_xy[0] + self.blu475_xy[0]),
                  0.5*(self.grn530_xy[1] + self.blu475_xy[1])]
        self.GBH = self.xy_to_hue(GBH_xy, r_angle)

        BRH_xy = [0.5*(self.red630_xy[0] + self.blu475_xy[0]),
                  0.5*(self.red630_xy[1] + self.blu475_xy[1])]
        self.BRH = self.xy_to_hue(BRH_xy, r_angle)

        # ## explore whether RGB primaries are any better, they aren't!
        # ## xy  coordinates for BEST RGB colorspace
        # self.red_pri = [ 0.73519164,  0.26480836]
        # self.grn_pri = [ 0.21533613,  0.77415966]
        # self.blu_pri = [ 0.13012295,  0.03483607]

        # self.r_abs_angle = self.xy_to_hue(self.red_pri, 0.0)
        # r_angle = self.r_abs_angle
        # self.r_angle = 0.
        # self.a_angle = self.xy_to_hue(self.org590_xy, r_angle)
        # self.g_angle = self.xy_to_hue(self.grn_pri, r_angle)
        # self.b_angle = self.xy_to_hue(self.blu_pri, r_angle)
        # ## hue halfway between red and amber
        # RAH_xy = [0.5*(self.red_pri[0] + self.org590_xy[0]),
        #           0.5*(self.red_pri[1] + self.org590_xy[1])]
        # self.RAH = self.xy_to_hue(RAH_xy, r_angle)
        # AGH_xy = [0.5*(self.grn_pri[0] + self.org590_xy[0]),
        #           0.5*(self.grn_pri[1] + self.org590_xy[1])]
        # self.AGH = self.xy_to_hue(AGH_xy, r_angle)
        # GBH_xy = [0.5*(self.grn_pri[0] + self.blu_pri[0]),
        #           0.5*(self.grn_pri[1] + self.blu_pri[1])]
        # self.GBH = self.xy_to_hue(GBH_xy, r_angle)

        # BRH_xy = [0.5*(self.red_pri[0] + self.blu_pri[0]),
        #           0.5*(self.red_pri[1] + self.blu_pri[1])]
        # self.BRH = self.xy_to_hue(BRH_xy, r_angle)

        
        
    def gamut_RGBA_xy(self):
        return([self.red630_xy, self.grn530_xy, self.blu475_xy, self.org590_xy])

    def matrix_mult(self, M, inp):
        # matrix multiply input triple by M
        # obv easier to do with a library like numpy but doing
        # it long form here for easy translation
        a =  M[0][0] * inp[0] + M[0][1] * inp[1] + M[0][2] * inp[2] 
        b =  M[1][0] * inp[0] + M[1][1] * inp[1] + M[1][2] * inp[2] 
        c =  M[2][0] * inp[0] + M[2][1] * inp[1] + M[2][2] * inp[2] 
        return(a, b, c)

    
    def HSV_to_RGBA(self, h, s, v):
        """ Naive (original) implementation of HSV colorspace to red, green, blue and
        amber channels red and amber squashed into 1/3 the hue range
        (instead of 1/2 as in naive) HSV inputs and rgby outputs all
        floats between 0 and 1
        this code is described as in http://rotormind.com/blog/2015/Generating-RGBA-from-Hue/
        """
        # offset h so 0 is pure red, needed to match naive HSV-RGBA
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

    def HSV_to_RGBA_CIE(self, h, s, v, warp=True):
        """implementation of HSV colorspace to red, green, blue and
        amber channels, based on (rough!) angles in CIE space.
        # hue breakpoints are tweaked to be (roughly!) linear with
        warped input
        HSV inputs and rgby outputs all floats between 0 and 1;
        not checked or clipped.
        this code is described as in 
        http://rotormind.com/blog/2015/Generating-RGBA-from-Hue/
        """

        # h = h - 25/360.
        # if h < 0:
        #     h += 1.0

        # prewarp hue input so it's more linear with RGBA_to_hue() inverse
        if warp:
            h = self.prewarp(h)
            
        fudge = 0.0


        # wrap around
        h = h%1.0

        # hand-tweaked for linearity
        RH = 0.
        AH = 20./360.
        GH = 120./360.
        BH = 240./360.

        # halfway points between above hue angles
        RAH = AH/2
        #AGH = 90./360.
        AGH = 100./360.
        GBH = 180./360
        BRH = 300./360.

        
        # print(360*(RAH - RAH1))
        # print(360*(AGH - AGH1))
        # print(360*(GBH - GBH1))
        # print(360*(BRH - BRH1))

        if s == 0.0:
            return [v, v, v, v]

        # i = int(h*6.0) # what hue range are we in?

        #                         # v is top flat
        # b = v*(1.0 - s)         # bottom flat


        if h < RAH : # hue below red-amber boundary
            #print(f"{h}: octant 0")
            h_frac = h/RAH 
            return [v, b, b, v*(1.0 - s*(1-h_frac))]  # max r, a up 
        elif h < AH : # hue below amber boundary
            #print(f"{h}: octant 1")
            h_frac = (h - RAH)/(AH - RAH)
            return [v*(1.0 - s*h_frac), b, b, v]  # max a, r down 
        elif h < AGH : # hue below amber-green boundary
            #print(f"{h}: octant 2")
            h_frac = (h - AH)/(AGH - AH)
            return [b, v*(1.0 - s*(1-h_frac)), b, v]  # max a, g up 
        elif h < GH  : # hue below green boundary
            #print(f"{h}: octant 3")
            h_frac = (h - AGH)/(GH - AGH)
            return [b, v, b, v*(1.0 - s*h_frac)]   # max g, a down 
        elif h < GBH  : # hue below green-blue boundary
            #print(f"{h}: octant 4")
            h_frac = (h - GH)/(GBH - GH)
            return [b, v, v*(1.0 - s*(1-h_frac)), b]   # max g, b up 
        elif h < BH : # hue below blue boundary
            #print(f"{h}: octant 5")
            h_frac = (h - GBH)/(BH - GBH)
            return [b, v*(1.0 - s*h_frac), v, b]   # max b, g down 
        elif h < BRH : # hue below blue -red boundary
            #print(f"{h}: octant 6")
            h_frac = (h - BH)/(BRH - BH)
            return [v*(1.0 - s*(1-h_frac)), b, v, b]   # max b, r up 
        else : # # hue below red boundary (1.0 = 0.0)
            #print(f"{h}: octant 7")
            h_frac = (h - BRH)/(1.0 + RH - BRH)
            return [v, b, v*(1.0 - s*h_frac), b]  # max r, b down 


    def xyY_to_XYZ(self, x, y, Y=1.0):
        # from http://www.brucelindbloom.com/index.html?Eqn_xyY_to_XYZ.html

        if y == 0:
            return(0., 0., 0.)

        return( x*Y/y, Y, (1 - x - y)*Y/y)


    def hue_angle(self, xy_vec):
        import math
        ''' calculate hue angle in degrees from CIE space vector
        result is  0. - 1. '''
        # numpy version
        #angle = (360/(2*math.pi))*math.arctan2(xy_vec[0], xy_vec[1])
        angle = (1/(2*math.pi))*math.atan2(xy_vec[0], xy_vec[1])
        #sat = np.linalg.norm(xy_vec)
        return((angle + 1.)%1.)



    
    def RGBA_to_xy_geom(self, R, G, B, A):

        red = self.red630_xy
        grn = self.grn530_xy
        blu = self.blu475_xy
        org = self.org590_xy

        csum = (R + G + B + A)

        if csum < 0.001:
            return(cx.white_xy)

        x = (R*red[0] + G*grn[0] + B*blu[0] + A*org[0])/csum
        y = (R*red[1] + G*grn[1] + B*blu[1] + A*org[1])/csum

        return([x, y])

        
    # def RGBA_to_xy(self, R, G, B, A):
    #     X, Y, Z = self.RGBA_to_XYZ(R, G, B, A)
    #     x, y, Y = self.XYZ_to_xyY(X, Y, Z)
    #     return([x, y])

    # def RGBA_to_XYZ(self, R, G, B, A):
    #     red_XYZ= self.RGB_to_XYZ(R,  0.,  0.)
    #     grn_XYZ= self.RGB_to_XYZ(0., G,   0.)
    #     blu_XYZ= self.RGB_to_XYZ(0., 0.,  B )

    #     #org590_xy = [0.575151311, 0.424232235]
    #     org_XYZ = self.xyY_to_XYZ(self.org590_xy[0], self.org590_xy[1], A)


    #     X = red_XYZ[0] + grn_XYZ[0] + blu_XYZ[0] + org_XYZ[0]
    #     Y = red_XYZ[1] + grn_XYZ[1] + blu_XYZ[1] + org_XYZ[1]
    #     Z = red_XYZ[2] + grn_XYZ[2] + blu_XYZ[2] + org_XYZ[2]

    #     #return(X, Y, Z)
    #     return(X/Y, 1, Z/Y)

    # def RGBA_to_XYZ2(self, R, G, B, A):
    #     rgb_XYZ= self.RGB_to_XYZ(R,  G,  B)

    #     #org590_xy = [0.575151311, 0.424232235]
    #     org_XYZ = self.xyY_to_XYZ(self.org590_xy[0], self.org590_xy[1], A)


    #     X = rgb_XYZ[0] + org_XYZ[0]
    #     Y = rgb_XYZ[1] + org_XYZ[1]
    #     Z = rgb_XYZ[2] + org_XYZ[2]

    #     return(X, Y, Z)

                               
    def XYZ_to_hue(self, X, Y, Z, red_offset=None):
        ''' convert color/luminance specification in CIE XYZ 
        coordinates to hue angle '''
        x, y, Y = self.XYZ_to_xyY(X, Y, Z)
        return(self.xy_to_hue([x, y], red_offset))


    def prewarp(self, x, coeff=None):
        # this is a horrible hack to prewarp the hue scale.
        # using a polynomial fit
        # only valid for if x is inside 0-1
        # coeffs determined by polynomal fit to inverse
        # hue_to_RGBA() -> RGBA_to_hue_CIE hue difference
        # see github linearized_RGBA_to_hue.ipynb for coeff calc
        coeff_w = [0.00577841,  0.68687651, -2.50694912, 17.56006685
                   -37.34857261, 34.06893111 -11.45908063]
        
        # default to local warping coefficients if none specified
        if coeff is None:
            coeff = coeff_w

        # limit to 0-1 range so we don't blow sky high
        #x = x%1.0    
        xn = 1
        res = 0.
        for c in coeff:
            res += c * xn
            xn = xn *x
        return(res)
    

        
    def xy_to_hue(self, xy_vec, red_offset=None):
        import math
        
        #D50 colorpoint from Best RGB colorspace
        #white_xy = self.D65_xy

        white_xy = self.LED_white_xy



        tmp_vec = [white_xy[0] - xy_vec[0], white_xy[1] - xy_vec[1]]
        # calculate angle and offset so red is = zero = 360
        if red_offset is None:
            # horrible hack constant but it works
            #red_offset = 5/36
            red_offset = self.r_abs_angle
        return((2.0 - red_offset - self.hue_angle(tmp_vec))%1.)

    def RGBA_to_HSV_CIE(self, R, G, B, A):
        ''' convert RGBA to HSV values approximately'''
        maxc = max((R, G, B, A))
        minc = min((R, G, B, A))

        # X, Y, Z = self.RGBA_to_XYZ(R,G,B,A)
        # hue = self.XYZ_to_hue(X, Y, Z)
        xy = self.RGBA_to_xy_geom(R,G,B,A)
        hue = self.xy_to_hue(xy)

        # yes is is that simple, but yes it is perceptually wrong
        val = maxc
        sat = 1. - minc
        return(hue, sat, val)

    
    def RGB_to_XYZ(self, r, g, b, gamma=True):
        gval = self.gamma
        if gamma:
            r = r**gval
            g = g**gval
            b = b**gval
        X, Y, Z = self.matrix_mult(self.M_XYZ_to_RGB,
                                   [r, g, b])
        return(X, Y, Z)

    def clip(self, x, clipval=1.0):
        if x > clipval:
            return clipval
        elif x < 0:
            return 0.
        else:
            return x
        
    def XYZ_to_RGB(self, x, y, z, gamma=True, clip=True):
        gval = 1/self.gamma
        R, G, B = self.matrix_mult(self.M_RGB_to_XYZ,[x, y, z])
        
        if gamma:
            if R > 0.:
                R = R**gval
            if G > 0.:
                G = G**gval
            if B > 0.:
                B = B**gval
        # clip because some of these results may be out of [0-1] range
        if clip:
            R = self.clip(R)
            G = self.clip(G)
            B = self.clip(B)
                
        return(R, G, B)


    
    def RGBA_to_rgbOLD(self, r, g, b, a):
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
        RGBA_RGB = self.matrix_mult(self.M_XYZ_to_RGB,
                                  [X, Y, Z])


        maxc = max(RGBA_RGB)
        #print(f"maxc is {maxc}")
        R = RGBA_RGB[0]
        G = RGBA_RGB[1]
        B = RGBA_RGB[2]
            
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

    def rgb_to_HSV(self, r, g, b):
        # from https://www.niwa.nu/2013/05/math-behind-colorspace-conversions-rgb-hsl/
        #https://stackoverflow.com/questions/359612/how-to-convert-rgb-color-to-HSV/1626175#1626175
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



    
    def print_HSV(self, hue, sat, val, extra=""):
        print(f" HSV: {hue:.3f} {sat:.3f} {val:.3f} {extra}")            

    def print_RGBA(self, r, g, b, a, extra=""):
        print(f"RGBA: {r:.3f} {g:.3f} {b:.3f} {a:.3f} {extra}")        

        
if __name__ == '__main__':

    cx = ColorXform()

    steps = 16

    maxhdiff = 0.0
    for hue_degrees in range(0,360,10):
        hue = float(hue_degrees)/360.
        sat = 0.7
        val = 1.0
        print(f" HSV: {hue:.3f} {sat:.3f} {val:.3f}")        
        r, g, b, a = cx.HSV_to_RGBA(hue, sat, val)
        #print(f"RGBA: {r:.3f} {g:.3f} {b:.3f} {a:.3f}")        
        R, G, B = cx.RGBA_to_rgb(r,g,b,a)
        #print(f"RGB: {R:.3f} {G:.3f} {B:.3f}")
        H, S, V = cx.rgb_to_HSV(R, G, B)
        print(f" HSV: {H:.3f} {S:.3f} {V:.3f}")
        huediff = abs(hue-H)
        if huediff > 0.5:
            huediff = abs(1-huediff)
        
        if huediff > maxhdiff:
            maxhdiff = huediff
        print(f"hue diff: {huediff}")

        print("===============")
        print(f"max hue diff: {maxhdiff}")

        
        
