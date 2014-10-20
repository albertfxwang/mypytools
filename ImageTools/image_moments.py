#!/usr/bin/env python
# 
# Routines for computing image moments and things from them
#
from numpy import *

class moments:
    """ Define an object to compute and retain the image moments 
        x = first-moment center of the image in x
        y = first-moment center of the image in y
        r = a 2-d image holding the radial distance from the center
        x2 = second-moment in x
        y2 = second-moment in y
        xy = second-moment cross term
        cxx, cyy, cxy = ellipse parameters (see SExtractor manual)
        rel = 2-d image holding the elliptical radii corresponding
          to each pixel...rel=2 is ellipse at roughly the Kron Radius;
          r=2.5 is the default SExtractor MAG_AUTO aperture...I think.
          To plot this ellipse in matplotlib: contour(self.rel,[2.5])
    """
    def __init__(self,image,section=None,threshold=None,radius=None):
       """ Initialize the image moments 
           Inputs: 
             image 
             section = optional subsection [x0,x1,y0,y1] (python list)
             threshold = ignore pixels below this value
             radius = ignore pixels further than this distance from the center
           Computes:
             image first & second moments
             radius from the center (2-d image)
             SExtractor Kron ellipse parameters
           Returns:
             Instance of the class with these values already computed
       """
       if section != None:
           x0,x1,y0,y1 = section
	   self.img = image[x0:x1,y0:y1]
       else:
           self.img = image
       if radius != None:
           # Zero out pixels beyond the desired radius
           r2 = radius**2
           ylen,xlen = shape(self.img)
	   x0 = xlen/2.
	   y0 = ylen/2.
           mg = mgrid[0:ylen,0:xlen]
#	   print shape(mg[0])
#	   print shape(mg[1])
#	   print shape(self.img)
	   dd = choose( (mg[0]-x0)**2 + (mg[1]-y0)**2 < r2, (0,self.img))
           self.img = dd
       if threshold != None:
           self.img = choose(self.img > threshold,(0,self.img))
       nx = shape(self.img)[0]
       ny = shape(self.img)[1]
       self.g = mgrid[0:nx,0:ny]
       self.x = self.g[0]
       self.y = self.g[1]
       self.firstorder()
       self.secondorder()
       self.ellipse()
       self.radius()
       self.r_ellipse()

    # Compute radius from a given center
    def radius(self,xc=None,yc=None):
       """ Compute the radius from a given center
            Uses the first-moment center as the default
            Optional args: xc,yc allow specifying a different center
       """
       if xc == None:
           xc = self.x1
       if yc == None:
           yc = self.y1
       self.r = sqrt((self.x-xc)**2+(self.y-yc)**2)

    # Compute angle relative to the center in pixel coords
    def xyangle(self,xc=None,yc=None):
       """ Compute the position angle of each pixel relative to the center
            Uses the first-moment center as the default
            Optional args: xc,yc allow specifying a different center
       """
       if xc == None:
           xc = self.x1
       if yc == None:
           yc = self.y1
       dx = self.x-xc
       dy = self.y-yc
       self.angle = arctan2(dx,dy) # in radians
       self.sin = sin(self.angle)
       self.cos = cos(self.angle)
    
    # Compute first-order moments
    def firstorder(self):
       """ Compute the first-order moments and x,y center of an image """
       f = self.img
       x = self.x
       y = self.y
       self.x1 = sum(f*x)/sum(f)
       self.y1 = sum(f*y)/sum(f)
    
    # Compute second-order moments
    def secondorder(self):
       """ Compute the second-order moments an image """
       f = self.img
       x = self.x
       y = self.y
       self.x2 = sum(f*x**2)/sum(f) - self.x1**2
       self.y2 = sum(f*y**2)/sum(f) - self.y1**2
       self.xy = sum(f*x*y)/sum(f) - self.x1*self.y1

    # Compute ellipse parametrs (from SExtractor manual)
    def ellipse(self):
       """ Compute the ellipse parameters moments of an image
           See SExtractor manual. 
           a = major axis length
           b = minor axis length
           theta = position angle (0 == aligned along x axis I think)
           cxx, cyy, cxy = SExtractor ellipse parameters
       """
       f = self.img
       x = self.x
       y = self.y
       x2 = self.x2
       y2 = self.y2
       xy = self.xy
       self.a2 = (x2+y2) + sqrt(((x2-y2)/2.)**2 + xy**2)
       self.b2 = (x2+y2) - sqrt(((x2-y2)/2.)**2 + xy**2)
       self.a = sqrt(self.a2)
       self.b = sqrt(self.b2)
       tan2theta = 2* (xy/(x2-y2))
       self.theta = arctan(tan2theta)/2.
       denominator = sqrt(((x2-y2)/2)**2+xy**2)
       self.cxx = y2/denominator
       self.cyy = x2/denominator
       self.cxy = -2*xy/denominator

    # Compute ellipse radius parameter (SExtractor manual)
    def r_ellipse(self,xc=None,yc=None):
       """ Compute elliptical radius of each pixel (see SExtractor manual).
          To plot an ellipse at 2.5 x the ellipse computed from
          the second moments, in matplotlib: contour(self.rel,[2.5])
       """
       x = self.x
       y = self.y
       if xc == None:
           xc = self.x1
       if yc == None:
           yc = self.y1
       self.rel = sqrt(self.cxx*(x-xc)**2 +
		       self.cyy*(y-yc)**2 +
		       self.cxy*(x-xc)*(y-yc)
		       )

