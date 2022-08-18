
__all__ = ["Point", "Box", "Sphere", "Naca"]

from .base_statement import Statement

class Point(Statement):
    """
    Point primitive
    """
    def __init__(self, xloc:float=0, yloc:float=0, zloc:float=0):
        super(Point,self).__init__(content=f"point {xloc} {yloc} {zloc}")

class Box(Statement):
    """
    Box primitive
    """
    def __init__(self, xbase:float=0, ybase:float=0, zbase:float=0, dx:float=0, dy:float=0, dz:float=0):
        super(Box,self).__init__(content=f"box {xbase} {ybase} {zbase} {dx} {dy} {dz}")

class Sphere(Statement):
    """
    Sphere primitive
    """
    def __init__(self, xbeg:float=0, ybeg:float=0, zbeg:float=0, xend:float=0, yend:float=0, zend:float=0, radius:float=1.0):
        super(Sphere,self).__init__(content=f"sphere {xbeg} {ybeg} {zbeg} {xend} {yend} {zend} {radius}")

class Naca(Statement):
    """
    Naca primitive
    """
    def __init__(self, series:int=None):
        #super(Naca,self).__init__(content=f"naca {xbeg} {ybeg} {zbeg} {xend} {yend} {zend} {radius}")
        pass
        #naca series thickness camber maxloc offset sharpte