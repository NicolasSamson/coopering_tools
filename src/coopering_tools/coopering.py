from courbe import ExpectedShape
from pentagone import trapezoid

import shapely as sh
from shapely.geometry import Polygon, Point
import numpy as np
from scipy.spatial.transform import Rotation as R
from matplotlib import pyplot as plt
from shapely.geometry import Polygon, Point
from shapely.ops import unary_union

class ExpectedCoopering:
    def __init__(self, courbe):
        self.courbe = courbe


    def compute_from_saw_angle(self,):
        pass

    def compute_from_number_of_block(self, number_of_block):
        """
        Compute the expected coopering shape based on the number of blocks.
        """
        if not isinstance(self.courbe, ExpectedShape):
            raise TypeError("Expected courbe to be an instance of ExpectedShape")

        self.angle_pentagone = (self.courbe.angle / number_of_block)
        
        # test
        c = self.courbe.inside_r 
        a = c * np.cos(self.angle_pentagone / 2)
        thickness_p  = self.courbe.outside_r - a 
        r_cp =  a +thickness_p / 2
        h = thickness_p
        b = 2 * c * np.sin(self.angle_pentagone / 2)
        B = 2 * self.courbe.outside_r * np.tan(self.angle_pentagone / 2)

        # Assume that the angle is already in radians
        #b = 2 * self.courbe.inside_r * np.sin(self.angle_pentagone / 2)
        #B = 2 * self.courbe.outside_r * np.tan(self.angle_pentagone / 2)
        #h = self.courbe.thickness +self.courbe.radius * (1- np.cos(self.angle_pentagone / 2))
        print(self.angle_pentagone)
        print(h,B,b)
        self.trapezoid = trapezoid(b, B, thickness_p)
        print(self.trapezoid.polygon)
        list_pentagons = []
        list_pentagons_shapely = []
        for i in range(number_of_block):

            test =  trapezoid(b, B, h)
            #print(test.polygon)
            angle_spec = i * self.angle_pentagone + self.courbe.init_angle
            print("Angle spec:", angle_spec)
            #print(self.courbe.radius)
            position_center_trapezoid = self.courbe.center + r_cp * np.array([np.cos(angle_spec),
                                                                        np.sin(angle_spec)])
            #print(test.polygon)
            test.set_center(position_center_trapezoid, rotation=angle_spec+np.deg2rad(-90),degrees=False)
            list_pentagons.append(test)
            list_pentagons_shapely.append(Polygon(test.polygon_abs.T))
            #print(test.polygon_abs)
        self.list_pentagons = list_pentagons
        self.list_pentagons_shapely = list_pentagons_shapely
    def plot_coopering(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
        
        
        self.courbe.plot_update(ax=ax)
        
        for pentagon in self.list_pentagons:
            
            pentagon.plot_trapezoid(ax=ax)
    
        ax.set_aspect('equal')
        return ax
    
    def compute_area_diff(self):

        polygon = unary_union(self.list_pentagons_shapely)
        circle = Polygon(self.courbe.polygon_abs.T)

        print("Area of Coopering Shape:", polygon.area)
        print("Area of Circle:", circle.area)
        print("Difference in Area:", polygon.area - circle.area)        

if __name__ == "__main__":

    center = [0, 0]
    radius = 0.0254*10
    angle = 2*np.pi  # 90 degrees in radians
    init_angle = 0
    thickness = 0.0254
    nb_block = 18
    expected_shape = ExpectedShape(center, radius, angle, init_angle, thickness)
    
    coopering = ExpectedCoopering(expected_shape)
    coopering.compute_from_number_of_block(nb_block)
    coopering.compute_area_diff()

    ax = coopering.plot_coopering()
    ax.set_title("Expected Coopering Shape")

    plt.show()