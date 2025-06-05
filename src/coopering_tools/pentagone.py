import shapely as sh
from shapely.geometry import Polygon, Point
import numpy as np
from scipy.spatial.transform import Rotation as R
from matplotlib import pyplot as plt


class trapezoid:
    def __init__(self, b, B, h):

        self.center = np.array([0, 0])
        self.polygon = np.array([[-B/2,h/2], 
                                [B/2,h/2], 
                                [b/2,-h/2],
                                [-b/2,-h/2],
                                [-B/2,h/2]]).T
        self.polygon_abs = self.polygon.copy()
        self.transformation_matrix = np.eye(3)

    
    def set_center(self, center, rotation=0):
        if isinstance(center, (list, tuple)):
            center = np.array(center)
        elif isinstance(center, np.ndarray) and center.ndim == 2:
            center = center.flatten()


        rotation_matrix = R.from_euler('z', rotation, degrees=True).as_matrix()
        self.transformation_matrix = np.eye(3)
        self.transformation_matrix[:2, :2] = rotation_matrix[:2, :2]
        self.transformation_matrix[:2, 2] = center
        self.center = center
        self.polygon_abs = self.transformation_matrix @ np.vstack((self.polygon, np.ones(self.polygon.shape[1])))

        
    def area(self):
        return self.polygon.area
    



def test_sans_rotation(test_translation=True,test_rotation=True):
    t = trapezoid(1, 2, 3)
    
    t.polygon_abs = t.polygon
    print("Polygon:", t.polygon)
    print("Polygon Absolute Coordinates:", t.polygon_abs)


    fig, ax = plt.subplots(1,1)

    ax.plot(t.polygon_abs[0, :], t.polygon_abs[1, :], 'r-')
    ax.scatter(t.center[0], t.center[1], color='red', label='(0,0, 0)',marker='+' )
    ax.set_aspect('equal')

    if test_translation:
        t.set_center([1, 2], rotation=0)
        ax.plot(t.polygon_abs[0, :], t.polygon_abs[1, :], 'b-')
        ax.scatter(t.center[0], t.center[1], color='blue', label='(1,2, 0)',marker='+' )
    if test_rotation:
        t.set_center([1, 2], rotation=45)
        ax.plot(t.polygon_abs[0, :], t.polygon_abs[1, :], 'g-')
        ax.scatter(t.center[0], t.center[1], color='green', label='(0,0, 45)',marker='+' )
    ax.legend()
    ax.set_title("Trapezoid Polygon")
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_ylim(-2, 4)
    ax.set_xlim(-2, 4)
    plt.show()

if __name__ == "__main__":
    

    test_sans_rotation()
  