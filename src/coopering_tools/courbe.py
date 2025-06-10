import shapely as sh
from shapely.geometry import Polygon, Point
import numpy as np
from scipy.spatial.transform import Rotation as R
from matplotlib import pyplot as plt



class ExpectedShape:
    def __init__(self, center, radius, angle, init_angle,thickness,delta_angle = 0.001, degree=False):
        if isinstance(center, (list, tuple)):
            center = np.array(center)
        elif isinstance(center, np.ndarray) and center.ndim == 2:
            center = center.flatten()
        if isinstance(radius, (list, tuple)):
            radius = radius
        elif isinstance(radius, np.ndarray) and radius.ndim == 2:
            radius = radius

        if degree:
            self.angle = np.deg2rad(angle)
            self.init_angle = np.deg2rad(init_angle)
        else:
            self.angle = angle
            self.init_angle = init_angle
        self.center = center
        self.radius = radius
        self.thickness = thickness
        self.delta_angle = delta_angle
        self.tfransformation_matrix = np.eye(3)
        self.inside_r = radius - thickness/2
        self.outside_r = radius + thickness/2
    
        self.compute_center_lines()

        
    def compute_curve(self,radius):

        aprox_delta_angle = self.angle/((self.angle / self.delta_angle))

        angles = np.arange(self.init_angle, self.angle , aprox_delta_angle/2)

        points = np.zeros((2, len(angles)))
        for i, angle in enumerate(angles):
            #print(angle)
            rotation_matrix = R.from_euler('z', angle, degrees=False).as_matrix()[:2, :2]
            point = rotation_matrix @ np.array([radius, 0])
            points[:, i] = point
        
        return points
    def compute_center_lines(self):

        self.center_line = self.compute_curve(self.radius)
        self.inner_line = self.compute_curve(self.inside_r)
        self.outer_line = self.compute_curve(self.outside_r)
        #print(self.outer_line.shape, self.inner_line.shape, self.center_line.shape)
        print(self.inner_line.shape)
        print(self.inner_line[:,0])
        print(self.inner_line[:,-1])
        
        self.polygon = np.hstack((self.outer_line, np.flipud(self.inner_line.T).T, np.array([[self.outer_line[0, 0]], [self.outer_line[1, 0]]]).reshape(2, 1)))
        print(self.polygon)
        self.polygon_abs = self.polygon.copy()
        self.center_line_abs = self.center_line.copy()
        self.inner_line_abs = self.inner_line.copy()
        self.outer_line_abs = self.outer_line.copy()

    def set_center(self, center, rotation=0):
        if isinstance(center, (list, tuple)):
            center = np.array(center)
        elif isinstance(center, np.ndarray) and center.ndim == 2:
            center = center.flatten()

        transformation_matrix = np.eye(3)
        transformation_matrix[:2, 2] = center
        transformation_matrix[:2, :2] = R.from_euler('z', rotation, degrees=True).as_matrix()[:2, :2]
        
        self.polygon_abs = transformation_matrix[:2, :2] @ self.polygon + center.reshape(2, 1)
        self.center_line_abs = transformation_matrix[:2, :2] @ self.center_line + center.reshape(2, 1)
        self.inner_line_abs = transformation_matrix[:2, :2] @ self.inner_line + center.reshape(2, 1)
        self.outer_line_abs = transformation_matrix[:2, :2] @ self.outer_line + center.reshape(2, 1)
        
        self.center = center

        self.tfransformation_matrix = transformation_matrix

    def plot_lines(self, ax):
        ax.plot(self.center_line_abs[0, :], self.center_line_abs[1, :], 'r-', label='Center Line')
        ax.plot(self.inner_line_abs[0, :], self.inner_line_abs[1, :], 'b-', label='Inner Line')
        ax.plot(self.outer_line_abs[0, :], self.outer_line_abs[1, :], 'g-', label='Outer Line')
        #ax.plot(self.polygon_abs[0, :], self.polygon_abs[1, :], 'k-', label='Polygon')
        
    def plot_polygon(self, ax):
        ax.plot(self.polygon_abs[0, :], self.polygon_abs[1, :], 'k-', label='Polygon')
        #ax.fill(self.polygon_abs[0, :], self.polygon_abs[1, :], alpha=0.5, label='Filled Polygon')

    def plot_center(self, ax):
        ax.scatter(self.center[0], self.center[1], color='red', label='Center', marker='+')
        #ax.text(self.center[0], self.center[1], f'({self.center[0]:.2f}, {self.center[1]:.2f})', fontsize=12, ha='right')

    def plot_frame(self, ax):
        
        frame = np.array([[1, 0,],
                         [0, 1,],
                         [1, 1,]])
        
        self.frame_abs = self.tfransformation_matrix @frame 

        #ax.scatter(self.center[0, 2], self.frame_abs[1, 2], color='black', label='Frame Origin', marker='o')
        ax.quiver(self.tfransformation_matrix[0, 2], self.tfransformation_matrix[1, 2],
                   self.tfransformation_matrix[0, 0], self.tfransformation_matrix[1, 0],
                   angles='xy', scale_units='xy', scale=1, color='red', label='Frame X')
        ax.quiver(self.tfransformation_matrix[0, 2], self.tfransformation_matrix[1, 2],
                   self.tfransformation_matrix[0, 1], self.tfransformation_matrix[1, 1],
                   angles='xy', scale_units='xy', scale=1, color='blue', label='Frame Y')

    def plot_update(self, ax):
        self.plot_lines(ax)
        self.plot_polygon(ax)
        self.plot_center(ax)
        self.plot_frame(ax)
        #ax.set_title("Expected Shape")
        ax.set_xlabel("X-axis")
        ax.set_ylabel("Y-axis")
        

if __name__ == "__main__":

    center = [0, 0]
    radius = 1
    angle = np.pi 
    init_angle = 0.0
    thickness = 0.1
    expected_shape = ExpectedShape(center, radius, angle, init_angle, thickness)
    
    fig, ax = plt.subplots(1, 1)
    expected_shape.plot_update(ax)
    expected_shape.set_center([1, 2], rotation=0)
    expected_shape.plot_update(ax)

    expected_shape.set_center([2, 4], rotation=45)
    expected_shape.plot_update(ax)
    
    ax.set_aspect('equal')
    #ax.legend()
    plt.show()
    