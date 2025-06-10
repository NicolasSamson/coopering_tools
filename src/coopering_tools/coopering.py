from courbe import ExpectedShape
from pentagone import trapezoid

import shapely as sh
from shapely.geometry import Polygon, Point
import numpy as np
from scipy.spatial.transform import Rotation as R
from matplotlib import pyplot as plt
from shapely.geometry import Polygon, Point
from shapely.ops import unary_union
from copy import deepcopy
class ExpectedCoopering:
    def __init__(self, courbe):
        self.courbe = courbe


    def compute_from_saw_angle(self,):
        pass

    def compute_from_number_of_block(self, number_of_block, jeu_epaisseur=0.01):
        """
        Compute the expected coopering shape based on the number of blocks.
        """
        if not isinstance(self.courbe, ExpectedShape):
            raise TypeError("Expected courbe to be an instance of ExpectedShape")

        self.angle_pentagone = (self.courbe.angle / number_of_block)
        self.number_block = number_of_block
        # test
        c = self.courbe.radius - self.courbe.thickness / 2 - jeu_epaisseur 
        a = c * np.cos(self.angle_pentagone / 2)
        thickness_p  = self.courbe.radius + self.courbe.thickness / 2 + jeu_epaisseur  - a 
        r_cp =  a +thickness_p / 2
        h = thickness_p
        b = 2 * c * np.sin(self.angle_pentagone / 2)
        B = 2 * (self.courbe.radius + self.courbe.thickness / 2 + jeu_epaisseur )* np.tan(self.angle_pentagone / 2)

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
    def plot_coopering(self, ax=None, unit="in"):
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

    def plot_h_u_cote(self, ax, cote, under=False,unit="in"):
        
        if under:
            vertical_decalage = -0.001
            vertical_shift = -0.005
        else:
            vertical_decalage = 0.001
            vertical_shift = 0.005   
        if unit == 'm':
            scale = 1
        elif unit == 'cm':
            scale = 100
        elif unit == 'mm':
            scale = 1000
        elif unit == 'in':
            scale = 39.3701

        straigth_line = np.array([[0,vertical_decalage ], [0, vertical_shift],[0, 2*vertical_shift]]).T + cote[:2,0].reshape(2,1)

        straigth_line_2 = np.array([[0,vertical_decalage ], [0, vertical_shift],[0, 2*vertical_shift]]).T + cote[:2,1].reshape(2,1)

        cote_central = deepcopy(cote) 
        cote_central[1,:] = cote_central[1,:] + vertical_shift
        distance = np.linalg.norm(cote[:,0] - cote[:,1])

        ax.plot(straigth_line[0,:]*scale, straigth_line[1,:]*scale, color='black', linewidth=2)
        ax.plot(straigth_line_2[0,:]*scale, straigth_line_2[1,:]*scale, color='black', linewidth=2)
        ax.plot(cote_central[0,:]*scale, cote_central[1,:]*scale, color='black', linewidth=2)
        #ax.plot(cote_central[0,:], cote_central[1,:], color='black', linewidth=2)
        diff_cote = np.diff(cote_central[0,:])
        ax.text(straigth_line[0,1]*scale+ diff_cote[0]/2*scale, straigth_line[1,1]*scale , f"{np.round(distance*scale,4)}", fontsize=12, ha='center', va='bottom')

        return ax
    


    def plot_v_u_cote(self, ax, cote, left=False,unit="in"):
        
        if left:
            horizontal_decalage = -0.001
            horizontal_shift = -0.005
        else:
            horizontal_decalage = 0.001
            horizontal_shift = 0.005   
        if unit == 'm':
            scale = 1
        elif unit == 'cm':
            scale = 100
        elif unit == 'mm':
            scale = 1000
        elif unit == 'in':
            scale = 39.3701

        straigth_line = np.array([[horizontal_decalage,0], [horizontal_shift,0],[2*horizontal_shift,0]]).T + cote[:2,0].reshape(2,1)

        straigth_line_2 = np.array([[horizontal_decalage,0 ], [horizontal_shift,0],[2*horizontal_shift,0]]).T 
        straigth_line_2[0,:] = straigth_line_2[0,:] + cote[0,0] 
        straigth_line_2[1,:] = straigth_line_2[1,:] + cote[1,1]

        cote_central = cote 
        cote_central[0,1] = cote[0,0]
        cote_central[0,:] = cote_central[0,:]+ horizontal_shift
        distance = np.linalg.norm(cote[1,0] - cote[1,1])

        print(cote*scale)
        print(straigth_line*scale)
        print(straigth_line_2*scale)
        print(cote_central*scale)
        ax.plot(straigth_line[0,:]*scale, straigth_line[1,:]*scale, color='black', linewidth=2)
        ax.plot(straigth_line_2[0,:]*scale, straigth_line_2[1,:]*scale, color='black', linewidth=2)
        ax.plot(cote_central[0,:]*scale, cote_central[1,:]*scale, color='black', linewidth=2)
        #ax.plot(cote_central[0,:], cote_central[1,:], color='black', linewidth=2)
        diff_cote = np.diff(cote_central[1,:])
        ax.text((straigth_line[0,1]+horizontal_shift)*scale, straigth_line[1,1]*scale+ diff_cote[0]/2*scale , f"{np.round(distance*scale,4)}", fontsize=12, ha='center', va='bottom')

        ax.text(0,0,f"Angle chanfrein {np.round(self.angle_pentagone*180/np.pi/2,2)}°", fontsize=12, ha='center', va='bottom')
        return ax
    def plot_the_draft(self,ax, unit="in"):
        if unit == 'm':
            scale = 1
        elif unit == 'cm':
            scale = 100
        elif unit == 'mm':
            scale = 1000
        elif unit == 'in':
            scale = 39.3701

        if ax is None:
            fig,ax = plt.subplots(2,1)
        
        self.trapezoid.plot_trapezoid(ax=ax,unit=unit)
        ax.set_aspect('equal')
        
        ax.set_title("Draft of Coopering Shape")

        
        b = self.trapezoid.polygon_abs[:2, 0:2]
        h = self.trapezoid.polygon_abs[:2, 1:3]
        B = self.trapezoid.polygon_abs[:2, 2:4]
        
        plot_h_u_cote = self.plot_h_u_cote(ax, b, under=False,unit=unit)
        plot_h_u_cote = self.plot_h_u_cote(ax, B, under=True,unit=unit)
        plot_v_u_cote = self.plot_v_u_cote(ax, h, left=False,unit=unit)
        ax.set_xlabel(f"X ({unit})")
        ax.set_ylabel(f"Y ({unit})")
        #ax.set_title("Draft of Coopering Shape with Dimensions")
        
        return ax

    def plot_all(self, unit="in"):
        """
        Plot the coopering shape with dimensions.
        """
        fig, ax = plt.subplots(2, 1, figsize=(10, 10))
        
        ax[0] = self.plot_coopering(unit=unit, ax=ax[0])
        ax[1] = self.plot_the_draft(unit=unit, ax=ax[1])
        ax[0].set_title("Expected Coopering Shape {nb_block} blocks, {angle} degre".format(nb_block=self.number_block, angle=np.round(np.rad2deg(self.angle_pentagone/2),3)))
   
        plt.show()

if __name__ == "__main__":

    center = [0, 0]
    radius = 0.0254*7.5
    angle = 2*np.pi  # 90 degrees in radians
    init_angle = 0
    thickness = 0.0254
    nb_block = 18
    jeu_epaisseur= 0#1/8*0.0254
    expected_shape = ExpectedShape(center, radius, angle, init_angle, thickness)
    coopering = ExpectedCoopering(expected_shape)
    coopering.compute_from_number_of_block(nb_block,jeu_epaisseur=jeu_epaisseur)
    coopering.compute_area_diff()
    coopering.plot_all(unit="in")
    #ax = coopering.plot_coopering()
     #plt.show()

    