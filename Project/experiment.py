from solid import *
from solid.utils import *

import project_utils

SEGMENTS = 48

def box_with_hole(hole_radius = 0.3):
    """ Hole radius should be in [0, 0.5)  to give a connected shape. """
    height = 1
    width = 1
    thickness = 0.1
    hole_radius = width / 3  # [0, min(height, width)/2]

    box = cube([width, thickness, height], center=True)
    
    shape = box - rotate([90, 0, 0])(cylinder(r=hole_radius, h=thickness*2, center=True))

    return up(height/2)(shape)  # Align bottom with origin of z axis



if __name__ == '__main__':
    n_samples = 10

    for i in range(n_samples):
        geom = box_with_hole((i+1) / n_samples * 0.5)
        project_utils.save_to_scad_and_stl(geom, 'test_models/hole_' + str(i))
