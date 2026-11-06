"""This program is a five-shell vesicle geometry calculator"""

import itertools
import math

vesicle = { 
    "core_radius": 305.0,
    "thicknesses": [10, 8, 18, 8, 10] # t1 = t5, t2 = t4, t3 = independent
}

def validate_symmetry(thicknesses):
    condition1 = thicknesses[0] == thicknesses[4]
    condition2 = thicknesses[1] == thicknesses[3]

    return condition1 and condition2

def calculate_radii(core_radius, thicknesses):
    """This function returns a list of six radii: core radius
      and the five shell boundaries"""
    
    # Accumulate starting with core radius
    result = list(
        itertools.accumulate(
            thicknesses, 
            initial=core_radius
        )
    )
    return result # returns a list with six radii: core radius and five boundaries

def sphere_volume(radius):
    """This function calculates the volume of a sphere"""
    return (4/3)*math.pi*(radius ** 3)

def shell_volume(inner_radius, outer_radius):
    """This function calculates the difference between two sphere volumes"""
    return sphere_volume(outer_radius) - sphere_volume(inner_radius)

def sphere_area(radius):
    """This function calculates the area of a sphere"""
    return 4*math.pi*(radius ** 2)

if __name__ == "__main__":
    print(calculate_radii(vesicle["core_radius"], vesicle["thicknesses"]))
    print(validate_symmetry(vesicle["thicknesses"]))
                   
                