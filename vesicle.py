""" This is a small geometry calculator which converts
fitted SAXS parameters into geometric quantities that can be
used later in constraint equations
"""
import math

# Parameters from the original SAXS fit (Angstroms)
Rc = 305.18261
t1 = 12.52633
t2 = 18.273004
t3 = 8.5981313

# Calculate each shell boundary
Ri = Rc + t1
Rf = Ri + t2
Ro = Rf + t3

def sphere_volume(radius):
    """ This function returns the volume of a sphere """
    return (4 / 3) * math.pi * radius**3

def sphere_area(radius):
    """This function returns the area of a sphere"""
    return 4 * math.pi * radius**2

# Calculate volumes of the three shells
V_inner = sphere_volume(Ri) - sphere_volume(Rc)
V_tails = sphere_volume(Rf) - sphere_volume(Ri)
V_outer = sphere_volume(Ro) - sphere_volume(Rf)

# Calculate reference interfacial areas
A_inner = sphere_area(Ri)
A_outer = sphere_area(Rf)

print("VESICLE GEOMETRY")
print(f"Core radius: {Rc:.3f} Å")
print(f"Outer radius: {Ro:.3f} Å")
print(f"Membrane thickness: {Ro - Rc:.3f} Å")

print("\nSHELL VOLUMES")
print(f"Inner headgroups: {V_inner:.3f} Å³")
print(f"Hydrocarbon region: {V_tails:.3f} Å³")
print(f"Outer headgroups: {V_outer:.3f} Å³")

print("\nINTERFACIAL AREAS")
print(f"Inner interface: {A_inner:.3f} Å²")
print(f"Outer interface: {A_outer:.3f} Å²")
