"""This program is a five-shell vesicle geometry calculator.

It calculates:
1. Shell radii
2. Shell volumes
3. Inner and outer surface areas
4. Shell symmetry
5. Total membrane volume consistency
"""

import itertools
import math

# Input parameters (units: Angstroms)
vesicle = {
    "core_radius": 305.0,
    "thicknesses": [10, 8, 18, 8, 10]
    # t1 = t5, t2 = t4, t3 = independent
}


def validate_symmetry(thicknesses):
    """Validate the thickness symmetry of the five-shell model."""

    # Compare shell 1 with shell 5
    condition1 = math.isclose(
        thicknesses[0], thicknesses[4], rel_tol=1e-6
    )

    # Compare shell 2 with shell 4
    condition2 = math.isclose(
        thicknesses[1], thicknesses[3], rel_tol=1e-6
    )

    return condition1 and condition2


def calculate_radii(core_radius, thicknesses):
    """Return the core radius and five shell boundaries."""

    # Accumulate each thickness starting at the core radius
    result = list(
        itertools.accumulate(
            thicknesses,
            initial=core_radius
        )
    )

    return result


def sphere_volume(radius):
    """Calculate the volume of a sphere"""
    return (4 / 3) * math.pi * (radius ** 3)


def shell_volume(inner_radius, outer_radius):
    """Calculate the volume between two spherical boundaries"""
    return sphere_volume(outer_radius) - sphere_volume(inner_radius)


def sphere_area(radius):
    """Calculate the surface area of a sphere"""
    return 4 * math.pi * (radius ** 2)


def calculate_shell_volumes(radii):
    """Calculate the volume of each of the five shells"""

    volumes = []

    # Six radii define five shells
    for i in range(len(radii) - 1):

        # Find the boundaries of the current shell
        inner_radius = radii[i]
        outer_radius = radii[i + 1]

        # Calculate the shell volume
        volume = shell_volume(inner_radius, outer_radius)

        # Store the result
        volumes.append(volume)

    return volumes


def calculate_shell_areas(radii):
    """Calculate inner and outer surface areas for each shell."""

    areas = []

    for i in range(len(radii) - 1):

        inner_radius = radii[i]
        outer_radius = radii[i + 1]

        # Calculate the two boundary areas
        inner_area = sphere_area(inner_radius)
        outer_area = sphere_area(outer_radius)

        # Store both areas as a tuple
        areas.append((inner_area, outer_area))

    return areas


def validate_volume(radii, shell_volumes):
    """Check that shell volumes sum to total membrane volume."""

    # Sum the individual shell volumes
    calculated_volume = sum(shell_volumes)

    # Calculate total volume between core and outer boundary
    total_membrane_volume = shell_volume(
        radii[0], radii[-1]
    )

    # Account for floating-point rounding
    return math.isclose(
        calculated_volume,
        total_membrane_volume,
        rel_tol=1e-9,
        abs_tol=1e-6
    )


def analyze_vesicle(core_radius, thicknesses):
    """Calculate geometry and validate the five-shell vesicle."""

    # Step 1: Check the inputs
    if not math.isfinite(core_radius) or core_radius <= 0:
        raise ValueError("Core radius must be positive and finite.")

    if len(thicknesses) != 5:
        raise ValueError("Exactly five shell thicknesses are required.")

    if any(not math.isfinite(t) or t <= 0 for t in thicknesses):
        raise ValueError("Shell thicknesses must be positive and finite.")

    # Calculate shell boundaries
    radii = calculate_radii(core_radius, thicknesses)

    # Calculate shell volumes
    volumes = calculate_shell_volumes(radii)

    # Calculate shell surface areas
    areas = calculate_shell_areas(radii)

    # Validate the constraints
    symmetry_valid = validate_symmetry(thicknesses)
    volume_valid = validate_volume(radii, volumes)

    # Display the results
    print("\nFive shell vesicle geometry")
    print("-" * 70)

    print(f"Core radius: {core_radius:.3f} Å")
    print(f"Outer radius: {radii[-1]:.3f} Å")
    print(f"Total membrane thickness: {radii[-1] - core_radius:.3f} Å")

    print("\nShell Geometry")
    print("-" * 70)

    print(
        f"{'Shell':<8}"
        f"{'Inner R':<12}"
        f"{'Outer R':<12}"
        f"{'Thickness':<12}"
        f"{'Volume (Å³)':<18}"
    )

    for i in range(5):
        print(
            f"{i + 1:<8}"
            f"{radii[i]:<12.3f}"
            f"{radii[i + 1]:<12.3f}"
            f"{thicknesses[i]:<12.3f}"
            f"{volumes[i]:<18.3f}"
        )

    print("\nShell Surface Areas")
    print("-" * 70)

    print(
        f"{'Shell':<8}"
        f"{'Inner Area (Å²)':<20}"
        f"{'Outer Area (Å²)':<20}"
    )

    for i in range(5):
        inner_area, outer_area = areas[i]

        print(
            f"{i + 1:<8}"
            f"{inner_area:<20.3f}"
            f"{outer_area:<20.3f}"
        )

    print("\nVolume information")
    print("-" * 70)

    print(f"Sum of shell volumes: {sum(volumes):.3f} Å³")
    print(
        f"Total membrane volume: "
        f"{shell_volume(radii[0], radii[-1]):.3f} Å³"
    )

    print("\nValidation")
    print("-" * 70)

    print(f"Shell 1 = Shell 5: {math.isclose(thicknesses[0], thicknesses[4], rel_tol=1e-6)}")
    print(f"Shell 2 = Shell 4: {math.isclose(thicknesses[1], thicknesses[3], rel_tol=1e-6)}")
    print(f"Overall symmetry: {'PASS' if symmetry_valid else 'FAIL'}")
    print(f"Volume accounting: {'PASS' if volume_valid else 'FAIL'}")

    overall_valid = symmetry_valid and volume_valid

    print(f"\nOverall geometry validation: {'PASS' if overall_valid else 'FAIL'}")

    # Return calculated data for future optimization
    return {
        "radii": radii,
        "volumes": volumes,
        "areas": areas,
        "symmetry_valid": symmetry_valid,
        "volume_valid": volume_valid,
        "overall_valid": overall_valid
    }


if __name__ == "__main__":

    results = analyze_vesicle(
        vesicle["core_radius"],
        vesicle["thicknesses"]
    )
