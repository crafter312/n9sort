# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 18:12:19 2026

@author: h.s.webb
"""

from scipy.interpolate import BSpline
import numpy as np


def parse_atima_splines(filepath):
    """Parses a 17-line ATIMA file into a list of 4 SciPy BSpline objects."""
    splines = []

    with open(filepath, "r") as f:
        # Read non-empty lines, ignoring leading/trailing whitespace
        lines = [line.strip() for line in f if line.strip()]

    # Skip the 4-line header
    spline_lines = lines[4:]

    # Each spline block consists of 3 lines: metadata, knots, control points
    for i in range(0, 12, 3):
        # Line 1: Metadata (n_coeffs, order) -> e.g., "497 5"
        meta = spline_lines[i].split()
        n_coeffs = int(meta[0])
        order = int(meta[1])  # order = degree + 1 (5 for quartic)
        degree = order - 1

        # Line 2: Knots (502 values)
        knots = np.fromstring(spline_lines[i + 1], sep=" ")

        # Line 3: Coefficients / Control Points (497 values)
        coeffs = np.fromstring(spline_lines[i + 2], sep=" ")

        # Sanity check array lengths against formula: len(knots) = len(coeffs) + order
        if len(coeffs) != n_coeffs:
            raise ValueError(
                f"Expected {n_coeffs} coefficients, but read {len(coeffs)}"
            )
        if len(knots) != len(coeffs) + order:
            raise ValueError(
                f"Knot count ({len(knots)}) does not match coeffs count + order ({len(coeffs) + order})"
            )

        # Construct the SciPy BSpline instance
        spl = BSpline(t=knots, c=coeffs, k=degree)
        splines.append(spl)

    return splines


# --- Example Usage & Testing ---
if __name__ == "__main__":
    file_path = "../LossFiles/ATIMA/H1-Al.dat/H1-Al.dat"

    # Extract all 4 quartic splines
    splines = parse_atima_splines(file_path)

    print(f"Successfully loaded {len(splines)} splines.\n")

    # ==========================================
    # Spline 1 Evaluation
    # ==========================================
    spline1 = splines[0]
    min_energy1 = spline1.t[spline1.k]
    max_energy1 = spline1.t[-spline1.k - 1]

    print("==========================================")
    print(f"Spline 1 Valid Energy Range: {min_energy1:.4e} to {max_energy1:.4e}")
    print("==========================================")

    test_energies1 = np.linspace(min_energy1, max_energy1, 10)
    results1 = spline1(test_energies1)

    print("Sample Evaluation Points:")
    for E, val in zip(test_energies1, results1):
        print(f"  E = {E:12.4f}  ->  f(E) = {val:14.6e}")
    print()

    # ==========================================
    # Spline 2 Evaluation
    # ==========================================
    spline2 = splines[1]
    min_energy2 = spline2.t[spline2.k]
    max_energy2 = spline2.t[-spline2.k - 1]

    print("==========================================")
    print(f"Spline 2 Valid Energy Range: {min_energy2:.4e} to {max_energy2:.4e}")
    print("==========================================")

    test_energies2 = np.linspace(min_energy2, max_energy2, 10)
    results2 = spline2(test_energies2)

    print("Sample Evaluation Points:")
    for E, val in zip(test_energies2, results2):
        print(f"  E = {E:12.4f}  ->  f(E) = {val:14.6e}")
    print()

    # ==========================================
    # Spline 3 Evaluation
    # ==========================================
    spline3 = splines[2]
    min_energy3 = spline3.t[spline3.k]
    max_energy3 = spline3.t[-spline3.k - 1]

    print("==========================================")
    print(f"Spline 3 Valid Energy Range: {min_energy3:.4e} to {max_energy3:.4e}")
    print("==========================================")

    test_energies3 = np.linspace(min_energy3, max_energy3, 10)
    results3 = spline3(test_energies3)

    print("Sample Evaluation Points:")
    for E, val in zip(test_energies3, results3):
        print(f"  E = {E:12.4f}  ->  f(E) = {val:14.6e}")
    print()

    # ==========================================
    # Spline 4 Evaluation
    # ==========================================
    spline4 = splines[3]
    min_energy4 = spline4.t[spline4.k]
    max_energy4 = spline4.t[-spline4.k - 1]

    print("==========================================")
    print(f"Spline 4 Valid Energy Range: {min_energy4:.4e} to {max_energy4:.4e}")
    print("==========================================")

    test_energies4 = np.linspace(min_energy4, max_energy4, 10)
    results4 = spline4(test_energies4)

    print("Sample Evaluation Points:")
    for E, val in zip(test_energies4, results4):
        print(f"  E = {E:12.4f}  ->  f(E) = {val:14.6e}")
    print()
    
    # ==========================================
    # SRIM Comparison Block
    # ==========================================
    test_energies = np.array([0.02, 0.1])  # Input MeV/u directly
    
    print("==========================================")
    print("      SRIM COMPARISON TEST POINTS         ")
    print("==========================================")
    
    for j, E in enumerate(test_energies):
        print(f"\nEnergy: {E:.4f} MeV/u")
        print("-" * 45)
        for i, spl in enumerate(splines):
            # Evaluate spline directly at E, then exponentiate the y-result
            log_val = spl(E)
            linear_val = np.exp(log_val) / 2
            print(
                f"  Spline {i+1}: loge_y*2 = {log_val:10.6f} -> Linear = {linear_val:12.6e}"
            )