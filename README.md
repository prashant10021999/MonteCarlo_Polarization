# MonteCarlo_Polarization
 
Monte Carlo Simulation for Polarized Light Transport

Overview

This repository contains a Monte Carlo simulation code for polarized light transport in scattering media. 
Compilation and Execution

Compilation

To compile the program, use the following command:

make

Execution

To run the program, execute:

./iquv

Output Description

The program outputs the total reflectance and transmittance for each polarization state. Additionally, it generates 100x100 pixel images to construct the Mueller matrix.

Mueller Matrix Components

The generated matrices are represented as:

HH, HV, HP, HR

VH, VV, VP, VR

PH, PV, PP, PR

RH, RV, RP, RR

Where:

H: Horizontally polarized (parallel to reference frame)

V: Vertically polarized (perpendicular to reference frame)

P: 45-degree polarized

R: Right circularly polarized

Example Test Cases

Case 1: Sphere Diameter = 2.00000 μm

Parameters:

Scattering coefficient (μₛ) = 827.81840 cm⁻¹

Anisotropy factor (g) = 0.76861

Particle density (ρ) = 0.01000 particles/μm³

Slab thickness = 0.00483 cm

Number of photons = 100,000

Medium refractive index (n_media) = 1

Sphere refractive index (n_sphere) = 1.59

Results:

Horizontal polarization (H) [1 1 0 0]

R = [0.29299, 0.08992, -0.00054, -0.00014]

T = [0.70701, 0.51252, -0.00094, 0.00133]

Vertical polarization (V) [1 -1 0 0]

R = [0.29306, -0.09181, -0.00049, -0.00004]

T = [0.70694, -0.51268, -0.00114, 0.00011]

45-degree polarization (P) [1 0 1 0]

R = [0.29254, 0.00313, -0.08724, -0.00070]

T = [0.70746, -0.00033, 0.51334, 0.00084]

Right circular polarization (R) [1 0 0 1]

R = [0.29285, -0.00081, -0.00037, -0.00875]

T = [0.70715, 0.00039, 0.00078, 0.51589]

Case 2: Sphere Diameter = 0.01000 μm

Parameters:

Scattering coefficient (μₛ) = 0.00000 cm⁻¹

Anisotropy factor (g) = 0.00051

Particle density (ρ) = 0.01000 particles/μm³

Slab thickness = 275959772.13944 cm

Medium refractive index (n_media) = 1

Sphere refractive index (n_sphere) = 1.59

Results:

Horizontal polarization (H) [1 1 0 0]

R = [0.68934, 0.28741, -0.00059, -0.00000]

T = [0.31066, 0.07206, 0.00035, -0.00000]

Vertical polarization (V) [1 -1 0 0]

R = [0.68858, -0.28715, 0.00055, -0.00000]

T = [0.31142, -0.07305, -0.00013, -0.00000]

45-degree polarization (P) [1 0 1 0]

R = [0.68850, -0.00040, -0.28743, 0.00000]

T = [0.31150, 0.00017, 0.07283, 0.00000]

Right circular polarization (R) [1 0 0 1]

R = [0.68870, -0.00065, 0.00032, -0.20367]

T = [0.31130, -0.00029, -0.00061, 0.04982]

Notes

Reflectance (R) and transmittance (T) values are normalized by the number of photons launched.

The results vary based on input parameters such as sphere diameter, scattering coefficient, and polarization state.
