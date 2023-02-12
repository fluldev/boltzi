# Fitting Constants

This is an overview of fitting constants which are used in the literature
for approximations of the scattering rate tau^-1 dependency on temperature
and angular frequency. Due to the simulation using wave vector power dependcies
the constants must be appropriately modified with the group velocity using k=w/v_g.
Constants are always in the respective SI-Units. T is the temperature and
w the angular frequency of the phonon. L ist the abbreviation of longitudinal
and T of transversal. N stands for normal three phonon processes and
U for umklapp thee phonon scattering. I is impurity scattering and
IS isotope scattering.

## Silicon

### Scattering

**Type** | **Term** | **Constants** | **Note/Source**
--- | --- | --- | ---
LN | B_L T^3 w^2 | B_L=2e-24 | [^1]
TU | B_T T exp(-A/T) w^2 | B_T=1.73e-19, A=137.3 | [^1]
I | C w^4 | C=8e-45 | [^1]
--- | --- | --- | ---
U | B w^2 T exp(-C/T) | B=2.8e-19, C=140 | B\*4.5 for Debye Modell [^2]
I | D w^4 | D=1.32e-45 | [^2]

### Group Velocities

**Type** | **Group Velocity** | **Note/Source**
--- | --- | ---
L | 9.01e3 | [^1]
T | 5.23e3 | [^1]

### Surface Roughness

**RMS** | **Correlation Length** | **Note/Source**
--- | --- | ---
3-50e-10 | 60e-10 | Nanowire [^1]

[1] 10.1103/PhysRevLett.102.125503
[2] http://link.aip.org/link/doi/10.1063/1.4710993
