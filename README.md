# CQ-MFS
Convolution Quadrature methods with the Method of Fundamental Solutions for exterior acoustic scattering

**Author**: Ignacio Labarca, Seminar for Applied Mathematics, ETH Zurich.\
**email**:  ignacio.labarca@sam.math.ethz.ch\
**date**:   September 2020.  


Details can be found in the [report](https://math.ethz.ch/sam/research/reports.html?id=922).
## Model
We solve the exterior acoustic scattering problem: the wave equation with Dirichlet boundary conditions.

## Code
* `incident_field.m` is a function that computes a planewave incident field.
* `main_MFS_CQ.m` solves an exterior problem by multistep CQ and MFS. Plots the solution.
* `main_MFS_RKCQ.m` solves an exterior problem by multistage CQ and MFS. Plots the solution.
* `test_error_exterior.m` runs the test for the Method of Fundamental Solutions in the frequency domain (Helmholtz problems) for different wavenumbers.
* `test_interior_BDF_RK.m` runs the test for the CQ+MFS for interior problems with exact solutions by planewaves.
* `test_exterior_BDF_RK_overkill.m` runs the test for the CQ+MFS for exterior problems with an overkill solution.
