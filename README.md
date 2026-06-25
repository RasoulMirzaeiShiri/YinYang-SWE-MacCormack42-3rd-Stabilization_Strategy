# Shallow Water Equations on a Yin–Yang Grid: Strategy 3

This repository contains the implementation of Strategy 3 for the numerical solution of the shallow water equations (SWEs) on a Yin–Yang overset grid.

The solver employs a fourth-order compact MacCormack scheme for spatial discretization and the classical fourth-order Runge–Kutta (RK4) method for time integration. The implementation corresponds to the third stabilization strategy proposed and analyzed in the associated study.

---

## Overview

High-order numerical methods on overset grids may suffer from stability issues caused by repeated interpolation and derivative evaluations across overlapping interfaces. Strategy III was developed to improve the robustness of the Yin–Yang framework while preserving high-order accuracy.

The present repository contains the complete implementation of this approach together with benchmark test cases used for validation.

---

## Numerical Method

* Governing equations: Shallow Water Equations (SWEs)
* Grid topology: Yin–Yang overset grid
* Spatial discretization: Fourth-order compact MacCormack scheme
* Time integration: Fourth-order Runge–Kutta (RK4)
* Grid coupling: High-order interpolation between Yin and Yang subgrids

---

## Strategy 3

Strategy 3 introduces two modifications to the standard Yin–Yang solution procedure:

### Enhanced Overlap Region

The overlap between the Yin and Yang subgrids is increased in order to improve information exchange across grid interfaces and reduce interpolation-induced inconsistencies.

### Yin-Centered Derivative Evaluation

At every Runge–Kutta stage:

1. Prognostic variables from both subgrids are transferred to the Yin grid.
2. Spatial derivatives are computed exclusively on the Yin grid.
3. Derivative quantities required by the Yang grid are interpolated from the Yin grid to the Yang grid.
4. The RK4 update is then performed using the resulting derivative fields.

This procedure eliminates the need for independent derivative evaluations on both component grids and improves numerical stability near overset interfaces.

---

## Validation

The implementation is verified using two classical benchmark test cases for the shallow water equations on the sphere.

The benchmark problems are used to evaluate:

* Numerical stability
* Solution accuracy
* Conservation properties
* Performance of Strategy III

---

## Related Repositories

This repository contains only the implementation of Strategy 3.

Implementations of Strategy 1 and Strategy 2 are provided in separate repositories, while a detailed comparison of all three approaches can be found in the associated publication.

---

## Repository Structure

```text
.
├── src/          # Source code
├── grid/         # Yin–Yang grid and interpolation utilities
├── tests/        # Benchmark test cases
├── results/      # Numerical results
└── README.md
```

---

## Citation

If you use this code in your research, please cite the associated publication describing the three stabilization strategies and their comparative assessment.

---

## License

See the LICENSE file for details.
