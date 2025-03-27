# ORCA Tutorial: Geometry Optimization and Vertical Excitation Spectrum of Uracil

This tutorial outlines the steps to optimize the geometry of uracil and compute its vertical excitation spectrum using ORCA.

## Prerequisites

- **ORCA Installation**: Ensure ORCA is installed on your system. [Installation Guide](https://www.faccts.de/docs/orca/6.0/tutorials/installing-orca.html)
- **Input Files**: Prepare the molecular structure of uracil in an appropriate format (e.g., `.xyz` or `.pdb`).

## Step 1: Geometry Optimization

Create an input file named `uracil_opt.inp` with the following content:

```plaintext
! B3LYP def2-SVP Opt

* xyz 0 1
  C      0.000000    1.396792    0.000000
  C      1.209657    0.698396    0.000000
  N      1.209657   -0.698396    0.000000
  C      0.000000   -1.396792    0.000000
  C     -1.209657   -0.698396    0.000000
  N     -1.209657    0.698396    0.000000
  O      2.209657    1.298396    0.000000
  O     -2.209657   -1.298396    0.000000
  H      0.000000    2.482582    0.000000
  H      2.148774    2.482582    0.000000
  H     -2.148774   -2.482582    0.000000
  H      0.000000   -2.482582    0.000000
*
```

**Explanation**:

- `! B3LYP def2-SVP Opt`: Specifies the use of the B3LYP functional with the def2-SVP basis set for geometry optimization.
- The `* xyz 0 1` block defines the molecular geometry in Cartesian coordinates with a charge of 0 and a multiplicity of 1.

Run the optimization:

```bash
orca uracil_opt.inp > uracil_opt.out
```

After completion, the optimized geometry can be found in the output file `uracil_opt.out`.

## Step 2: Vertical Excitation Spectrum Calculation

Using the optimized geometry, create a new input file named `uracil_spectrum.inp`:

```plaintext
! B3LYP def2-SVP TDDFT

* xyzfile 0 1 uracil_opt.xyz
```

**Explanation**:

- `! B3LYP def2-SVP TDDFT`: Specifies the use of the B3LYP functional with the def2-SVP basis set to perform a Time-Dependent Density Functional Theory (TDDFT) calculation for excited states.
- `* xyzfile 0 1 uracil_opt.xyz`: Imports the optimized geometry from the `uracil_opt.xyz` file.

Run the TDDFT calculation:

```bash
orca uracil_spectrum.inp > uracil_spectrum.out
```

The vertical excitation energies and oscillator strengths will be reported in the `uracil_spectrum.out` file.

## References

For more detailed information, refer to the [ORCA 6.0 Tutorials](https://www.faccts.de/docs/orca/6.0/tutorials/).
