# Vibrational Analysis

This example runs the PQAnalysis input-file based vibrational analysis command.

```bash
pqanalysis vibrations input.in
```

The standalone entry point is equivalent:

```bash
vibrations input.in
```

The input file reads a restart structure, a Cartesian Hessian matrix and a moldescriptor file with partial charges:

```text
structure_file = h2o.rst
hessian_file = hessian.dat
moldescriptor_file = moldescriptor.dat
out_file = wavenumbers.dat
normal_modes_file = normal_modes.dat
modes_prefix = mode
modes_file = modes.xyz
modes = positive
modes_frames = 30
modes_amplitude = 0.25
modes_threshold = 1.0e-6
unit = kcal
hessian_sign = auto
```

`wavenumbers.dat` contains frequencies, IR intensities, reduced masses and force constants. `normal_modes.dat` contains the normal-mode matrix. The `mode-*.xyz` files are sinusoidal XYZ animations for the selected modes. `modes.xyz` is an extended XYZ file with all selected mode vectors and metadata.
