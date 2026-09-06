# Water tutorial fixture

A 25-frame isolated water molecule in a 10 Å cubic cell, plus the H₂O Hessian
used by the vibrational example. It is small enough to run every analysis in
the documentation. It is not a bulk-liquid benchmark.

From the repository root:

```bash
cd examples/water
pqanalysis rdf rdf.in
pqanalysis msd msd.in
pqanalysis vacf vacf.in
pqanalysis vibrations vibrations.in
pqanalysis check_momentum trajectory.vel --selection all --output momentum.dat
```

`window` values in `msd.in` and `vacf.in` are sized to this 25-frame file, not
to a production trajectory. See the documentation page *Example data*.
