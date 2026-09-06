# Water tutorial fixture

A 25-frame isolated water molecule in a 10 Å cubic cell, plus the H₂O Hessian
used by the vibrational example. Every analysis in the documentation runs on
it; the numbers are tutorial output, not bulk-liquid results.

From the repository root:

```bash
cd examples/water
pqanalysis rdf rdf.in
pqanalysis msd msd.in
pqanalysis vacf vacf.in
pqanalysis vibrations vibrations.in
pqanalysis check_momentum trajectory.vel --selection all --output momentum.dat
```

`window` in `msd.in` and `vacf.in` is sized to these 25 frames. See the
documentation page *Example data*.
