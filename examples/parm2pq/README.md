# Solvated Trp-cage with parm2pq

Files for the tutorial *Classical MD of a solvated protein* in the
documentation. They take Trp-cage (PDB 1L2Y) with AMBER ff14SB and TIP3P from
the PDB file to a PQ NVT and NPT run.

| File | Role |
|---|---|
| `build_system.py` | OpenMM: hydrogens, water box, Cl-, minimisation; ParmEd writes `system.prmtop` + `system.inpcrd` |
| `nvt.in` | PQ input: 10 ps NVT at 300 K, 2 fs, SHAKE, Wolf Coulomb |
| `npt.in` | PQ input: 20 ps NPT at 1 bar from the NVT restart |
| `rdf.in` | PQAnalysis input: water O-O RDF of the NPT trajectory |
| `tleap.in`, `min.in` | Alternative build with AmberTools (tleap + sander) |

```bash
pip install "PQAnalysis[ff]" openmm
curl -O https://files.rcsb.org/download/1L2Y.pdb
python build_system.py 1L2Y.pdb
parm2pq system.prmtop -c system.inpcrd -d md --shake-h
cp nvt.in npt.in rdf.in md/
cd md
PQ nvt.in
PQ npt.in
pqanalysis convert npt.en -o npt.xvg --timestep 2
pqanalysis rdf rdf.in
```

The generated files are not tracked (see `.gitignore`).
