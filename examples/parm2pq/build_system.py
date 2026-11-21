"""Build solvated Trp-cage (PDB 1L2Y) with AMBER ff14SB + TIP3P for parm2pq.

usage: python build_system.py 1L2Y.pdb

Writes system.prmtop and system.inpcrd (minimised coordinates and box).
"""

import sys

import openmm
import parmed
from openmm import app, unit

pdb_file = sys.argv[1] if len(sys.argv) > 1 else "1L2Y.pdb"

# ff14SB for the protein, TIP3P for water and ions
forcefield = app.ForceField("amber14/protein.ff14SB.xml", "amber14/tip3p.xml")

# 1L2Y is an NMR ensemble; take the first model and add missing hydrogens
pdb = app.PDBFile(pdb_file)
modeller = app.Modeller(pdb.topology, pdb.getPositions(frame=0))
modeller.addHydrogens(forcefield, pH=7.0)
n_protein = modeller.topology.getNumAtoms()

# cubic water box with 1 nm padding; neutralise the +1 charge with Cl-
modeller.addSolvent(
    forcefield,
    model="tip3p",
    padding=1.0 * unit.nanometer,
    boxShape="cube",
    neutralize=True,
)

# no constraints: parm2pq needs every bond force constant (--shake-h adds SHAKE)
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0 * unit.nanometer,
    constraints=None,
    rigidWater=False,
)

# minimise with OpenMM, so PQ starts from a relaxed structure
context = openmm.Context(
    system,
    openmm.VerletIntegrator(0.001),
    openmm.Platform.getPlatformByName("CPU"),
)
context.setPositions(modeller.positions)
openmm.LocalEnergyMinimizer.minimize(context, 10.0, 5000)
state = context.getState(getPositions=True, getEnergy=True)
energy = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)

# convert the OpenMM system into an AMBER prmtop + inpcrd with ParmEd
structure = parmed.openmm.load_topology(
    modeller.topology, system, xyz=state.getPositions()
)
structure.save("system.prmtop", overwrite=True)
structure.save("system.inpcrd", format="rst7", overwrite=True)

residues = [residue.name for residue in modeller.topology.residues()]
box = structure.box[0]
print(f"protein atoms   {n_protein}")
print(f"water molecules {residues.count('HOH')}")
print(f"chloride ions   {residues.count('CL')}")
print(f"total atoms     {modeller.topology.getNumAtoms()}")
print(f"cubic box       {box:.2f} Angstrom")
print(f"minimised E     {energy:.1f} kcal/mol (OpenMM, PME)")
print("wrote system.prmtop and system.inpcrd")
