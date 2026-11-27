"""Prepare a protein in water for AMBER-to-PQ conversion."""

from dataclasses import dataclass
from math import isfinite
from pathlib import Path



@dataclass(frozen=True)
class PreparedProtein:

    """AMBER output paths and summary of the solvated model."""

    topology: Path
    coordinates: Path
    protein_atoms: int
    water_molecules: int
    ions: int
    total_atoms: int
    box_angstrom: tuple[float, float, float]
    energy_kcal_mol: float



def prepare_protein(
    pdb_file: str | Path,
    output_dir: str | Path,
    *,
    padding_nm: float = 1.0
) -> PreparedProtein:
    """Solvate the first PDB model with ff14SB/TIP3P and write an AMBER pair.

    Adds hydrogens at pH 7 and neutralizing ions, then minimizes with OpenMM.
    The output directory must already exist; existing outputs are not replaced.
    Missing heavy atoms and unsupported residues are not repaired.
    """
    if not isfinite(padding_nm) or padding_nm < 1:
        raise ValueError("Water padding must be at least 1 nm.")
    output_dir = Path(output_dir)
    topology = output_dir / "system.prmtop"
    coordinates = output_dir / "system.inpcrd"
    if topology.exists() or coordinates.exists():
        raise FileExistsError(
            "The output directory already contains an AMBER output."
        )
    try:
        import openmm
        import parmed
        from openmm import app, unit
    except ImportError as error:
        raise ImportError(
            "Install PQAnalysis[protein] for protein preparation."
        ) from error

    # ff14SB for the protein, TIP3P for water and ions
    forcefield = app.ForceField(
        "amber14/protein.ff14SB.xml", "amber14/tip3p.xml"
    )

    # 1L2Y is an NMR ensemble; take the first model and add missing hydrogens
    pdb = app.PDBFile(str(pdb_file))
    modeller = app.Modeller(pdb.topology, pdb.getPositions(frame=0))
    modeller.addHydrogens(forcefield, pH=7.0)
    n_protein = modeller.topology.getNumAtoms()

    # cubic water box with 1 nm padding; neutralise the +1 charge with Cl-
    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=padding_nm * unit.nanometer,
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
    energy = state.getPotentialEnergy().value_in_unit(
        unit.kilocalorie_per_mole
    )

    # convert the OpenMM system into an AMBER prmtop + inpcrd with ParmEd
    structure = parmed.openmm.load_topology(
        modeller.topology, system, xyz=state.getPositions()
    )
    structure.save(str(topology))
    structure.save(str(coordinates), format="rst7")

    residues = [residue.name for residue in modeller.topology.residues()]
    return PreparedProtein(
        topology,
        coordinates,
        n_protein,
        residues.count("HOH"),
        residues.count("NA") + residues.count("CL"),
        modeller.topology.getNumAtoms(),
        tuple(float(x) for x in structure.box[:3]),
        energy,
    )
