import os

import numpy as np

from . import pytestmark

from PQAnalysis.io.nep.nep_writer import NEPWriter
from PQAnalysis.traj import Trajectory
from PQAnalysis.core import Atom, Cell
from PQAnalysis.atomic_system import AtomicSystem

INFO_FILE_CONTENT = """\
-----------------------------------------------------------------------------------------
|                                  PQ info file                                  |
-----------------------------------------------------------------------------------------
|   SIMULATION-TIME         5.00000 ps       TEMPERATURE           294.45308 K          |
|   PRESSURE             2875.71714 bar      E(TOT)            -186008.40141 kcal/mol   |
|   E(QM)             -186197.10854 kcal/mol N(QM-ATOMS)             0.00000 -          |
|   E(KIN)                188.70714 kcal/mol E(INTRA)                0.00000 kcal/mol   |
|   VOLUME               4310.41014 A^3      DENSITY                 0.83455 g/cm^3     |
|   MOMENTUM                2.5e-10 amuA/fs  LOOPTIME                0.90443 s          |
-----------------------------------------------------------------------------------------
"""



class TestNEPWriter:

    def test_write_from_files_with_test_ratio(self, tmpdir):
        xyz_frame = (
            "2 10.0 10.0 10.0 90.0 90.0 90.0\n#\n"
            "C 0.0 0.0 0.0\nH 1.0 1.0 1.0\n"
        )
        with open("md.xyz", "w", encoding="utf-8") as file:
            file.write(xyz_frame * 4)

        with open("md.info", "w", encoding="utf-8") as file:
            file.write(INFO_FILE_CONTENT)

        with open("md.instant_en", "w", encoding="utf-8") as file:
            for i in range(1, 5):
                file.write(
                    f"{i} 298.0 23000.0 -185898.0 -18590{i}.0 0.0 5.6 0.0 "
                    "3729.1 0.83 1.1e-16 1.0\n"
                )

        writer = NEPWriter("out.xyz")
        writer.write_from_files(["md"], test_ratio=0.25)

        assert np.isclose(writer.test_ratio, 0.25)
        assert np.isclose(writer.validation_ratio, 0.0)
        assert writer.n_train_frames == 3
        assert writer.n_test_frames == 1
        assert os.path.exists("out.xyz_train")
        assert os.path.exists("out.xyz_test")

    def test_write_from_files_with_total_ratios(self, tmpdir):
        xyz_frame = (
            "2 10.0 10.0 10.0 90.0 90.0 90.0\n#\n"
            "C 0.0 0.0 0.0\nH 1.0 1.0 1.0\n"
        )
        with open("md.xyz", "w", encoding="utf-8") as file:
            file.write(xyz_frame * 4)

        with open("md.info", "w", encoding="utf-8") as file:
            file.write(INFO_FILE_CONTENT)

        with open("md.instant_en", "w", encoding="utf-8") as file:
            for i in range(1, 5):
                file.write(
                    f"{i} 298.0 23000.0 -185898.0 -18590{i}.0 0.0 5.6 0.0 "
                    "3729.1 0.83 1.1e-16 1.0\n"
                )

        writer = NEPWriter("out.xyz")
        writer.write_from_files(["md"], total_ratios="3:1")

        assert np.isclose(writer.test_ratio, 0.25)
        assert np.isclose(writer.validation_ratio, 0.0)
        assert writer.n_train_frames == 3
        assert writer.n_test_frames == 1

    def test_write_from_trajectory(self, tmpdir):
        system = AtomicSystem(
            atoms=[Atom("C"), Atom("H")],
            pos=np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]),
            cell=Cell(10.0, 10.0, 10.0),
        )
        system.energy = -1.5
        trajectory = Trajectory([system, system])

        writer = NEPWriter("traj_out.xyz")
        writer.write_from_trajectory(trajectory)

        with open("traj_out.xyz", "r", encoding="utf-8") as file:
            lines = file.read().splitlines()

        assert len(lines) == 8
        assert lines[0] == "2"
        assert lines[4] == "2"
        assert lines[1].startswith("energy=")
        assert lines[2].startswith("C")
        assert lines[3].startswith("H")
