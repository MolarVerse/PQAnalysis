import pytest
import numpy as np

from PQAnalysis.atomic_system import AtomicSystem
from PQAnalysis.core import Atom
from PQAnalysis.io import gen2xyz, read_gen_file, write_gen_file, xyz2gen
from PQAnalysis.io.gen_file.exceptions import GenFileReaderError

from . import pytestmark  # pylint: disable=unused-import



class TestGenFileConversion:

    @pytest.mark.usefixtures("tmpdir")
    def test_xyz2gen_preserves_atom_name_case(self):
        with open("input.xyz", "w", encoding="utf-8") as file:
            print("3", file=file)
            print("", file=file)
            print("C 0.0 0.0 0.0", file=file)
            print("Cl 1.0 0.0 0.0", file=file)
            print("Br 0.0 1.0 0.0", file=file)

        xyz2gen("input.xyz", output="output.gen")

        with open("output.gen", "r", encoding="utf-8") as file:
            lines = file.read().splitlines()

        assert lines[0] == "3 C"
        assert lines[1] == "C Cl Br"
        assert [line.split()[1] for line in lines[2:]] == ["1", "2", "3"]

        system = read_gen_file("output.gen")
        assert [atom.name for atom in system.atoms] == ["C", "Cl", "Br"]

    @pytest.mark.usefixtures("tmpdir")
    def test_write_gen_uses_atom_names_for_repeated_elements(self):
        atoms = [Atom("C1", "C"), Atom("C2", "C"), Atom("C1", "C")]
        system = AtomicSystem(
            atoms=atoms,
            pos=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                ]
            ),
        )

        write_gen_file(system, "output.gen", periodic=False)

        with open("output.gen", "r", encoding="utf-8") as file:
            lines = file.read().splitlines()

        assert lines[1] == "C1 C2"
        assert [line.split()[1] for line in lines[2:]] == ["1", "2", "1"]

    @pytest.mark.usefixtures("tmpdir")
    def test_xyz2gen2xyz_preserves_periodic_pq_atom_labels(self):
        with open("input.xyz", "w", encoding="utf-8") as file:
            print("2 10.0 11.0 12.0 90.0 90.0 90.0", file=file)
            print("", file=file)
            print("C1 0.0 0.0 0.0", file=file)
            print("C2 1.0 0.0 0.0", file=file)

        xyz2gen("input.xyz", output="output.gen")
        gen2xyz("output.gen", output="output.xyz")

        with open("output.xyz", "r", encoding="utf-8") as file:
            lines = file.read().splitlines()

        header = lines[0].split()
        assert header[0] == "2"
        assert np.allclose([float(value) for value in header[1:]], [
            10.0,
            11.0,
            12.0,
            90.0,
            90.0,
            90.0,
        ])
        assert [line.split()[0] for line in lines[2:]] == ["C1", "C2"]

    @pytest.mark.usefixtures("tmpdir")
    def test_gen2xyz_preserves_atom_name_case(self):
        with open("input.gen", "w", encoding="utf-8") as file:
            print("3 C", file=file)
            print("C Cl Br", file=file)
            print("1 1 0.0 0.0 0.0", file=file)
            print("2 2 1.0 0.0 0.0", file=file)
            print("3 3 0.0 1.0 0.0", file=file)

        gen2xyz("input.gen", output="output.xyz")

        with open("output.xyz", "r", encoding="utf-8") as file:
            lines = file.read().splitlines()

        assert [line.split()[0] for line in lines[2:]] == ["C", "Cl", "Br"]

    @pytest.mark.usefixtures("tmpdir")
    def test_read_gen_truncated_coordinates_raises(self):
        with open("truncated.gen", "w", encoding="utf-8") as file:
            print("4 C", file=file)
            print("C H", file=file)
            print("1 1 0.0 0.0 0.0", file=file)
            print("2 2 1.0 0.0 0.0", file=file)

        with pytest.raises(GenFileReaderError) as exception:
            read_gen_file("truncated.gen")

        assert "truncated.gen" in str(exception.value)
        assert "4" in str(exception.value)
        assert "2" in str(exception.value)

    @pytest.mark.usefixtures("tmpdir")
    def test_read_gen_truncated_cell_raises(self):
        with open("truncated_cell.gen", "w", encoding="utf-8") as file:
            print("2 S", file=file)
            print("C H", file=file)
            print("1 1 0.0 0.0 0.0", file=file)
            print("2 2 1.0 0.0 0.0", file=file)
            print("0.0 0.0 0.0", file=file)
            print("10.0 0.0 0.0", file=file)
            print("0.0 10.0 0.0", file=file)

        with pytest.raises(GenFileReaderError) as exception:
            read_gen_file("truncated_cell.gen")

        assert "truncated_cell.gen" in str(exception.value)
        assert "4 cell lines" in str(exception.value)
        assert "3" in str(exception.value)

    @pytest.mark.usefixtures("tmpdir")
    def test_read_gen_out_of_range_atom_id_raises(self):
        for filename, bad_id in [("id_high.gen", 5), ("id_low.gen", 0)]:
            with open(filename, "w", encoding="utf-8") as file:
                print("2 C", file=file)
                print("C H", file=file)
                print("1 1 0.0 0.0 0.0", file=file)
                print(f"2 {bad_id} 1.0 0.0 0.0", file=file)

            with pytest.raises(GenFileReaderError) as exception:
                read_gen_file(filename)

            assert filename in str(exception.value)
            assert str(bad_id) in str(exception.value)

    @pytest.mark.usefixtures("tmpdir")
    def test_read_gen_well_formed_periodic_parses(self):
        with open("periodic.gen", "w", encoding="utf-8") as file:
            print("2 S", file=file)
            print("C H", file=file)
            print("1 1 0.0 0.0 0.0", file=file)
            print("2 2 1.0 0.0 0.0", file=file)
            print("0.0 0.0 0.0", file=file)
            print("10.0 0.0 0.0", file=file)
            print("0.0 10.0 0.0", file=file)
            print("0.0 0.0 10.0", file=file)

        system = read_gen_file("periodic.gen")

        assert [atom.name for atom in system.atoms] == ["C", "H"]
        assert np.allclose(system.pos, [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        assert np.allclose(system.cell.box_lengths, [10.0, 10.0, 10.0])
