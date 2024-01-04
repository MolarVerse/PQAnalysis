import numpy as np
import pytest

from PQAnalysis.topology import Selection, Topology
from PQAnalysis.core import Atom, Element


class TestSelection:
    def test__selection(self):
        selection = Selection()
        assert selection.selection_object is None

        topology = Topology([Atom("C1", 6), Atom("C2", 6), Atom("C1", 6)])
        indices = selection.select(topology)
        assert np.all(indices == np.arange(0, 3))

        selection = Selection("C1")
        indices = selection.select(topology)
        assert np.all(indices == np.array([0, 2]))

        selection = Selection(np.array([1, 2]))
        indices = selection.select(topology)
        assert np.all(indices == np.array([1, 2]))

        selection = Selection([Atom("C1", 6), Atom("C6", 6)])
        indices = selection.select(topology)
        assert np.all(indices == np.array([0, 1, 2]))

        selection = Selection([Atom("C1", 6), Atom("C6", 6)])
        indices = selection.select(topology, use_full_atom_info=True)
        assert np.all(indices == np.array([0, 2]))

        selection = Selection(Atom("C6", 6))
        indices = selection.select(topology, use_full_atom_info=True)
        assert np.all(indices == np.array([]))

        selection = Selection(Atom("C2", 6))
        indices = selection.select(topology, use_full_atom_info=True)
        assert np.all(indices == np.array([1]))

        selection = Selection(Element(6))
        indices = selection.select(topology)
        assert np.all(indices == np.array([0, 1, 2]))

        selection = Selection([Element(7), Element('h')])
        indices = selection.select(topology)
        assert np.all(indices == np.array([]))

        selection = Selection([Element(7), Element('h')])
        with pytest.raises(ValueError) as exception:
            selection.select(topology, use_full_atom_info=True)
        assert str(
            exception.value) == "The use_full_atom_info parameter is not supported for Element objects."

        selection = Selection("C1")
        indices = selection.select(topology)
        assert np.all(indices == np.array([0, 2]))

        selection = Selection(selection)
        indices = selection.select(topology)
        assert np.all(indices == np.array([0, 2]))

        selection = Selection(["C1", "C2"])
        indices = selection.select(topology)
        assert np.all(indices == np.array([0, 1, 2]))

        selection = Selection("1")
        indices = selection.select(topology)
        assert np.all(indices == np.array([1]))

        selection = Selection("Atom(C1, 6)")
        indices = selection.select(topology, use_full_atom_info=True)
        assert np.all(indices == np.array([0, 2]))

        selection = Selection("Atom(C1, c)")
        indices = selection.select(topology, use_full_atom_info=True)
        assert np.all(indices == np.array([0, 2]))

        selection = Selection("Atom(C1, c)")
        indices = selection.select(topology)
        assert np.all(indices == np.array([0, 1, 2]))

        selection = Selection("Elem(6)")
        indices = selection.select(topology)
        assert np.all(indices == np.array([0, 1, 2]))

        selection = Selection("Elem(c)")
        indices = selection.select(topology)
        assert np.all(indices == np.array([0, 1, 2]))

        selection = Selection("0..2")
        indices = selection.select(topology)
        assert np.all(indices == np.array([0, 1, 2]))

        selection = Selection("0-2")
        indices = selection.select(topology)
        assert np.all(indices == np.array([0, 1, 2]))

        selection = Selection("0..2..2")
        indices = selection.select(topology)
        assert np.all(indices == np.array([0, 2]))

        selection = Selection("C1, 1")
        indices = selection.select(topology)
        assert np.all(indices == np.array([0, 1, 2]))

        selection = Selection("Elem(6) & 2")
        indices = selection.select(topology)
        assert np.all(indices == np.array([2]))

        selection = Selection("Elem(6) | 2")
        indices = selection.select(topology)
        assert np.all(indices == np.array([0, 1]))

        selection = Selection("Elem(6) | 2 & 1")
        indices = selection.select(topology)
        assert np.all(indices == np.array([1]))

        selection = Selection("Elem(6) | (2 & 1)")
        indices = selection.select(topology)
        assert np.all(indices == np.array([0, 1, 2]))
