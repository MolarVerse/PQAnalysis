from .. import pytestmark

from PQAnalysis.topology import JCoupling



class TestJCoupling:

    def test__init__(self):
        j_coupling = JCoupling(index1=1, index2=2, index3=3, index4=4)
        assert j_coupling.index1 == 1
        assert j_coupling.index2 == 2
        assert j_coupling.index3 == 3
        assert j_coupling.index4 == 4
        assert j_coupling.j_coupling_type is None
        assert j_coupling.comment is None

        j_coupling = JCoupling(
            index1=1,
            index2=2,
            index3=3,
            index4=4,
            j_coupling_type=1,
            comment="comment"
        )
        assert j_coupling.index1 == 1
        assert j_coupling.index2 == 2
        assert j_coupling.index3 == 3
        assert j_coupling.index4 == 4
        assert j_coupling.j_coupling_type == 1
        assert j_coupling.comment == "comment"

    def test_copy(self):
        j_coupling = JCoupling(
            index1=1, index2=2, index3=3, index4=4, j_coupling_type=1
        )
        j_coupling_copy = j_coupling.copy()
        assert j_coupling_copy.index1 == 1
        assert j_coupling_copy.index2 == 2
        assert j_coupling_copy.index3 == 3
        assert j_coupling_copy.index4 == 4
        assert j_coupling_copy.j_coupling_type == 1

        j_coupling_copy.index1 = 2
        j_coupling_copy.index2 = 3
        j_coupling_copy.index3 = 4
        j_coupling_copy.index4 = 5
        j_coupling_copy.j_coupling_type = 2

        assert j_coupling_copy.index1 != j_coupling.index1
        assert j_coupling_copy.index2 != j_coupling.index2
        assert j_coupling_copy.index3 != j_coupling.index3
        assert j_coupling_copy.index4 != j_coupling.index4
        assert j_coupling_copy.j_coupling_type != j_coupling.j_coupling_type

    def test__eq__(self):
        j_coupling = JCoupling(index1=1, index2=2, index3=3, index4=4)
        assert j_coupling == j_coupling

        j_coupling_copy = j_coupling.copy()
        assert j_coupling == j_coupling_copy

        j_coupling_copy.index1 = 2
        assert j_coupling != j_coupling_copy

        assert j_coupling != 1
        assert j_coupling != "a"
        assert j_coupling != None
