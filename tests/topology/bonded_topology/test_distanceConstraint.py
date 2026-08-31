import numpy as np

from .. import pytestmark

from PQAnalysis.topology import DistanceConstraint



class TestDistanceConstraint:

    def test__init__(self):
        constraint = DistanceConstraint(
            index1=1,
            index2=2,
            lower_distance=1.0,
            upper_distance=2.0,
            spring_constant=100.0,
            d_spring_constant_dt=-0.1
        )
        assert constraint.index1 == 1
        assert constraint.index2 == 2
        assert np.isclose(constraint.lower_distance, 1.0)
        assert np.isclose(constraint.upper_distance, 2.0)
        assert np.isclose(constraint.spring_constant, 100.0)
        assert np.isclose(constraint.d_spring_constant_dt, -0.1)
        assert constraint.comment is None

    def test_copy(self):
        constraint = DistanceConstraint(
            index1=1,
            index2=2,
            lower_distance=1.0,
            upper_distance=2.0,
            spring_constant=100.0,
            d_spring_constant_dt=-0.1,
            comment="comment"
        )
        constraint_copy = constraint.copy()
        assert constraint_copy.index1 == 1
        assert constraint_copy.index2 == 2
        assert np.isclose(constraint_copy.lower_distance, 1.0)
        assert np.isclose(constraint_copy.upper_distance, 2.0)
        assert np.isclose(constraint_copy.spring_constant, 100.0)
        assert np.isclose(constraint_copy.d_spring_constant_dt, -0.1)
        assert constraint_copy.comment == "comment"

        constraint_copy.index1 = 2
        constraint_copy.index2 = 3
        constraint_copy.lower_distance = 2.0
        constraint_copy.upper_distance = 3.0
        constraint_copy.spring_constant = 200.0
        constraint_copy.d_spring_constant_dt = -0.2

        assert constraint_copy.index1 != constraint.index1
        assert constraint_copy.index2 != constraint.index2
        assert constraint_copy.lower_distance != constraint.lower_distance
        assert constraint_copy.upper_distance != constraint.upper_distance
        assert constraint_copy.spring_constant != constraint.spring_constant
        assert (
            constraint_copy.d_spring_constant_dt !=
            constraint.d_spring_constant_dt
        )

    def test__eq__(self):
        constraint = DistanceConstraint(
            index1=1,
            index2=2,
            lower_distance=1.0,
            upper_distance=2.0,
            spring_constant=100.0,
            d_spring_constant_dt=-0.1
        )
        assert constraint == constraint

        constraint_copy = constraint.copy()
        assert constraint == constraint_copy

        constraint_copy.index1 = 2
        assert constraint != constraint_copy

        assert constraint != 1
        assert constraint != "a"
        assert constraint != None
