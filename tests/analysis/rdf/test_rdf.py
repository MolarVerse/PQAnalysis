# import numpy as np
# import pytest

# from PQAnalysis.analysis.rdf.rdf import _calculate_n_bins, _infer_r_max, _check_r_max, _calculate_r_max, _integration, _norm, _setup_bin_middle_points, _add_to_bins
# from PQAnalysis.analysis.rdf.exceptions import RDFError, RDFWarning
# from PQAnalysis.analysis import RDF
# from PQAnalysis.traj import Trajectory, Frame
# from PQAnalysis.core import AtomicSystem, Cell


# def test__calculate_n_bins():
#     r_min = 1.0
#     r_max = 101.5
#     delta_r = 1.0

#     n_bins, r_max = _calculate_n_bins(r_min, r_max, delta_r)

#     assert n_bins == 100
#     assert np.isclose(r_max, 101.0)


# def test__infer_r_max():

#     system1 = AtomicSystem(cell=Cell(10, 10, 10, 90, 90, 90))
#     system2 = AtomicSystem(cell=Cell(16, 13, 12, 90, 90, 90))

#     traj = Trajectory([Frame(system1), Frame(system2)])

#     r_max = _infer_r_max(traj)

#     assert np.isclose(r_max, 5.0)

#     system3 = AtomicSystem()
#     traj.append(Frame(system3))

#     print(traj.check_vacuum())

#     with pytest.raises(RDFError) as exception:
#         r_max = _infer_r_max(traj)
#     assert str(exception.value) == "To infer r_max of the RDF analysis, the trajectory cannot be a vacuum trajectory. Please specify r_max manually or use the combination n_bins and delta_r."


# def test__check_r_max():
#     r_max = 5.0
#     traj = Trajectory()

#     assert np.isclose(_check_r_max(r_max, traj), r_max)

#     system1 = AtomicSystem(cell=Cell(10, 10, 10, 90, 90, 90))
#     system2 = AtomicSystem(cell=Cell(16, 13, 12, 90, 90, 90))

#     traj = Trajectory([Frame(system1), Frame(system2)])

#     assert np.isclose(_check_r_max(r_max, traj), r_max)

#     r_max = 10.0

#     warning_message = f"The calculated r_max {r_max} is larger than the maximum allowed radius \
# according to the box vectors of the trajectory 5.0. \
# r_max will be set to the maximum allowed radius."

#     with pytest.warns(RDFWarning, match=warning_message):
#         r_max = _check_r_max(r_max, traj)
#         assert np.isclose(r_max, 5.0)


# def test__calculate_r_max():
#     n_bins = 50
#     delta_r = 0.1
#     r_min = 0.0
#     traj = Trajectory()

#     r_max = _calculate_r_max(n_bins, delta_r, r_min, traj)

#     assert np.isclose(r_max, 5.0)

#     r_min = 3.0
#     r_max = _calculate_r_max(n_bins, delta_r, r_min, traj)

#     assert np.isclose(r_max, 8.0)

#     system1 = AtomicSystem(cell=Cell(10, 10, 10, 90, 90, 90))
#     system2 = AtomicSystem(cell=Cell(16, 13, 12, 90, 90, 90))

#     traj = Trajectory([Frame(system1), Frame(system2)])

#     warning_message = f"The calculated r_max {r_max} is larger than the maximum allowed radius \
# according to the box vectors of the trajectory 5.0. \
# r_max will be set to the maximum allowed radius."

#     print(warning_message)

#     with pytest.warns(RDFWarning, match=warning_message):
#         r_max = _calculate_r_max(n_bins, delta_r, r_min, traj)
#         assert np.isclose(r_max, 5.0)


# def test__setup_bin_middle_points():
#     n_bins = 5
#     r_min = 3.0
#     r_max = 8.0
#     delta_r = 1.0

#     bin_middle_points = _setup_bin_middle_points(n_bins, r_min, r_max, delta_r)

#     assert np.allclose(bin_middle_points, np.array([3.5, 4.5, 5.5, 6.5, 7.5]))


# def test__integration():
#     bins = np.array([1, 2, 3, 4, 5])
#     len_reference_indices = 3
#     len_frames = 10

#     integration = _integration(bins, len_reference_indices, len_frames)

#     n_total = len_reference_indices * len_frames
#     assert np.allclose(integration, np.array(
#         [1 / n_total, 3 / n_total, 6 / n_total, 10 / n_total, 15 / n_total]))


# def test__norm():
#     n_bins = 5
#     n_frames = 10
#     n_reference_indices = 3
#     delta_r = 1.0
#     target_density = 2.0

#     norm = _norm(n_bins, delta_r, target_density,
#                  n_reference_indices, n_frames)

#     help_1 = np.arange(0, n_bins)
#     help_2 = np.arange(1, n_bins + 1)
#     norm_ref = (help_2**3 - help_1**3)*delta_r**3 * 4 / 3 * np.pi

#     assert np.allclose(norm, norm_ref * target_density *
#                        n_reference_indices * n_frames)


# def test__add_to_bins():
#     n_bins = 5
#     r_min = 3.0
#     delta_r = 1.0

#     distances = np.array([1.5, 2.5, 3.5, 3.6, 3.7, 4.5, 4.6, 5.5, 6.5, 8.5])

#     assert np.allclose(_add_to_bins(distances, r_min, delta_r,
#                        n_bins), np.array([3, 2, 1, 1, 0]))


# class TestRDF:
#     def test__init__(self):
#         with pytest.raises(RDFError) as exception:
#             RDF(Trajectory(), ["h"], [
#                 "h"], r_max=8.0, r_min=3.0)
#         assert str(
#             exception.value) == "Trajectory cannot be of length 0."

#         system1 = AtomicSystem(cell=Cell(10, 10, 10, 90, 90, 90))
#         system2 = AtomicSystem(cell=Cell())
#         traj = Trajectory([Frame(system1), Frame(system2)])

#         with pytest.raises(RDFError) as exception:
#             RDF(
#                 traj, ["h"], ["h"], r_max=8.0, r_min=3.0)
#         assert str(exception.value) == "The provided trajectory is not fully periodic or in vacuum, meaning that some frames are in vacuum and others are periodic. This is not supported by the RDF analysis."

#         system1 = AtomicSystem(cell=Cell(10, 10, 10, 90, 90, 90))
#         system2 = AtomicSystem(cell=Cell(16, 13, 12, 90, 90, 90))

#         traj = Trajectory([Frame(system1), Frame(system2)])

#         with pytest.raises(RDFError) as exception:
#             RDF(
#                 traj, ["h"], ["h"], r_max=8.0, r_min=3.0)
#         assert str(
#             exception.value) == "Either n_bins or delta_r must be specified."

#         with pytest.raises(RDFError) as exception:
#             RDF(
#                 traj, ["h"], ["h"], r_max=8.0, r_min=3.0, delta_r=0.1, n_bins=5)
#         assert str(exception.value) == "It is not possible to specify all of n_bins, delta_r and r_max in the same RDF analysis as this would lead to ambiguous results."

#         # initialize rdf only with n_bins and delta_r

#         n_bins = 5
#         delta_r = 1.0

#         rdf = RDF(
#             traj, ["h"], ["h"], delta_r=delta_r, n_bins=n_bins)

#         assert np.isclose(rdf.r_max, 5.0)
#         assert np.isclose(rdf.r_min, 0.0)
#         assert len(rdf.bins) == 5
#         assert np.allclose(rdf.bin_middle_points,
#                            np.array([0.5, 1.5, 2.5, 3.5, 4.5]))
#         assert rdf.n_bins == 5
#         assert np.isclose(rdf.delta_r, 1.0)

#         # r_max has to be taken from trajectory

#         n_bins = 10

#         with pytest.warns(RDFWarning):
#             rdf = RDF(
#                 traj, ["h"], ["h"], delta_r=delta_r, n_bins=n_bins)

#         assert np.isclose(rdf.r_max, 5.0)

#         system1 = AtomicSystem(cell=Cell())
#         system2 = AtomicSystem(cell=Cell())

#         traj = Trajectory([Frame(system1), Frame(system2)])

#         with pytest.raises(RDFError) as exception:
#             RDF(
#                 traj, ["h"], ["h"], delta_r=delta_r)
#         assert str(exception.value) == "To infer r_max of the RDF analysis, the trajectory cannot be a vacuum trajectory. Please specify r_max manually or use the combination n_bins and delta_r."

#         with pytest.raises(RDFError) as exception:
#             RDF(
#                 traj, ["h"], ["h"], n_bins=n_bins)
#         assert str(exception.value) == "To infer r_max of the RDF analysis, the trajectory cannot be a vacuum trajectory. Please specify r_max manually or use the combination n_bins and delta_r."

#         r_max = 5.0

#         rdf = RDF(
#             traj, ["h"], ["h"], delta_r=delta_r, r_max=r_max)

#         assert np.isclose(rdf.r_max, 5.0)
#         assert np.isclose(rdf.r_min, 0.0)
#         assert len(rdf.bins) == 5
#         assert np.allclose(rdf.bin_middle_points,
#                            np.array([0.5, 1.5, 2.5, 3.5, 4.5]))
#         assert rdf.n_bins == 5
#         assert np.isclose(rdf.delta_r, 1.0)

#         n_bins = 5

#         rdf = RDF(
#             traj, ["h"], ["h"], n_bins=n_bins, r_max=r_max)

#         assert np.isclose(rdf.r_max, 5.0)
#         assert np.isclose(rdf.r_min, 0.0)
#         assert len(rdf.bins) == 5
#         assert np.allclose(rdf.bin_middle_points,
#                            np.array([0.5, 1.5, 2.5, 3.5, 4.5]))
#         assert rdf.n_bins == 5
#         assert np.isclose(rdf.delta_r, 1.0)
