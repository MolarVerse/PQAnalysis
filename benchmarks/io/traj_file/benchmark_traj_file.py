import sys

from pympler import asizeof
from PQAnalysis.io import TrajectoryReader



def main():

    filename = sys.argv[1]

    reader = TrajectoryReader(filename)
    traj = reader.read()
    atomic_system = traj[0]

    traj2 = reader.read()

    print(asizeof.asized([traj, traj2], stats=1))

    print(asizeof.asized(traj))
    print(asizeof.asized(traj.logger))
    print(asizeof.asized(atomic_system))
    print(asizeof.asized(atomic_system.topology))
    print(asizeof.asized(atomic_system[0:0].topology))
    print(asizeof.asized(atomic_system[0].topology))
    print(asizeof.asized(atomic_system[0:2].topology))
    print(asizeof.asized(atomic_system[0:3].topology.atoms[1]))



if __name__ == "__main__":
    main()
