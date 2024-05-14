from PQAnalysis.io.traj_file import TrajectoryReader

class BenchmarkTrajReader:

    def benchmark_read_2frames(self, benchmark):
        traj_reader = TrajectoryReader(
            "tests/data/traj2box/test.xyz"
        )  # 2 frames
        benchmark(traj_reader.read)

    def benchmark_read_10frames(self, benchmark):
        traj_reader = TrajectoryReader(
            ["tests/data/traj2box/test.xyz" for _ in range(5)]
        )  # 10 frames
        benchmark(traj_reader.read)
