import pytest

from PQAnalysis.io.traj_file import TrajectoryReader



@pytest.mark.benchmark(group="TrajectoryReader")
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

    def benchmark_windows(self, benchmark):
        traj_reader = TrajectoryReader(
            ["tests/data/traj2box/test.xyz" for _ in range(5)]
        )

        def read_windows():
            for window in traj_reader.window_generator(window_size=1):
                pass

        benchmark(read_windows)

    def benchmark_window_half(self, benchmark):
        traj_reader = TrajectoryReader(
            ["tests/data/traj2box/test.xyz" for _ in range(5)]
        )

        def read_windows():
            for window in traj_reader.window_generator(
                window_size=1, trajectory_start=2, trajectory_stop=7
            ):  # 5 frames
                pass

        benchmark(read_windows)

    def benchmark_frames(self, benchmark):
        traj_reader = TrajectoryReader(
            ["tests/data/traj2box/test.xyz" for _ in range(5)]
        )

        def read_frames():
            for frame in traj_reader.frame_generator():
                pass

        benchmark(read_frames)

    def benchmark_frame_half(self, benchmark):
        traj_reader = TrajectoryReader(
            ["tests/data/traj2box/test.xyz" for _ in range(5)]
        )

        def read_frames():
            for frame in traj_reader.frame_generator(
                trajectory_start=2, trajectory_stop=7
            ):
                pass

        benchmark(read_frames)
