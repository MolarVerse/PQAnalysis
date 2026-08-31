import pytest

from PQAnalysis.io.traj_file import RawTrajectoryReader, TrajectoryReader



def _write_raw_trajectory(path, n_frames=2000, n_atoms=100):
    lines = []

    for frame_index in range(n_frames):
        lines.append(f"{n_atoms} 20.0 20.0 20.0\n\n")

        for atom_index in range(n_atoms):
            value = frame_index * 0.0001 + atom_index * 0.001
            lines.append(
                f"X {value:.12f} {value + 0.25:.12f} "
                f"{value - 0.5:.12f}\n"
            )

    path.write_text("".join(lines), encoding="utf-8")

    return str(path)



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



@pytest.mark.benchmark(group="RawTrajectoryReader")
class BenchmarkRawTrajReader:

    @staticmethod
    def _consume(reader):
        checksum = 0.0

        for values, _ in reader.raw_frame_generator():
            checksum += values[0, 0]

        return checksum

    def benchmark_float32(self, benchmark, tmp_path):
        filename = _write_raw_trajectory(tmp_path / "raw-f32.xyz")
        reader = RawTrajectoryReader(filename, dtype="float32")

        benchmark(self._consume, reader)

    def benchmark_float64(self, benchmark, tmp_path):
        filename = _write_raw_trajectory(tmp_path / "raw-f64.xyz")
        reader = RawTrajectoryReader(filename, dtype="float64")

        benchmark(self._consume, reader)

    def benchmark_float64_batch(self, benchmark, tmp_path):
        filename = _write_raw_trajectory(tmp_path / "raw-batch-f64.xyz")
        reader = RawTrajectoryReader(filename, dtype="float64")

        def read_batch():
            batch = reader.try_read_all_frames(
                expected_n_atoms=100,
                expected_n_frames=2000,
                max_bytes=512 * 1024 * 1024,
            )
            assert batch is not None
            values, _cells = batch
            return values[0, 0, 0]

        benchmark(read_batch)
