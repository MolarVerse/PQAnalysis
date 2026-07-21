"""
Generator for the Green-Kubo test trajectory.

An Ornstein-Uhlenbeck (OU / AR(1)) velocity process is integrated to
produce a position trajectory (.xyz) and the consistent velocity
trajectory (.vel) that share the same dynamics, so that the
Einstein-MSD and the Green-Kubo velocity-ACF self-diffusion
coefficients refer to the same trajectory and must agree.

Units
-----
* time step dt in ps
* correlation time tau in ps  => phi = exp(-dt / tau)
* per-component velocity std sigma in Angstrom / ps (dynamics unit)
* positions x_{n+1} = x_n + v_n * dt  in Angstrom
* the .vel file stores velocities in Angstrom / s, i.e. v_n * 1e12
  (1 Angstrom / ps = 1e12 Angstrom / s), reproducing the PQ velocity
  output magnitude of ~1e12 - 1e13 Angstrom / s.

Analytic diffusion coefficient (per component and total scalar):
    D = sigma^2 * tau        [Angstrom^2 / ps]
      = sigma^2 * tau * 1e-8  [m^2 / s]
"""

import numpy as np

# --- parameters ---------------------------------------------------------
SEED = 20240721
N_FRAMES = 3000
N_ATOMS = 64
DT_PS = 0.002            # ps
TAU_PS = 0.1             # ps
SIGMA_ANG_PER_PS = 2.0   # Angstrom / ps
BOX = 100.0              # Angstrom (large: no wrapping for tiny drift)
PS_TO_S = 1.0e12         # 1 Angstrom/ps = 1e12 Angstrom/s

ANGSTROM2_PER_PS_TO_M2_PER_S = 1.0e-8


def generate():
    rng = np.random.default_rng(SEED)

    phi = np.exp(-DT_PS / TAU_PS)
    noise_scale = SIGMA_ANG_PER_PS * np.sqrt(1.0 - phi * phi)

    # OU velocities in Angstrom / ps, shape (N_FRAMES, N_ATOMS, 3)
    vel = np.zeros((N_FRAMES, N_ATOMS, 3))
    vel[0] = SIGMA_ANG_PER_PS * rng.standard_normal((N_ATOMS, 3))
    for n in range(1, N_FRAMES):
        vel[n] = phi * vel[n - 1] + noise_scale * rng.standard_normal(
            (N_ATOMS, 3)
        )

    # positions in Angstrom: x_{n+1} = x_n + v_n * dt
    pos = np.zeros_like(vel)
    pos[0] = rng.uniform(-3.0, 3.0, (N_ATOMS, 3))
    for n in range(1, N_FRAMES):
        pos[n] = pos[n - 1] + vel[n - 1] * DT_PS

    # .vel velocities in Angstrom / s
    vel_file = vel * PS_TO_S

    _write_xyz("traj.xyz", pos)
    _write_vel("traj.vel", vel_file)

    d_analytic_m2_s = (
        SIGMA_ANG_PER_PS ** 2 * TAU_PS * ANGSTROM2_PER_PS_TO_M2_PER_S
    )

    print(f"phi = {phi:.10f}")
    print(f"velocity file magnitude (rms) = {np.sqrt(np.mean(vel_file**2)):.3e} Angstrom/s")
    print(f"D_analytic = {d_analytic_m2_s:.6e} m^2/s")
    print(f"D_analytic = {d_analytic_m2_s * 1e4:.6e} cm^2/s")


def _write_xyz(filename, pos):
    with open(filename, "w", encoding="utf-8") as file:
        for frame_index, frame in enumerate(pos):
            file.write(f"{N_ATOMS} {BOX} {BOX} {BOX}\n")
            file.write(f"frame {frame_index + 1}\n")
            for x, y, z in frame:
                file.write(f"Ar {x:14.7f} {y:14.7f} {z:14.7f}\n")


def _write_vel(filename, vel):
    with open(filename, "w", encoding="utf-8") as file:
        for frame in vel:
            file.write(f"{N_ATOMS} {BOX} {BOX} {BOX}\n\n")
            for x, y, z in frame:
                file.write(f"Ar {x:16.7e} {y:16.7e} {z:16.7e}\n")


if __name__ == "__main__":
    generate()
