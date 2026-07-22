"""
Numerical routines and file writers for vibrational analysis.
"""

from contextlib import contextmanager
from dataclasses import dataclass
import sys
from pathlib import Path

import numpy as np

from PQAnalysis.atomic_system import AtomicSystem

from .exceptions import VibrationalAnalysisError

LINEAR_ROTATION_RTOL = 1e-6
SPEED_OF_LIGHT_CM_S = 2.99792458e10
BOLTZMANN_EV_K = 8.617333262145e-5
WAVENUMBER_TO_EV = 1.2398419843320026e-4
MODE_THRESHOLD_CM = 1.0e-8



@contextmanager
def _output_stream(filename: str | None):
    """
    Open a named output file or yield stdout.
    """
    if filename is None:
        yield sys.stdout
        return

    with open(filename, "w", encoding="utf-8") as file:
        yield file



@dataclass(frozen=True)
class VibrationalAnalysisResult:

    """
    Result container for a vibrational analysis.
    """

    wavenumbers: np.ndarray
    force_constants: np.ndarray
    reduced_masses: np.ndarray
    normal_modes: np.ndarray
    intensities: np.ndarray | None = None



def read_hessian_file(filename: str) -> np.ndarray:
    """
    Read a plain square Hessian matrix.

    Parameters
    ----------
    filename : str
        The Hessian file.

    Returns
    -------
    np.ndarray
        The square Hessian matrix.
    """
    path = Path(filename)
    if not path.is_file():
        raise VibrationalAnalysisError(f"Hessian file '{filename}' not found.")

    try:
        hessian_data = np.loadtxt(path, dtype=float)
    except ValueError as exception:
        raise VibrationalAnalysisError(
            f"Hessian file '{filename}' contains non-numeric data."
        ) from exception

    hessian = np.atleast_2d(hessian_data)

    if hessian.shape[0] != hessian.shape[1]:
        raise VibrationalAnalysisError(
            f"Hessian matrix must be square, got shape {hessian.shape}."
        )

    return hessian



def calculate(
    atom_masses: np.ndarray,
    atom_coords: np.ndarray,
    hessian: np.ndarray,
    atom_charges: np.ndarray | None = None,
    unit: str = "kcal",
    hessian_sign: str | float = "auto",
) -> VibrationalAnalysisResult:
    """
    Calculate wavenumbers, force constants, reduced masses and normal modes.

    Parameters
    ----------
    atom_masses : np.ndarray
        Atomic masses in amu.
    atom_coords : np.ndarray
        Atomic coordinates with shape ``(n_atoms, 3)``.
    hessian : np.ndarray
        Cartesian Hessian with shape ``(3 * n_atoms, 3 * n_atoms)``.
    atom_charges : np.ndarray | None, optional
        Atomic charges used to calculate IR intensities, by default None.
    unit : str, optional
        Hessian energy unit. Supported values are ``kcal``, ``hartree`` and
        ``ev``, by default ``kcal``.
    hessian_sign : str | float, optional
        Hessian sign convention. Supported values are ``auto``, ``positive``,
        ``negative``, ``1`` and ``-1``, by default ``auto``.
    """
    atom_masses = np.asarray(atom_masses, dtype=float)
    atom_coords = np.asarray(atom_coords, dtype=float)
    hessian = np.asarray(hessian, dtype=float)

    _check_shapes(atom_masses, atom_coords, hessian)

    sign_factor = hessian_sign_factor(
        atom_coords,
        atom_masses,
        hessian,
        hessian_sign,
    )
    hessian_mw = mass_weighted_hessian(
        hessian,
        atom_masses,
        sign=sign_factor,
    )

    eigenvalues, normal_modes, normalization = internal_coordinates(
        atom_coords,
        atom_masses,
        hessian_mw,
    )

    wavenumbers, omega = wavenumber(eigenvalues, unit=unit)
    reduced_masses = reduced_mass(normalization)
    force_constants = force_constant(omega, reduced_masses)

    intensities = None
    if atom_charges is not None:
        atom_charges = np.asarray(atom_charges, dtype=float)
        if atom_charges.shape != atom_masses.shape:
            raise VibrationalAnalysisError(
                "The number of atom charges must match the number of atoms."
            )

        intensities = infrared_intensity(
            normal_modes,
            atom_charges,
            reduced_masses,
        )

    return VibrationalAnalysisResult(
        wavenumbers=wavenumbers,
        intensities=intensities,
        force_constants=force_constants,
        reduced_masses=reduced_masses,
        normal_modes=normal_modes,
    )



def calculate_from_system(
    system: AtomicSystem,
    hessian: np.ndarray,
    unit: str = "kcal",
    hessian_sign: str | float = "auto",
    atom_charges: np.ndarray | None = None,
) -> VibrationalAnalysisResult:
    """
    Calculate vibrational data from an AtomicSystem and Hessian.
    """
    return calculate(
        system.atomic_masses,
        system.pos,
        hessian,
        atom_charges=atom_charges,
        unit=unit,
        hessian_sign=hessian_sign,
    )



def _check_shapes(
    atom_masses: np.ndarray,
    atom_coords: np.ndarray,
    hessian: np.ndarray,
) -> None:
    """
    Validate array shapes before running the numerical analysis.
    """
    if atom_masses.ndim != 1:
        raise VibrationalAnalysisError(
            "Atom masses must be a one-dimensional array."
        )

    if atom_coords.shape != (atom_masses.size, 3):
        raise VibrationalAnalysisError(
            "Atom coordinates must have shape (n_atoms, 3)."
        )

    expected_hessian_shape = (3 * atom_masses.size, 3 * atom_masses.size)
    if hessian.shape != expected_hessian_shape:
        raise VibrationalAnalysisError(
            "Hessian shape must be "
            f"{expected_hessian_shape}, got {hessian.shape}."
        )



def center_to_com(
    atom_coords: np.ndarray, atom_masses: np.ndarray
) -> np.ndarray:
    """
    Translate coordinates to the center of mass.
    """
    center_of_mass = np.sum(atom_coords * atom_masses[:, None], axis=0)
    center_of_mass = center_of_mass / np.sum(atom_masses)

    return atom_coords - center_of_mass



def inertia_tensor(
    atom_coords: np.ndarray, atom_masses: np.ndarray
) -> np.ndarray:
    """
    Calculate the inertia tensor.
    """
    x = atom_coords[:, 0]
    y = atom_coords[:, 1]
    z = atom_coords[:, 2]

    i_xx = np.sum(atom_masses * (y**2 + z**2))
    i_yy = np.sum(atom_masses * (x**2 + z**2))
    i_zz = np.sum(atom_masses * (x**2 + y**2))
    i_xy = -np.sum(atom_masses * (x * y))
    i_xz = -np.sum(atom_masses * (x * z))
    i_yz = -np.sum(atom_masses * (y * z))

    return np.array(
        [
            [i_xx, i_xy, i_xz],
            [i_xy, i_yy, i_yz],
            [i_xz, i_yz, i_zz],
        ],
        dtype=float,
    )



def translational_modes(atom_masses: np.ndarray) -> np.ndarray:
    """
    Calculate normalized translational modes.
    """
    translation = np.zeros((3 * atom_masses.size, 3), dtype=float)

    for atom_index, mass in enumerate(atom_masses):
        start = 3 * atom_index
        translation[start:start + 3, :] = np.sqrt(mass) * np.eye(3)

    norms = np.sqrt(np.sum(translation**2, axis=0))
    return translation / norms



def rotational_modes(
    atom_coords: np.ndarray, atom_masses: np.ndarray
) -> np.ndarray:
    """
    Calculate normalized rotational modes.
    """
    atom_coords_cm = center_to_com(atom_coords, atom_masses)
    _, eigenvectors = np.linalg.eigh(inertia_tensor(atom_coords, atom_masses))

    x_frame = eigenvectors
    p_frame = atom_coords_cm @ eigenvectors
    sqrt_masses = np.sqrt(atom_masses)[:, None]

    rotation = np.zeros((atom_coords.shape[0] * 3, 3), dtype=float)

    rotation[:, 0] = (
        (
            p_frame[:, 1, None] * x_frame[2, :][None, :] -
            p_frame[:, 2, None] * x_frame[1, :][None, :]
        ) * sqrt_masses
    ).ravel(order="F")
    rotation[:, 1] = (
        (
            p_frame[:, 2, None] * x_frame[0, :][None, :] -
            p_frame[:, 0, None] * x_frame[2, :][None, :]
        ) * sqrt_masses
    ).ravel(order="F")
    rotation[:, 2] = (
        (
            p_frame[:, 0, None] * x_frame[1, :][None, :] -
            p_frame[:, 1, None] * x_frame[0, :][None, :]
        ) * sqrt_masses
    ).ravel(order="F")

    rotation_norms = np.sqrt(np.sum(rotation**2, axis=0))
    threshold = LINEAR_ROTATION_RTOL * max(1.0, float(np.max(rotation_norms)))
    keep = rotation_norms > threshold

    rotation = rotation[:, keep]

    if rotation.size:
        rotation = rotation / rotation_norms[keep]

    return rotation



def transformation_matrix(
    atom_coords: np.ndarray, atom_masses: np.ndarray
) -> np.ndarray:
    """
    Calculate the external-mode transformation matrix.
    """
    return np.hstack(
        [
            translational_modes(atom_masses),
            rotational_modes(atom_coords, atom_masses),
        ]
    )



def internal_subspace(
    atom_coords: np.ndarray, atom_masses: np.ndarray
) -> np.ndarray:
    """
    Calculate the internal-coordinate subspace.
    """
    transformation = transformation_matrix(atom_coords, atom_masses)
    total_modes = atom_coords.shape[0] * 3
    external_modes = transformation.shape[1]

    if external_modes >= total_modes:
        return np.zeros((total_modes, 0), dtype=float)

    q_matrix, _ = np.linalg.qr(transformation, mode="complete")
    return q_matrix[:, external_modes:total_modes]



def internal_coordinates(
    atom_coords: np.ndarray,
    atom_masses: np.ndarray,
    hessian_mw: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate normal modes from a mass-weighted Hessian.
    """
    transformation = transformation_matrix(atom_coords, atom_masses)
    transformation, _ = np.linalg.qr(transformation, mode="complete")

    internal_hessian = transformation.T @ hessian_mw @ transformation
    internal_hessian = symmetrize_addition(internal_hessian)

    eigenvalues, eigenvectors = np.linalg.eigh(internal_hessian)

    mass_factors = np.repeat(1.0 / np.sqrt(atom_masses), 3)
    mass_matrix = np.diag(mass_factors)
    eigenvectors_internal = mass_matrix @ transformation @ eigenvectors

    normalization = np.sqrt(1.0 / np.sum(eigenvectors_internal**2, axis=0))
    normal_modes = eigenvectors_internal * normalization[None, :]

    return eigenvalues, normal_modes, normalization



def symmetrize_addition(hessian: np.ndarray) -> np.ndarray:
    """
    Symmetrize a matrix by averaging it with its transpose.
    """
    return (hessian + hessian.T) / 2.0



def masses_matrix(atom_masses: np.ndarray) -> np.ndarray:
    """
    Create the mass-scaling matrix sqrt(m_i * m_j).
    """
    masses_repeat = np.repeat(atom_masses, 3)
    return np.sqrt(np.outer(masses_repeat, masses_repeat))



def mass_weighted_hessian(
    hessian: np.ndarray,
    atom_masses: np.ndarray,
    sign: float = 1.0,
) -> np.ndarray:
    """
    Symmetrize and mass-weight a Hessian.
    """
    hessian = symmetrize_addition(float(sign) * hessian)
    return hessian / masses_matrix(atom_masses)



def hessian_sign_factor(
    atom_coords: np.ndarray,
    atom_masses: np.ndarray,
    hessian: np.ndarray,
    hessian_sign: str | float,
) -> float:
    """
    Resolve a Hessian sign setting into a numeric factor.
    """
    if isinstance(hessian_sign, (int, float)):
        if hessian_sign in {1, -1}:
            return float(hessian_sign)
    else:
        value = hessian_sign.lower()
        if value == "positive":
            return 1.0
        if value == "negative":
            return -1.0
        if value == "auto":
            hessian_mw = mass_weighted_hessian(hessian, atom_masses, sign=1.0)
            subspace = internal_subspace(atom_coords, atom_masses)

            if subspace.shape[1] == 0:
                return 1.0

            eigenvalues = np.linalg.eigvalsh(
                subspace.T @ hessian_mw @ subspace
            )
            threshold = np.sqrt(np.finfo(float).eps)
            threshold *= max(1.0, float(np.max(np.abs(eigenvalues))))
            positive_count = np.count_nonzero(eigenvalues > threshold)
            negative_count = np.count_nonzero(eigenvalues < -threshold)

            if negative_count > positive_count:
                return -1.0
            if positive_count > negative_count:
                return 1.0

            positive_weight = np.sum(
                np.abs(eigenvalues[eigenvalues > threshold])
            )
            negative_weight = np.sum(
                np.abs(eigenvalues[eigenvalues < -threshold])
            )
            return -1.0 if negative_weight > positive_weight else 1.0

    raise VibrationalAnalysisError(
        "hessian_sign must be auto, positive, negative, 1, or -1."
    )



def wavenumber(eigenvalues: np.ndarray,
               unit: str = "kcal") -> tuple[np.ndarray, np.ndarray]:
    """
    Convert Hessian eigenvalues to wavenumbers.
    """
    unit = unit.lower()

    if unit == "kcal":
        factor = 4184.0 * 1.0e23
    elif unit == "hartree":
        factor = 2625500.2 * (0.188972598857892e11**2) * 1000.0
    elif unit == "ev":
        factor = 96485.307499 * 1.0e23
    else:
        raise VibrationalAnalysisError(
            "Invalid unit. Options are kcal, hartree and ev."
        )

    omega = signed_sqrt(eigenvalues * factor)
    wavenumbers = omega / (2.0 * np.pi * SPEED_OF_LIGHT_CM_S)

    return wavenumbers, omega



def signed_sqrt(values: np.ndarray) -> np.ndarray:
    """
    Calculate a sign-preserving square root.
    """
    return np.sign(values) * np.sqrt(np.abs(values))



def reduced_mass(normalization: np.ndarray) -> np.ndarray:
    """
    Calculate reduced masses from normal-mode normalization factors.
    """
    return normalization**2



def force_constant(
    omega: np.ndarray, reduced_masses: np.ndarray
) -> np.ndarray:
    """
    Calculate force constants in mdyn A^-1.
    """
    return omega**2 * reduced_masses / 6.022 / 1.0e28



def infrared_intensity(
    normal_modes: np.ndarray,
    atom_charges: np.ndarray,
    reduced_masses: np.ndarray,
) -> np.ndarray:
    """
    Calculate infrared intensities in km mol^-1.
    """
    intensities = np.zeros(normal_modes.shape[1], dtype=float)

    for index in range(normal_modes.shape[1]):
        mode = normal_modes[:, index].reshape((-1, 3))
        intensity = np.sum(
            (
                np.sum(mode * atom_charges[:, None], axis=0) / 0.2081943 /
                np.linalg.norm(mode)
            )**2 / reduced_masses[index] * 42.2561
        )
        intensities[index] = intensity

    return intensities



def write_calculate_output(
    result: VibrationalAnalysisResult,
    filename: str | None = None,
) -> None:
    """
    Write the tabular vibrational analysis output.
    """
    with _output_stream(filename) as file:
        if result.intensities is None:
            print(
                "# Wavenumbers (cm-1)  Force constants (mdyn A^-1)  "
                "Reduced masses (amu)",
                file=file,
            )
            for wavenumber_value, force_const, reduced_mass_value in zip(
                result.wavenumbers,
                result.force_constants,
                result.reduced_masses,
            ):
                print(
                    f"{wavenumber_value:<8.8e}\t"
                    f"{force_const:<8.8e}\t"
                    f"{reduced_mass_value:<8.8e}",
                    file=file,
                )
        else:
            print(
                "# Wavenumbers (cm-1)  Intensities (km mol-1)  "
                "Force constants (mdyn A^-1)  Reduced masses (amu)",
                file=file,
            )
            for (
                wavenumber_value,
                intensity,
                force_const,
                reduced_mass_value,
            ) in zip(
                result.wavenumbers,
                result.intensities,
                result.force_constants,
                result.reduced_masses,
            ):
                print(
                    f"{wavenumber_value:<8.8e}\t"
                    f"{intensity:<8.8e}\t"
                    f"{force_const:<8.8e}\t"
                    f"{reduced_mass_value:<8.8e}",
                    file=file,
                )



def write_normal_modes(
    normal_modes: np.ndarray,
    filename: str | None = None,
) -> None:
    """
    Write normal modes in matrix form.
    """
    with _output_stream(filename) as file:
        for row in normal_modes:
            print(" ".join(str(value) for value in row), file=file)



def select_mode_indices(
    wavenumbers: np.ndarray,
    modes: str | list[int] | None = "all",
    threshold: float = MODE_THRESHOLD_CM,
) -> list[int]:
    """
    Select mode indices from a user-facing one-based mode selection.
    """
    wavenumbers = np.asarray(wavenumbers, dtype=float)
    n_modes = wavenumbers.size

    if modes is None or modes == "all":
        return list(range(n_modes))

    if isinstance(modes, str):
        modes = modes.lower()
        if modes == "nonzero":
            return [
                index for index, value in enumerate(wavenumbers)
                if abs(value) > threshold
            ]
        if modes == "positive":
            return [
                index for index, value in enumerate(wavenumbers)
                if value > threshold
            ]

        raise VibrationalAnalysisError(
            "modes must be all, nonzero, positive or one-based mode numbers."
        )

    selected_indices = []
    for mode in modes:
        if mode < 1 or mode > n_modes:
            raise VibrationalAnalysisError(
                f"Mode {mode} is outside the available range 1..{n_modes}."
            )

        mode_index = mode - 1
        if mode_index not in selected_indices:
            selected_indices.append(mode_index)

    return selected_indices



def write_xyz_modes(
    normal_modes: np.ndarray,
    atom_coords: np.ndarray,
    atom_names: list[str],
    filename: str = "modes",
    amplitude: float = 0.25,
    step: float | None = None,
    *,
    wavenumbers: np.ndarray | None = None,
    modes: str | list[int] | None = "all",
    n_frames: int | None = None,
    temperature: float | None = None,
    threshold: float = MODE_THRESHOLD_CM,
) -> None:
    """
    Write one sinusoidal XYZ trajectory per selected normal mode.
    """
    n_frames = _mode_frame_count(n_frames, amplitude, step)
    _check_mode_writer_settings(n_frames, amplitude, temperature, threshold)

    n_atoms = atom_coords.shape[0]
    wavenumbers = _mode_wavenumbers(wavenumbers, normal_modes)
    mode_indices = select_mode_indices(wavenumbers, modes, threshold)
    phases = np.linspace(0.0, 2.0 * np.pi, n_frames, endpoint=False)

    for mode_index in mode_indices:
        mode = normal_modes[:, mode_index].reshape((-1, 3))
        displacement = mode_displacement(
            mode,
            wavenumbers[mode_index],
            amplitude=amplitude,
            temperature=temperature,
            threshold=threshold,
        )

        with open(
            f"{filename}-{mode_index + 1}.xyz", "w", encoding="utf-8"
        ) as file:
            for frame, phase in enumerate(phases, start=1):
                print(n_atoms, file=file)
                print(
                    _mode_comment(
                        mode_index,
                        wavenumbers[mode_index],
                        frame,
                        n_frames,
                        phase,
                    ),
                    file=file,
                )

                for atom_name, atom_coord, atom_mode in zip(
                    atom_names,
                    atom_coords,
                    displacement,
                ):
                    coords = atom_coord + np.sin(phase) * atom_mode
                    symbol = _xyz_atom_symbol(atom_name)
                    print(
                        f"{symbol}    {coords[0]}   {coords[1]}   {coords[2]}",
                        file=file,
                    )



def write_extxyz_modes(
    normal_modes: np.ndarray,
    atom_coords: np.ndarray,
    atom_names: list[str],
    filename: str = "modes.xyz",
    wavenumbers: np.ndarray | None = None,
    intensities: np.ndarray | None = None,
    modes: str | list[int] | None = "all",
    threshold: float = MODE_THRESHOLD_CM,
) -> None:
    """
    Write selected normal modes to one extended XYZ file.
    """
    if threshold < 0.0:
        raise VibrationalAnalysisError("Mode threshold must not be negative.")

    n_atoms = atom_coords.shape[0]
    wavenumbers = _mode_wavenumbers(wavenumbers, normal_modes)
    mode_indices = select_mode_indices(wavenumbers, modes, threshold)

    with open(filename, "w", encoding="utf-8") as file:
        for mode_index in mode_indices:
            mode = normal_modes[:, mode_index].reshape((-1, 3))

            print(n_atoms, file=file)
            print(
                _extxyz_comment(
                    mode_index,
                    wavenumbers[mode_index],
                    None if intensities is None else intensities[mode_index],
                ),
                file=file,
            )

            for atom_name, atom_coord, atom_mode in zip(
                atom_names,
                atom_coords,
                mode,
            ):
                symbol = _xyz_atom_symbol(atom_name)
                print(
                    (
                        f"{symbol} "
                        f"{atom_coord[0]:.12e} {atom_coord[1]:.12e} "
                        f"{atom_coord[2]:.12e} {atom_mode[0]:.12e} "
                        f"{atom_mode[1]:.12e} {atom_mode[2]:.12e}"
                    ),
                    file=file,
                )



def mode_displacement(
    mode: np.ndarray,
    wavenumber_value: float,
    amplitude: float = 0.25,
    temperature: float | None = None,
    threshold: float = MODE_THRESHOLD_CM,
) -> np.ndarray:
    """
    Scale a mode for display.
    """
    if temperature is None:
        return _scale_mode_to_amplitude(mode, amplitude)

    if temperature <= 0.0:
        raise VibrationalAnalysisError("Mode temperature must be positive.")

    energy = abs(wavenumber_value) * WAVENUMBER_TO_EV
    if energy <= threshold * WAVENUMBER_TO_EV:
        return _scale_mode_to_amplitude(mode, amplitude)

    return mode * np.sqrt(BOLTZMANN_EV_K * temperature / energy)



def _scale_mode_to_amplitude(mode: np.ndarray, amplitude: float) -> np.ndarray:
    """
    Scale a mode so the largest atomic displacement equals amplitude.
    """
    if amplitude < 0.0:
        raise VibrationalAnalysisError("Mode amplitude must not be negative.")

    if amplitude == 0.0:
        return np.zeros_like(mode)

    max_displacement = np.max(np.linalg.norm(mode, axis=1))
    if max_displacement == 0.0:
        return np.zeros_like(mode)

    return mode * amplitude / max_displacement



def _mode_wavenumbers(
    wavenumbers: np.ndarray | None,
    normal_modes: np.ndarray,
) -> np.ndarray:
    """
    Return mode wavenumbers or a zero placeholder array.
    """
    if wavenumbers is None:
        return np.zeros(normal_modes.shape[1], dtype=float)

    wavenumbers = np.asarray(wavenumbers, dtype=float)
    if wavenumbers.shape != (normal_modes.shape[1], ):
        raise VibrationalAnalysisError(
            "The number of wavenumbers must match the number of normal modes."
        )

    return wavenumbers



def _check_mode_writer_settings(
    n_frames: int,
    amplitude: float,
    temperature: float | None,
    threshold: float,
) -> None:
    """
    Validate mode writer settings.
    """
    if n_frames < 1:
        raise VibrationalAnalysisError(
            "Number of mode frames must be positive."
        )
    if amplitude < 0.0:
        raise VibrationalAnalysisError("Mode amplitude must not be negative.")
    if temperature is not None and temperature <= 0.0:
        raise VibrationalAnalysisError("Mode temperature must be positive.")
    if threshold < 0.0:
        raise VibrationalAnalysisError("Mode threshold must not be negative.")



def _mode_frame_count(
    n_frames: int | None,
    amplitude: float,
    step: float | None,
) -> int:
    """
    Return an explicit frame count, preserving the legacy step setting.
    """
    if n_frames is not None:
        return n_frames

    if step is None:
        return 30

    if step <= 0.0:
        raise VibrationalAnalysisError("Mode step must be positive.")

    return max(1, int(np.floor(2.0 * amplitude / step)) + 1)



def _mode_comment(
    mode_index: int,
    wavenumber_value: float,
    frame: int,
    n_frames: int,
    phase: float,
) -> str:
    """
    Build a standard XYZ comment line for one animated mode frame.
    """
    return (
        f"mode={mode_index + 1} "
        f"frequency_cm-1={wavenumber_value:.8e} "
        f"frame={frame}/{n_frames} phase={phase:.8e}"
    )



def _extxyz_comment(
    mode_index: int,
    wavenumber_value: float,
    intensity: float | None,
) -> str:
    """
    Build an extended XYZ comment line for one mode vector image.
    """
    comment = (
        "Properties=species:S:1:pos:R:3:mode:R:3 "
        f"mode={mode_index + 1} frequency_cm-1={wavenumber_value:.8e}"
    )

    if intensity is not None:
        comment += f" IR_intensity={intensity:.8e}"

    return comment



def _xyz_atom_symbol(atom_name: str | None) -> str:
    """
    Return an XYZ/Jmol-friendly atom symbol.
    """
    if atom_name is None:
        return "X"

    atom_name = str(atom_name)
    if not atom_name:
        return "X"

    return atom_name[0].upper() + atom_name[1:].lower()
