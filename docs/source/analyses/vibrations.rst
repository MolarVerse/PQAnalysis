Vibrational Analysis
====================

Vibrational analysis diagonalizes the mass-weighted Cartesian Hessian of a
single structure [Wilson1955]_. Its eigenvectors are the harmonic normal modes;
its eigenvalues give signed wavenumbers, force constants and reduced masses.
When partial charges are supplied, point-charge infrared intensities are
reported as well. The whole calculation is the harmonic approximation applied
to one isolated structure: no dynamics, temperature or anharmonicity enters it.

Throughout this page :math:`N` is the number of atoms, :math:`i` indexes atoms,
:math:`a,b\in\{x,y,z\}` index Cartesian directions, and
:math:`\alpha,\beta\in\{1,\dots,3N\}` index Cartesian coordinates in the file
order :math:`x_1,y_1,z_1,x_2,\dots`. :math:`m_\alpha` is the mass in amu of the
atom owning coordinate :math:`\alpha`, and :math:`j` indexes normal modes.

Mass-weighted Hessian
---------------------

The Hessian read from ``hessian_file`` is symmetrized, multiplied by the sign
factor :math:`s\in\{+1,-1\}` and mass-weighted,

.. math::

   H^{\mathrm{mw}}_{\alpha\beta}
   = \frac{s}{2}\,
     \frac{H_{\alpha\beta}+H_{\beta\alpha}}{\sqrt{m_\alpha m_\beta}} .

Symmetrization is unconditional: a Hessian that is only approximately symmetric
is averaged with its transpose rather than rejected.

Normal-mode eigenproblem
------------------------

PQAnalysis first assembles a trial matrix :math:`D` of external modes. Its
translational columns are

.. math::

   D^{\mathrm{trans}}_{(i,a),b}
   = \frac{\sqrt{m_i}\,\delta_{ab}}{\sqrt{\sum_k m_k}} ,

and its rotational columns are mass-weighted rigid-body rotations, built from
the atomic positions relative to the center of mass and projected onto the
eigenvectors of the inertia tensor of the uncentered coordinates. Any
complete orthogonal basis leaves the spectrum unchanged, so this choice does
not affect the wavenumbers. A rotational column whose norm falls below
:math:`10^{-6}` times the larger of one and the biggest column norm is
discarded, which is what leaves a linear molecule with two rotations instead of
three.

A complete QR factorization of :math:`D` supplies an orthonormal basis
:math:`Q\in\mathbb{R}^{3N\times 3N}` whose leading columns span these external
modes. The symmetrized transformed matrix is then diagonalized:

.. math::

   \bigl(Q^{\mathsf{T}}H^{\mathrm{mw}}Q\bigr)\,\mathbf{c}_j
   = \lambda_j\,\mathbf{c}_j ,
   \qquad
   \lambda_1\le\lambda_2\le\dots\le\lambda_{3N} .

Because :math:`Q` is orthogonal, this is a similarity transformation and
:math:`\{\lambda_j\}` is exactly the spectrum of :math:`H^{\mathrm{mw}}`. The
external-mode basis orients the eigenvectors and defines the internal block
used by the sign heuristic below; it does not project translations and
rotations out of the reported spectrum. Those appear as the near-zero
eigenvalues of the full :math:`3N`-dimensional problem and are removed only by
the ``modes`` selection when mode files are written.

The eigenvectors are transformed back to Cartesian displacements and
normalized,

.. math::

   l_{\alpha j} = \frac{(Q\mathbf{c}_j)_\alpha}{\sqrt{m_\alpha}} ,
   \qquad
   e_{\alpha j}
   = \frac{l_{\alpha j}}{\bigl(\sum_\beta l_{\beta j}^{2}\bigr)^{1/2}} ,
   \qquad
   \sum_\alpha e_{\alpha j}^{2} = 1 ,

so the columns :math:`e_{\alpha j}` written to ``normal_modes_file`` are
dimensionless unit vectors. Modes are reported in order of increasing
:math:`\lambda_j`, so imaginary modes come first and the stiffest internal mode
comes last.

Eigenvalues, wavenumbers and the unit chain
-------------------------------------------

Each eigenvalue is converted to an angular frequency and then to a wavenumber,

.. math::

   \omega_j = \operatorname{sgn}(\lambda_j)\sqrt{\lvert\lambda_j\rvert\,C_u} ,
   \qquad
   \tilde{\nu}_j = \frac{\omega_j}{2\pi c} ,
   \qquad
   c = 2.99792458\times10^{10}\ \mathrm{cm\,s^{-1}} ,

with :math:`\omega_j` in rad·s⁻¹ and :math:`\tilde{\nu}_j` in cm⁻¹. The
constant :math:`C_u` is fixed by the ``unit`` key and carries the entire unit
chain from the Hessian file to s⁻²:

.. list-table:: Conversion constants selected by ``unit``
   :class: pq-record-table
   :header-rows: 1
   :widths: 14 30 56

   * - ``unit``
     - Expected Hessian unit
     - Factor :math:`C_u` taking :math:`\lambda_j` to s⁻²
   * - ``kcal``
     - kcal·mol⁻¹·Å⁻²
     - :math:`4184\times10^{23}`
   * - ``ev``
     - eV·Å⁻²
     - :math:`96485.307499\times10^{23}`
   * - ``hartree``
     - hartree·bohr⁻²
     - :math:`2\,625\,500.2\times(1.88972598857892\times10^{10})^{2}\times10^{3}`

Each constant is the product of three conversions. For ``kcal``, the factor
:math:`4184` converts kcal to J, the factor :math:`10^{20}` converts Å⁻² to
m⁻², and the remaining :math:`10^{3}` comes from combining the molar energy
with the reciprocal atomic mass unit: because
:math:`1\ \mathrm{amu} = 10^{-3}\ \mathrm{kg\,mol^{-1}}/N_\mathrm{A}`, the
Avogadro constant cancels and leaves :math:`10^{3}\ \mathrm{kg^{-1}}`. That
cancellation is why no value of :math:`N_\mathrm{A}` appears anywhere in the
conversion. ``ev`` follows the same chain with
:math:`96485.307499\ \mathrm{J\,mol^{-1}}` per eV. ``hartree`` replaces the
Å⁻² step with
:math:`(1/a_0)^{2} = (1.88972598857892\times10^{10}\ \mathrm{m^{-1}})^{2}`, so
this option expects the Hessian in bohr⁻², while ``kcal`` and ``ev`` expect
Å⁻². Choosing the energy unit therefore also chooses the length unit.

Imaginary modes and the ``hessian_sign`` heuristic
--------------------------------------------------

The square root is taken with a sign-preserving convention,
:math:`\operatorname{sgn}(x)\sqrt{\lvert x\rvert}`, so a negative eigenvalue
produces a negative :math:`\omega_j` and a negative wavenumber rather than an
imaginary number. A reported :math:`\tilde{\nu}_j < 0` therefore means an
imaginary mode of magnitude :math:`\lvert\tilde{\nu}_j\rvert` cm⁻¹, that is,
negative curvature of the potential-energy surface along that mode.

The sign convention of the input file matters because Hessians are written
either as second derivatives of the energy or as derivatives of the forces,
which differ by a factor of :math:`-1`. In an input file ``hessian_sign``
accepts ``positive`` (:math:`s=+1`), ``negative`` (:math:`s=-1`) and ``auto``.
The Python interface additionally accepts the numbers ``1`` and ``-1``.

``auto`` resolves the convention from the curvature statistics of the internal
subspace. Let :math:`U` be the block of :math:`Q` spanning the complement of
the external modes, and let :math:`\lambda^{\mathrm{int}}` be the eigenvalues
of :math:`U^{\mathsf{T}}H^{\mathrm{mw}}U` evaluated with :math:`s=+1`. With the
tolerance

.. math::

   \tau = \sqrt{\varepsilon_{\mathrm{mach}}}\;
          \max\bigl(1,\ \max_j\lvert\lambda^{\mathrm{int}}_j\rvert\bigr) ,
   \qquad
   \sqrt{\varepsilon_{\mathrm{mach}}}\approx 1.49\times10^{-8} ,

the counts :math:`n_+ = \#\{\lambda^{\mathrm{int}}_j > \tau\}` and
:math:`n_- = \#\{\lambda^{\mathrm{int}}_j < -\tau\}` decide the sign:
:math:`s=-1` when :math:`n_- > n_+`, and :math:`s=+1` when :math:`n_+ > n_-`.
A tie is broken by the larger of
:math:`\sum\lvert\lambda^{\mathrm{int}}\rvert` over the positive and the
negative set. If the structure has no internal subspace at all — a single atom,
whose three coordinates are exhausted by the translations — the heuristic
returns :math:`s=+1`.

The heuristic exists because a bound structure must have positive curvature
along most of its internal coordinates, so the sign that makes the majority of
internal eigenvalues positive is the physical one. This is also the only place
where the internal subspace is genuinely projected out. Being a majority vote,
it is reliable for minima and for transition states — one imaginary mode among
many real ones — and unreliable for structures far from any stationary point.
Set ``hessian_sign`` explicitly whenever the convention of the producing code
is known.

Force constants, reduced masses and IR intensities
--------------------------------------------------

The reduced mass of a mode is the inverse squared length of its unnormalized
Cartesian displacement,

.. math::

   \mu_j
   = \left(\sum_\alpha l_{\alpha j}^{2}\right)^{-1}
   = \left(\sum_\alpha
       \frac{(Q\mathbf{c}_j)_\alpha^{2}}{m_\alpha}\right)^{-1} ,

reported in amu. A purely translational mode reports the total mass divided
by the number of atoms, 6.01 amu for a water molecule; a localized hydrogen
stretch approaches 1 amu.

The force constant combines the frequency with the reduced mass,

.. math::

   k_j = \frac{\omega_j^{2}\,\mu_j}{6.022\times10^{28}} ,

in mdyn·Å⁻¹. The single constant applies both
:math:`1\ \mathrm{amu} = 1.66054\times10^{-27}\ \mathrm{kg}` and
:math:`1\ \mathrm{mdyn\,\AA^{-1}} = 100\ \mathrm{N\,m^{-1}}`. It uses the
rounded value :math:`6.022`, so force constants carry a systematic offset of
about :math:`2\times10^{-5}` relative to the CODATA constant. That is far below
any physical uncertainty in a Hessian, but it is visible when comparing digit
by digit against another program.

With a ``moldescriptor_file`` every atom carries a fixed partial charge
:math:`q_i` in units of the elementary charge, and the infrared intensity is

.. math::

   I_j = \frac{42.2561}{\mu_j}
         \sum_{a\in\{x,y,z\}}
         \left(\frac{1}{0.2081943\,\lVert e_j\rVert}
               \sum_i q_i\,e_{(i,a)\,j}\right)^{2} ,

in km·mol⁻¹. The inner sum :math:`\sum_i q_i e_{(i,a)j}` is the derivative of
the point-charge dipole moment along mode :math:`j` [Person1974]_; dividing by
:math:`0.2081943\ \mathrm{e\,\AA\,D^{-1}}` expresses it in D·Å⁻¹, and
:math:`42.2561` converts
:math:`(\mathrm{D\,\AA^{-1}})^{2}\,\mathrm{amu^{-1}}` to km·mol⁻¹. The norm
:math:`\lVert e_j\rVert` equals one by construction and acts only as a guard.

External modes and the ``modes`` selection
------------------------------------------

A non-linear structure has three translational and three rotational modes; a
linear one has three and two. At an exact stationary point these external modes
carry no curvature and appear at the bottom of the spectrum with
:math:`\tilde{\nu}\approx 0`. In practice they are small but non-zero, and they
may be slightly negative, because the Hessian was computed at a finite
convergence threshold and in finite precision.

``modes`` selects which modes reach ``modes_prefix`` and ``modes_file``,
comparing wavenumbers against the threshold :math:`\theta` set by
``modes_threshold`` (default :math:`10^{-8}` cm⁻¹):

.. list-table:: ``modes`` selection rules
   :class: pq-record-table
   :header-rows: 1
   :widths: 24 76

   * - Value
     - Selected modes
   * - ``all`` (default)
     - Every mode, in order of increasing wavenumber
   * - ``nonzero``
     - :math:`\lvert\tilde{\nu}_j\rvert > \theta`, keeping imaginary modes
   * - ``positive``
     - :math:`\tilde{\nu}_j > \theta`, dropping imaginary modes
   * - Numbers
     - Explicit one-based mode numbers: an integer, a list or a range

Because the default threshold is only :math:`10^{-8}` cm⁻¹, ``positive`` does
not by itself remove residual external modes. For the bundled H₂O fixture it
keeps the three translational modes at 0.02, 0.22 and 0.30 cm⁻¹ alongside the
three internal modes at 1493, 3669 and 3784 cm⁻¹, and drops the three
rotational modes, which come out imaginary between :math:`-50` and
:math:`-37` cm⁻¹. Raise ``modes_threshold`` to a few cm⁻¹, or higher, when the
intent is "internal modes only".

Internal-mode spectrum
----------------------

.. plot:: _plots/vibrations.py
   :alt: Infrared stick spectrum for the water validation fixture
   :caption: IR stick spectrum calculated by PQAnalysis from the bundled H₂O
      structure, Hessian and partial-charge fixtures. Only internal modes above
      100 cm⁻¹ are shown; translational and rotational modes are omitted.

Minimal input
-------------

.. code-block:: text

   structure_file = structure.rst
   hessian_file = hessian.dat
   out_file = wavenumbers.dat
   normal_modes_file = normal_modes.dat
   modes_file = modes.xyz
   modes = positive
   unit = kcal
   hessian_sign = auto

.. code-block:: console

   $ pqanalysis vibrations vibrations.in

``structure_file`` may be a PQ restart or a single-frame XYZ file. ``unit``
describes the Hessian energy unit and accepts ``kcal``, ``hartree`` or ``ev``.
IR intensities are written only when a ``moldescriptor_file`` supplies partial
charges.

Mode output
-----------

``normal_modes_file`` stores the dimensionless Cartesian mode matrix
:math:`e_{\alpha j}`. ``modes_file`` writes one extended-XYZ image per selected
mode, carrying the mode vectors and metadata. ``modes_prefix`` writes one
sinusoidal animation per selected mode, with frame :math:`k` of
``modes_frames`` (default 30) at

.. math::

   \mathbf{x}_i(\varphi_k) = \mathbf{x}_i^{0} + \sin(\varphi_k)\,\mathbf{d}_i ,
   \qquad
   \varphi_k = \frac{2\pi k}{n_{\mathrm{frames}}} .

By default the displacement :math:`\mathbf{d}_i` is the mode scaled so that the
largest atomic displacement equals ``modes_amplitude`` (default 0.25 Å). If
``modes_temperature`` :math:`T` is given, the mode is instead scaled by the
classical thermal factor

.. math::

   \mathbf{d}_i = \mathbf{e}_{ij}\sqrt{\frac{k_\mathrm{B}T}{E_j}} ,
   \qquad
   E_j = \lvert\tilde{\nu}_j\rvert\times
         1.2398419843320026\times10^{-4}\ \mathrm{eV\,cm} ,

with :math:`k_\mathrm{B} = 8.617333262145\times10^{-5}\ \mathrm{eV\,K^{-1}}`
and :math:`E_j` the mode energy obtained from the wavenumber. Modes whose
energy lies at or below the threshold fall back to the fixed amplitude, which
keeps near-zero modes from being drawn with a diverging excursion. Both
scalings are display conventions, not physical vibrational amplitudes.

Validity and interpretation
---------------------------

A trustworthy result has :math:`3N` modes in total, of which six — five for a
linear structure — are external and small compared with the softest internal
mode, and :math:`3N-6` are internal and positive at a minimum, or positive
except for exactly one imaginary mode at a transition state. Residual external
modes of a few tens of cm⁻¹, sometimes negative, indicate an incompletely
optimized geometry or a numerically noisy Hessian rather than physical soft
modes.

The method does not apply, or needs care, in these situations:

* **Away from a stationary point.** The harmonic expansion assumes vanishing
  gradients. PQAnalysis projects out neither the gradient nor the external
  modes, so residual forces leak into the translational and rotational modes
  and mix into the low-wavenumber internal modes.
* **Periodic systems.** The external-mode construction uses a center of mass
  and an inertia tensor, which presumes an isolated structure. A Hessian from a
  periodic calculation can still be diagonalized, but the rotational trial
  vectors and the interpretation of the near-zero modes are not meaningful.
* **Unit and ordering mismatches.** The Hessian coordinate order must match the
  structure atom order exactly, and the energy unit must match ``unit``,
  including the bohr length convention implied by ``hartree``. A mismatch
  yields a plausible-looking spectrum on the wrong scale, not an error.
* **Comparison with experiment.** These are harmonic wavenumbers. They
  systematically exceed observed fundamentals, and PQAnalysis applies no
  empirical scaling factor. Finite-temperature and anharmonic band shapes come
  from the time-correlation route in :doc:`vacf` instead [Thomas2013]_.
* **IR intensities.** Fixed atomic point charges carry no charge flux and no
  electronic polarization, so :math:`I_j` reproduces relative band strengths of
  strongly polar motions at best. Do not report them as quantitative
  absorption coefficients.
* **Degenerate modes.** Within a degenerate set the individual eigenvectors are
  arbitrary up to a rotation inside that subspace. Wavenumbers, force constants
  and the summed intensity are well defined; individual mode vectors are not.

Output and API
--------------

See :ref:`analysis-output-vibrations` for every table and file schema. The main
entry point is :func:`PQAnalysis.analysis.vibrational.api.vibrations`; direct
calculations use
:func:`PQAnalysis.analysis.vibrational.vibrational_analysis.calculate_from_system`.

References
----------

* [Wilson1955]_ is the standard treatment of the mass-weighted Hessian
  eigenvalue problem, normal coordinates and the separation of translation and
  rotation.
* [Person1974]_ defines infrared intensities through dipole moment
  derivatives and polar tensors, the quantity the partial-charge model
  approximates.
* [Thomas2013]_ compares this static normal-mode route with spectra obtained
  from molecular-dynamics time correlation functions, which PQAnalysis
  provides through :doc:`vacf`.

Full entries are listed in :doc:`../references`.
