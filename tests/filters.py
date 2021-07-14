def get_general_filter():
    """
    Returns filter for general settings
    """
    from runtest import get_filter

    f = [
        get_filter(string="Memory available for calculation:", abs_tolerance=1.0e-6),
    ]

    return f


def get_hf_filter(
    tolerance, convergence=False, restart=False, idempotent=True, z_matrix=False
):
    """
    Returns filters for a HF calculation.
    """

    from runtest import get_filter

    f = get_general_filter()

    g = [
        # Energy quantities
        get_filter(string="Nuclear repulsion energy:", abs_tolerance=tolerance),
        get_filter(string="HOMO-LUMO gap", abs_tolerance=tolerance),
        get_filter(string="Total energy:", abs_tolerance=tolerance),
    ]

    if idempotent:
        g.append(
            get_filter(
                from_string="Iteration       Energy (a.u.)      Max(grad.)    Delta E (a.u.)",
                num_lines=3,
                abs_tolerance=tolerance,
                mask=[2],
            )
        )

    if convergence:
        g.append(
            get_filter(string="Convergence criterion met in", abs_tolerance=1.0e-10)
        )

    if not restart:
        g.append(
            # non-idempotent SAD
            get_filter(string="Energy of initial guess:", abs_tolerance=1.0e-6)
        )

    if z_matrix:
        g.append(
            # first line in geometries
            get_filter(from_string="Z-matrix (", abs_tolerance=1.0e-6, num_lines=4)
        )

    f.extend(g)

    return f


def get_mlhf_filter(tolerance, convergence=True, restart=False, idempotent=True):
    """
    Returns filters for a MLHF calculation.
    """

    from runtest import get_filter

    f = get_hf_filter(tolerance, convergence, restart, idempotent)

    g = [
        # active energy
        get_filter(string="Active energy:", abs_tolerance=tolerance),
        # active-inactive energy
        get_filter(string="Active-inactive energy:", abs_tolerance=tolerance),
        # inactive energy
        get_filter(string="Inactive energy:", abs_tolerance=tolerance),
    ]

    f.extend(g)

    return f


def get_tdhf_filter(
    n_states,
    tolerance,
    convergence=False,
    restart=False,
    idempotent=True,
    z_matrix=False,
):
    """
    Returns filters for a TDHF calculation.
    """

    from runtest import get_filter

    f = get_hf_filter(tolerance, convergence, restart, idempotent, z_matrix)
    g = [
        get_filter(
            from_string="State                (Hartree)             (eV)",
            num_lines=2 + n_states,
            abs_tolerance=tolerance,
            mask=[2],
        )
    ]

    f.extend(g)

    return f


def get_gs_filter(tolerance, convergence=True, restart=False, Newton=False):
    """
    Returns filters for a GS calculation.
    """

    from runtest import get_filter

    f = get_hf_filter(tolerance, convergence, restart)

    g = [
        get_filter(string="Final ground state energy (a.u.):", abs_tolerance=tolerance),
        get_filter(
            from_string="Largest single amplitudes:",
            num_lines=5,
            abs_tolerance=tolerance,
            mask=[3],
            ignore_sign=True,
        ),
    ]

    if Newton:
        h = [
            get_filter(
                from_string="Iteration    Energy (a.u.)        |omega|       Delta E (a.u.)",
                num_lines=13,
                abs_tolerance=tolerance,
                mask=[3],
                ignore_sign=True,
            ),
        ]
        g.extend(h)

    f.extend(g)

    return f


def get_es_filter(n_states, tolerance, convergence=True, restart=False, Newton=False):
    """
    Returns filters for an ES calculation.
    """

    from runtest import get_filter

    f = get_gs_filter(tolerance, convergence, restart, Newton)

    g = [
        get_filter(
            from_string="State                (Hartree)             (eV)",
            num_lines=2 + n_states,
            abs_tolerance=tolerance,
            mask=[2],
        )
    ]

    f.extend(g)

    return f


def get_eom_filter(n_states, tolerance, convergence=True, restart=False, Newton=False):
    """
    Returns filters for an EOM calculation.
    """

    from runtest import get_filter

    f = get_es_filter(n_states, tolerance, convergence, restart, Newton)

    g = [
        get_filter(
            from_string="Comp. q     < n |q| m >       < m |q| n >        < n |q| m > < m |q| n >",
            num_lines=5,
            abs_tolerance=tolerance,
            mask=[1, 2],
            ignore_sign=True,
        ),
        get_filter(
            from_string="Comp. q     < n |q| m >       < m |q| n >        < n |q| m > < m |q| n >",
            num_lines=5,
            abs_tolerance=tolerance,
            mask=[3],
        ),
        get_filter(string="Oscillator strength:", abs_tolerance=tolerance),
    ]

    f.extend(g)

    return f


def get_td_filter(tolerance):
    """
    Returns filters for an timedependent calculation.
    """

    from runtest import get_filter

    f = [
        get_filter(
            from_string="Energy after propagation [au]:",
            num_lines=5,
            abs_tolerance=tolerance,
        ),
        get_filter(
            from_string="Dipole moment after propagation [au]:",
            num_lines=7,
            abs_tolerance=tolerance,
        ),
    ]

    return f


def get_dipole_filter(tolerance):
    """
    Returns filters for a calculation of dipoles.
    """

    from runtest import get_filter

    f = [
        get_filter(string=" x: ", abs_tolerance=tolerance),
        get_filter(string=" y: ", abs_tolerance=tolerance),
        get_filter(string=" z: ", abs_tolerance=tolerance),
    ]

    return f


def get_quadrupole_filter(tolerance):
    """
    Returns filters for a calculation of quadrupoles.
    """

    from runtest import get_filter

    f = [
        get_filter(string=" xx: ", abs_tolerance=tolerance),
        get_filter(string=" xy: ", abs_tolerance=tolerance),
        get_filter(string=" xz: ", abs_tolerance=tolerance),
        get_filter(string=" yy: ", abs_tolerance=tolerance),
        get_filter(string=" yz: ", abs_tolerance=tolerance),
        get_filter(string=" zz: ", abs_tolerance=tolerance),
    ]

    return f


def get_file_filter(tolerance):
    """
    Returns filters to check an entire file.
    """

    from runtest import get_filter

    f = [get_filter(abs_tolerance=tolerance, ignore_sign=True)]

    return f


def get_spin_filter(tolerance):
    """
    Returns filters for spin summary
    """

    from runtest import get_filter

    f = [
        get_filter(string="Sz:", abs_tolerance=tolerance),
        get_filter(string="Sz(Sz + 1):", abs_tolerance=tolerance),
        get_filter(string="S^2:", abs_tolerance=tolerance),
        get_filter(string="Spin contamination:", abs_tolerance=tolerance),
    ]

    return f


def get_polarizability_filter(tolerance):
    """
    Returns filter for polarizability summary
    """
    from runtest import get_filter

    f = [
        get_filter(string="<< mu_x, mu_x >>(0.20E-01):", abs_tolerance=tolerance),
        get_filter(string="<< mu_z, mu_x >>(0.20E-01):", abs_tolerance=tolerance),
        get_filter(string="<< mu_x, mu_x >>(0.40E-01):", abs_tolerance=tolerance),
        get_filter(string="<< mu_z, mu_x >>(0.40E-01):", abs_tolerance=tolerance),
        get_filter(string="<< mu_x, mu_x >>(0.60E-01):", abs_tolerance=tolerance),
        get_filter(string="<< mu_z, mu_x >>(0.60E-01):", abs_tolerance=tolerance),
    ]

    return f
