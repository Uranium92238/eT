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
    tolerance, convergence=False, restart=False, idempotent=True, HOMO_LUMO=True
):
    """
    Returns filters for a HF calculation.
    """

    from runtest import get_filter

    f = get_general_filter()

    g = [
        # Energy quantities
        get_filter(string="Nuclear repulsion energy:", abs_tolerance=tolerance),
        get_filter(string="Total energy:", abs_tolerance=tolerance),
    ]

    if HOMO_LUMO:
        g.append(
            get_filter(string="HOMO-LUMO gap", abs_tolerance=tolerance),
        )

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
        g.append(get_filter(string="Convergence criterion met in", abs_tolerance=1))

    if not restart:
        g.append(
            # non-idempotent SAD
            get_filter(string="Energy of initial guess:", abs_tolerance=1.0e-6)
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
):
    """
    Returns filters for a TDHF calculation.
    """

    from runtest import get_filter

    f = get_hf_filter(tolerance, convergence, restart, idempotent)
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
                num_lines=15,
                abs_tolerance=tolerance,
                mask=[2, 3],
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


def get_mean_value_filter(tolerance, n_components, norm=False):
    """
    Returns filters for a mean value calculation.
    """

    from runtest import get_filter

    f = [
        get_filter(
            from_string="Comp.         Electronic           Nuclear             Total",
            num_lines=n_components + 2,
            abs_tolerance=tolerance,
        ),
    ]

    if norm:
        f.append(
            get_filter(
                string="Norm of the total",
                abs_tolerance=tolerance,
            )
        )

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


def get_polarizability_filter(tolerance, components):
    """
    Returns filter for polarizability summary
    """
    from runtest import get_filter

    f = []

    for i in components:
        if i == 11:
            f.append(get_filter(string="<< mu_x, mu_x >>", abs_tolerance=tolerance))
        elif i == 12 or i == 21:
            f.append(get_filter(string="<< mu_y, mu_x >>", abs_tolerance=tolerance))
        elif i == 13 or i == 31:
            f.append(get_filter(string="<< mu_z, mu_x >>", abs_tolerance=tolerance))
        elif i == 22:
            f.append(get_filter(string="<< mu_y, mu_y >>", abs_tolerance=tolerance))
        elif i == 23 or i == 32:
            f.append(get_filter(string="<< mu_z, mu_y >>", abs_tolerance=tolerance))
        elif i == 33:
            f.append(get_filter(string="<< mu_z, mu_z >>", abs_tolerance=tolerance))

    return f


def get_fci_filter(
    tolerance, n_states, convergence=True, restart=False, HOMO_LUMO=True
):
    """
    Returns filters for a FCI calculation.
    """
    from runtest import get_filter

    f = get_hf_filter(tolerance, convergence, restart, HOMO_LUMO=HOMO_LUMO)

    g = [
        get_filter(
            from_string="State                (Hartree)             (eV)",
            num_lines=2 + n_states,
            abs_tolerance=tolerance,
        ),
        get_filter(string="Spin Multiplicity:", abs_tolerance=tolerance),
        get_filter(
            from_string="Largest CI amplitudes",
            to_re=r"^$",
            abs_tolerance=tolerance,
            ignore_sign=True,
            ignore_order=True,
        ),
    ]

    f.extend(g)
    return f


def get_z_matrix_filter():
    """
    Returns filter for the geometry print in z-matrix format
    """
    from runtest import get_filter

    f = [
        # first line in geometries
        get_filter(from_string="Z-matrix (", abs_tolerance=1.0e-6, num_lines=4)
    ]

    return f


def get_geoopt_filter(tolerance):
    """
    Returns filter for geometry optimization
    """
    from runtest import get_filter

    f = [
        get_filter(
            from_string="- Summary of geometry optimization:",
            num_lines=15,
            abs_tolerance=tolerance,
        ),
    ]

    return f


def get_lanczos_filter(tolerance, convergence=True, restart=False, Newton=False):
    """
    Returns filters for a Lanczos calculation
    """

    from runtest import get_filter

    f = get_gs_filter(tolerance, convergence, restart, Newton)

    f.append(
        get_filter(
            from_string="State.      energy [a.u]         energy [eV]         Osc. strength",
            num_lines=12,
            abs_tolerance=tolerance,
            mask=[2, 4],
        )
    )
    return f
