from runtest import get_filter
from typing import List


def get_general_filter():

    f = [
        get_filter(string="Memory available for calculation:", abs_tolerance=1.0e-6),
    ]

    return f


def get_hf_filter(
    tolerance: float,
    **kwargs,
):
    convergence = kwargs.get("convergence", False)
    restart = kwargs.get("restart", False)
    idempotent = kwargs.get("idempotent", True)
    HOMO_LUMO = kwargs.get("HOMO_LUMO", True)

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


def get_mlhf_filter(tolerance: float, **kwargs):
    f = get_hf_filter(tolerance, **kwargs)

    g = [
        get_filter(string="Active energy:", abs_tolerance=tolerance),
        get_filter(string="Active-inactive energy:", abs_tolerance=tolerance),
        get_filter(string="Inactive energy:", abs_tolerance=tolerance),
    ]
    f.extend(g)
    return f


def get_tdhf_filter(tolerance: float, n_states: int, **kwargs):
    f = get_hf_filter(tolerance, **kwargs)
    f.extend(get_es_filter(n_states, tolerance))
    return f


def get_gs_filter(tolerance: float, **kwargs):
    f = get_hf_filter(tolerance, **kwargs)

    g = [
        get_filter(
            string="Final ground state energy (a.u.):",
            abs_tolerance=tolerance,
        ),
        get_filter(
            from_string="Largest single amplitudes:",
            num_lines=5,
            abs_tolerance=tolerance,
            mask=[3],
            ignore_sign=True,
        ),
    ]

    if kwargs.get("newton", False):
        g.append(
            get_filter(
                from_string="Iteration    Energy (a.u.)        |omega|       Delta E (a.u.)",
                num_lines=15,
                abs_tolerance=tolerance,
                mask=[2, 3],
                ignore_sign=True,
            ),
        )
    f.extend(g)
    return f


def get_cc_es_filter(n_states: int, tolerance: float, **kwargs):
    f = get_gs_filter(tolerance, **kwargs)
    f.extend(get_es_filter(n_states, tolerance))
    return f


def get_es_filter(n_states: int, tolerance: float):
    f = [
        get_filter(
            from_string="State                (Hartree)             (eV)",
            num_lines=2 + n_states,
            abs_tolerance=tolerance,
            mask=[2],
        )
    ]
    return f


def get_eom_filter(n_states: int, tolerance: float, **kwargs):
    f = get_cc_es_filter(n_states, tolerance, **kwargs)

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


def get_td_filter(tolerance: float):
    """
    Returns filters for an timedependent calculation.
    """
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


def get_mean_value_filter(tolerance: float, n_components: int, **kwargs):
    """
    Returns filters for a mean value calculation.
    """
    f = [
        get_filter(
            from_string="Comp.         Electronic           Nuclear             Total",
            num_lines=n_components + 2,
            abs_tolerance=tolerance,
        ),
    ]

    if kwargs.get("norm", False):
        f.append(
            get_filter(
                string="Norm of the total",
                abs_tolerance=tolerance,
            )
        )
    return f


def get_file_filter(tolerance: float):
    """
    Returns filters to check an entire file.
    """
    return [
        get_filter(
            abs_tolerance=tolerance,
            ignore_sign=True,
        )
    ]


def get_spin_filter(tolerance: float):
    f = [
        get_filter(string="Sz:", abs_tolerance=tolerance),
        get_filter(string="Sz(Sz + 1):", abs_tolerance=tolerance),
        get_filter(string="S^2:", abs_tolerance=tolerance),
        get_filter(string="Spin contamination:", abs_tolerance=tolerance),
    ]
    return f


def get_polarizability_filter(tolerance: float, components: List[int]):
    """
    Returns filter for polarizabilities, where components is a list
    of integer 11, 21, 31, 22, 32, 33 specifying the components of
    the dipole polarizability analogous to the eT input.
    """
    f = []
    component_key = {"1": "x", "2": "y", "3": "z"}
    for number in components:
        lower, higher = sorted(str(number))
        component1 = component_key[higher]
        component2 = component_key[lower]
        f.append(
            get_filter(
                string=f"<< mu_{component1}, mu_{component2} >>",
                abs_tolerance=tolerance,
            ),
        )
    return f


def get_fci_filter(
    tolerance: float,
    n_states: int,
    **kwargs,
):
    f = get_hf_filter(tolerance, **kwargs)

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
    f = [
        # first line in geometries
        get_filter(from_string="Z-matrix (", abs_tolerance=1.0e-6, num_lines=4)
    ]

    return f


def get_geoopt_filter(tolerance: float):
    """
    Returns filter for geometry optimization
    """
    f = [
        get_filter(
            from_string="- Summary of geometry optimization:",
            num_lines=15,
            abs_tolerance=tolerance,
        ),
    ]
    return f


def get_lanczos_filter(tolerance: float, **kwargs):
    f = get_gs_filter(tolerance, **kwargs)

    f.append(
        get_filter(
            from_string="State.      energy [a.u]         energy [eV]         Osc. strength",
            num_lines=12,
            abs_tolerance=tolerance,
            mask=[2, 4],
        )
    )
    return f
