"""
grhverify
~~~~~~~~~
Algorithms for verifying the Generalized Riemann Hypothesis (GRH)
via the base-case and higher-power explicit inequalities.

Public re-exports:
    /* utils */
    - compute_zeros
    - compute_intervals
    - compute_lambda
    - compute_kronecker

    /* base_case.py */
    - iota
    - logarithmic_derivative
    - base_case_verify

    /* higher_power.py */

"""

from importlib.metadata import version, PackageNotFoundError

# ----------------------------------------------------------------------
# Package metadata
# ----------------------------------------------------------------------
try:
    __version__: str = version("grhverify")
except PackageNotFoundError:          # running from source, not pip-installed
    __version__ = "0.0.0"

__author__: str = "Uyen Nguyen"
__all__: list[str]                    # declared later (after imports)

# ----------------------------------------------------------------------
# Public API re-exports
# ----------------------------------------------------------------------
from .utils.generate_zeros import compute_zeros, compute_intervals
from .utils.von_mangoldt import compute_lambda
from .utils.kronecker_symbol import compute_kronecker

from .base_case import iota, logarithmic_derivative, base_case_verify
# from .higher_power import logarithmic_derivative, iota

__all__ = [
    "compute_zeros",
    "compute_intervals",
    "compute_lambda",
    "compute_kronecker",

    "iota",
    "logarithmic_derivative",
    "base_case_verify",
    
    "__version__",
]
