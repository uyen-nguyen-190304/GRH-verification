"""
grhverify: Python bindings for verifying the Generalized Riemann Hypothesis
using C++ numerical routines

Exposes:
    - EULER_CONSTANT
    - iota(eta)
    - C(Z)
    - log_derivative(chi, lam, K)
    - zero_contribution(gamma_minus, gamma_plus)
    - rh_verify(d, K, eta, intervals, chi, lam)
"""

# Import everything from the C++ extension module
from .grhverify import *
