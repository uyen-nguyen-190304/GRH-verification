import pytest
import numpy as np
import mpmath as mp

from grhverify.base_case import logarithmic_derivative
from grhverify.utils.von_mangoldt import compute_lambda
from grhverify.utils.kronecker_symbol import compute_kronecker
from driver import is_fundamental_discriminant

# =============================== CONSTANTS ===============================

mp.dps = 20             # Working precision 
delta  = -1             # Evaluate the logarithmic derivative at delta = -1 -> s = 2
K      = 100000         # Number of terms to compute the finite sum of the logarithmic derivative
d      = 10009          # Fundamental discriminants
tol    = mp.mpf("1e-5") # Current tolerance level              

# List of fundamental discriminants to test
d_list = [d for d in range(-100, 100) if is_fundamental_discriminant(d)]

# ================================= HELPER =================================

def _dirichlet_L_trunc(s: mp.mpf, K: int, chi: np.ndarray) -> mp.mpf:
    """Truncated Dirichlet series"""
    return mp.fsum(chi[k] / mp.power(k, s) for k in range(1, K + 1))


def _fd_log_derivative(delta: int, K: int, chi: np.ndarray, scheme: str = "central", h: mp.mpf = mp.mpf("1e-5")) -> mp.mpf:
    """
    Finite-difference estimate of L'(s)/L(s)

    scheme = "forward"  -> (log L(s + h) - log L(s)) / h
    scheme = "central"  -> (log L(s + h) - log L(s - h)) / 2h
    """
    s = mp.mpf(1 - delta)          # δ < 0 ⇒ s > 1
    if scheme == "forward":
        return (mp.log(_dirichlet_L_trunc(s + h, K, chi)) - mp.log(_dirichlet_L_trunc(s, K, chi))) / h
    if scheme == "central":
        return (mp.log(_dirichlet_L_trunc(s + h, K, chi)) - mp.log(_dirichlet_L_trunc(s - h, K, chi))) / (2 * h)
    raise ValueError("scheme must be 'forward' or 'central'")

# ================================= TESTS =================================

@pytest.mark.parametrize("d", d_list)
def test_log_derivative(d):
    """Analytic value (minus tail) approx forward & central finite differences"""
    chi = compute_kronecker(d, K)
    lam = compute_lambda(K)

    # Analytic term -> Lemma 6
    analytic = logarithmic_derivative(delta, K, chi, lam, remainder_bound=False)

    # Finite difference estimate
    fd_forward = _fd_log_derivative(delta, K, chi, scheme="forward")
    fd_central = _fd_log_derivative(delta, K, chi, scheme="central")

    # Test for value match-up
    assert abs(analytic - fd_central) < fd_central, (
        f"d={d}: central FD diff = {abs(analytic - fd_central)} exceeds the tolerance {fd_central}"
    )
    assert abs(analytic - fd_forward) < fd_forward, (
        f"d={d}: forward FD diff = {abs(analytic - fd_forward)} exceeds the tolerance {fd_forward}"
    )
