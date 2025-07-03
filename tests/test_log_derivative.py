import pytest
import subprocess
import numpy as np
import mpmath as mp

from grhverify.base_case import logarithmic_derivative
from grhverify.utils.discriminant import is_fundamental_discriminant
from grhverify.utils.von_mangoldt import compute_lambda
from grhverify.utils.kronecker_symbol import compute_kronecker

# ======================== WORKING CONSTANTS ========================

mp.dps    = 20         # Working precision
K         = 100000     # Truncating limit for L'/L
delta     = -1         # s = 2
h = mp.mpf("1e-5")      
tolerance = mp.mpf("1e-8")

# ==================== HELPER: CALL MATHEMATICA =====================

# Numerical differentiation to approximate L'/L
def mathematica_log_derivative(d: int, delta: int, K: int, accuracy: int=20) -> float:
    scheme = scheme.lower()
    code = f"""
    accuracy = {accuracy};
    d = {d};
    q = Abs[d];
    K = {K};
    h = {h};
    s = {1 - delta};

    Lval = N[DirichletL[q, (q + 1)/2, s], accuracy];
    Lval_plus = N[DirichletL[q, (q + 1)/2, s + h], accuracy];
    Lval_minus= N[DirichletL[q, (q + 1)/2, s - h], accuracy];

    centralFD = (Lval_plus - Lval_minus) / (2 h * Lval);
    forwardFD = (Lval_plus - Lval) / (h * Lval);

    centralFD := (logLplus - logLminus)/(2 h);
    forwardFD := (logLplus - logL)/(h);

    result = If["{scheme}" === "central", centralFD, forwardFD];
    N[result, accuracy]
    """
    result = subprocess.run(
        ["wolframscript", "-code", code],
        capture_output=True, text=True, check=True
    )
    return float(result.stdout.strip())

# ======================= CHOOSE SAMPLE d VALUES =======================

# Can change high and low here; current test size = 50
# Also, current not set seed -> dynamic d choices
d_list = [d for d in np.random.randint(-100000, 100000, size=50) if is_fundamental_discriminant(d)]

# ======================= TEST =======================

@pytest.mark.parametrize("d", d_list)
@pytest.mark.parametrize("scheme", ["central", "forward"])
def test_python_vs_mathematica_fd(d, scheme):
    chi = compute_kronecker(d, K)
    lam = compute_lambda(K)

    analytic_val = logarithmic_derivative(delta, K, chi, lam, remainder_bound=False)
    numerical_diff_val = mathematica_log_derivative(d, float(delta), K, h)

    diff = abs(analytic_val - numerical_diff_val)
    assert diff < tolerance, (
        f"Mismatch for d = {d}\n"
        f"Python analytic: {analytic_val}\n"
        f"Numerical Differentiation {scheme}: {numerical_diff_val}\n"
        f"Diff: {diff}"
    )
