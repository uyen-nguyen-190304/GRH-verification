"""
base_case.py

Perform the base case verification (k = 1) of GRH for a single quadratic Dirichlet L-function

Functions:
    - iota(eta): Compute the maximum missing-zeros contribution up to height η
    - logarithmic_derivative(...): Finite approximation of the logarithmic derivative L'/L(1 - delta, χ_d)
    - base_case_verify(...): Test flow for specific discriminant d. Record the number of zeros needed to verify the GRH up to height η

Constants:
    - mp.dps — mpmath precision
    - EULER  — Euler-Mascheroni constant
    - PI     — π
    - E      — e (base of natural logarithm)

Usage:
    success, eta, N_used = base_case_verify(
        d, K, eta, eps, lcalc_path, data_dir, log_path, chunk=10
    )
"""

from typing import Tuple, List
from pathlib import Path

import numpy as np 
import mpmath as mp

from .utils.generate_zeros import compute_zeros, compute_intervals, write_zeros, write_intervals
from .utils.von_mangoldt import compute_lambda, write_lambda
from .utils.kronecker_symbol import compute_kronecker, write_kronecker

# =============================== CONSTANTS ===============================

mp.dps = 50             # Working precision
EULER  = mp.euler       # Euler-Mascheroni constant
PI     = mp.pi          # π
E      = mp.e           # Base of natural logarithm 

# =========================== UTILITY FUNCTIONS ===========================

def iota(eta: mp.mpf) -> mp.mpf:
    """
    Purpose:
        Maximum contribution of missing zeros 
    Input:  
        eta (mp.mpf) - Width of window to verify the GRH
    Return:
        mp.mpf value of the maximum contribution
    """
    # Compute the two expressions
    eta2 = eta * eta
    term1 = 1.0 / (1.0 + eta2) + 2.0 / (4.0 + eta2)
    term2 = 12.0 / (9.0 + 4.0 * eta2)

    # Take the minimum and return as an mp.mpf
    return mp.mpf(min(term1, term2))

def logarithmic_derivative(delta: int, K: int, chi_arr: np.ndarray, lambda_arr: np.ndarray, remainder_bound: bool = True) -> mp.mpf:
    """
    Purpose: 
        Compute the partial-series approximation of L'/L at s = 1 - delta (delta < 0)
        Include an optional remainder bound
    Input:
        delta      - Negative integer (mostly used as -1 for the verification)
        K          - Truncation parameters (>= 18)
        chi_arr    - NumPy array of χ(k) for k=0..K (χ(0) unused)
        lambda_arr - NumPy array of Λ(k) for k=0..K (Λ(0) = 0)
        remainder_bound - whether to add the analytic tail bound
    Return:
        mp.mpf approximation to L'/L(s) at s = 1 - delta
    """
    # Input validation
    if not (isinstance(delta, int) and delta < 0):
        raise ValueError("delta must be a negative integer ")
    if not (isinstance(K, int) and K >= 18):
        raise ValueError("K must be an integer greater than or equal to 18")

    # Sum the explicit series
    total = mp.mpf("0")
    for k in range(1, K + 1):
        # Compute the general lambda value lambda_L
        lambda_L = mp.mpf(lambda_arr[k] * chi_arr[k])

        # Add the contribution of the k-th term to the sum
        total -= lambda_L / mp.power(k, 1 - delta)

    # Analytic upper bound of remainder term contribution
    if remainder_bound:
        total += (mp.power(K, delta) / delta) * (2.85 * (2 * delta - 1) / mp.log(K) - 1)

    return total

# =========================== BASE-CASE VERIFICATION ===========================

def base_case_verify(d: int, K: int, eta: float, eps: float, lcalc_path: str | Path, data_dir: str | Path, log_path: str | Path, chunk: int=10) -> Tuple[bool, int]:
    """
    Purpose:
        Verify the Generalized Riemann Hypothesis for the Dirichlet character χ_d with order k = 1 (base case)
    Input:
        d          - Fundamental discriminant
        K          - Truncation for chi/lambda arrays
        eta        - Window of interest to verify the GRH
        eps        - Small interval half-width around each ordinate
        lcalc_path — Path to lcalc executable
        data_dir   - Directory to store the data computed
        log_path   - Path to write any error logs
        chunk      - Number of zeros to process per internal batch
    Return:
        (success: bool, eta_used: float, N_used: int)
    """
    # --------- Pre-compute Kronecker and Λ arrays ---------
    chi_arr    = compute_kronecker(d, K)
    lambda_arr = compute_lambda(K)

    # ------------------------- RHS -------------------------

    # Compute the RHS constant based on the discriminant sign
    if d < 0:
        rhs_const = mp.mpf("0.5") * mp.log((abs(d) * E**2) / (4 * PI * E**EULER))
    else:
        rhs_const = mp.mpf("0.5") * mp.log(abs(d) / (PI * E**EULER))

    # Add the approximation of logarithmic derivative contribution to the RHS
    rhs = rhs_const + logarithmic_derivative(-1, K, chi_arr, lambda_arr, remainder_bound=True)

    # ----------------------- LHS -----------------------

    # Initialize the LHS with the iota(eta) of the missing zeros guard
    lhs = 2 * iota(mp.mpf(eta))

    # Check if we need contribution from any zeros at all to verify the GRH up to height η
    if lhs > rhs:
        return True, eta, 0     # Success to verify up to height eta without any zeros needed

    # Track the zeros and intervals used for verification
    zeros_acc: List[mp.mpf] = []
    intervals_acc: List[tuple[mp.mpf, mp.mpf]] = []

    # Contribution of the zeros
    N_used = 0
    start  = 0
    try:
        # Loop in chunks until we exceed the RHS or exhaust of zeros
        while True:
            # Load a new chunk of zeros and intervals 
            need = start + chunk 
            zeros = compute_zeros(d, need, lcalc_path)
            intervals = compute_intervals(d, need, eps, lcalc_path, zeros=zeros)

            if len(intervals) <= start:
                break   # Exhausted zero list

            # Process exactly the next 'chunk' intervals
            intervals = intervals[start:]
            for index, (gamma_minus, gamma_plus) in enumerate(intervals[start:]):
                gamma_minus = mp.mpf(gamma_minus)
                gamma_plus  = mp.mpf(gamma_plus)

                # Record the zeros and intervals used
                zeros_acc.append(mp.mpf(zeros[start + index]))
                intervals_acc.append((gamma_minus, gamma_plus))

                # Separate the contribution of the zeros by type
                if mp.almosteq(gamma_minus + gamma_plus, 0, rel_eps=0, abs_eps=1e-12):
                    # Type 2: symmetric [-gamma0, gamma0]
                    gamma0 = mp.fabs(gamma_plus)
                    lhs += 6 / (9 + 4 * gamma0 * gamma0)
                else:
                    # Type 1: asymmetric [gamma_minus, gamma_plus]
                    lhs += 12 / (9 + 4 * gamma_plus * gamma_plus)

                # Increment the number of used zeros
                N_used += 1

                # Check if the LHS exceeds the RHS
                if lhs > rhs:
                    raise StopIteration     # Success
            
            start += chunk 
                    
        # Loop exhaust without RH verified
        success = False

    except StopIteration:
        # Exited via success condition inside loop
        success = True

    except Exception as err:
        # Log the error to a file if computation fails for this d
        with open(log_path, "a") as log:
            log.write(f"Error: d = {d}, N = {N_used}, reason = {repr(err)}\n")
        success = False
    
    # Print out the result
    # print(f"Zeros used: {zeros_acc}")
    # print(f"Intervals used: {intervals_acc}")
    # print(f"Iota(eta): {iota(eta)}")
    # print(f"Logarithmic Derivative L'/L: {logarithmic_derivative(-1, K, chi_arr, lambda_arr, remainder_bound=True)}")
    # print(f"Constant RHS term: {rhs_const}")
    # print(f"2 * iota(eta): {2 * iota(mp.mpf(eta))}")
    # print(f"C(Z): {lhs - 2 * iota(mp.mpf(eta))}")
    # print(f"RHS: {rhs}")
    # print(f"LHS: {lhs}")

    # Save zeros, intervals, lambda, and kronecker values used to .txt files
    data_dir = Path(data_dir).expanduser().resolve()
    data_dir.mkdir(parents=True, exist_ok=True)

    write_zeros(d, [float(z) for z in zeros[:N_used]], data_dir)
    write_intervals(d, np.asarray(intervals[:N_used], dtype=float), data_dir)
    write_lambda(K, lambda_arr, data_dir)
    write_kronecker(d, K, chi_arr, data_dir)

    return success, eta, N_used
