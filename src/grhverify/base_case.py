from typing import Tuple, List
from pathlib import Path

import numpy as np 
import mpmath as mp

from utils.generate_zeros import compute_zeros, compute_intervals, write_zeros, write_intervals
from utils.von_mangoldt import compute_lambda, write_lambda
from utils.kronecker_symbol import compute_kronecker, write_kronecker

# =============================== CONSTANTS ===============================

mp.dps = 50             # Working precision
EULER  = mp.euler       # Euler-Mascheroni constant
PI     = mp.pi  
E      = mp.e

# =========================== UTILITY FUNCTIONS ===========================

def iota(eta: mp.mpf) -> mp.mpf:
    """
    Return iota(eta) = min(1 / (1 + eta^2) +  2 / (4 + eta^2), 12 / (9 + 4 * eta^2))
    """
    # Compute the two expressions
    eta2 = eta * eta
    term1 = 1.0 / (1.0 + eta2) + 2.0 / (4.0 + eta2)
    term2 = 12.0 / (9.0 + 4.0 * eta2)

    # Return the minimum value of the two expressions
    return mp.fmin(term1, term2)


def logarithmic_derivative(delta: int, K: int, chi_arr: np.ndarray, lambda_arr: np.ndarray, remainder_bound: bool = True) -> mp.mpf:
    """
    Approximate the logarithmic derivative of the L function at 1 - delta for negative delta
    Implementation of Lemma 6
    """
    # Input validation
    if not (isinstance(delta, int) and delta < 0):
        raise ValueError("delta must be a negative integer ")
    if not (isinstance(K, int) and K >= 18):
        raise ValueError("K must be an integer greater than or equal to 18")

    # Declare the sum variable
    total = mp.mpf("0")

    # Loop over all k from 1 to K
    for k in range(1, K + 1):
        # Compute the general lambda value lambda_L
        lambda_L = mp.mpf(lambda_arr[k]) * mp.mpf(chi_arr[k])

        # Add the contribution of the k-th term to the sum
        total -= lambda_L / mp.power(k, 1- delta)

    # Add the upper bound of remainder term contribution
    if remainder_bound:
        total += (mp.power(K, delta) / delta) * (2.85 * (2 * delta - 1) / mp.log(K) - 1)

    return total

# =========================== BASE CASE VERIFICATION ===========================

def base_case_verify(d: int, K: int, eta: float, eps: float, lcalc_path: str | Path, data_dir: str | Path, log_path: str | Path, chunk: int=10) -> Tuple[bool, int]:
    """
    Verify the inequality for discriminant d
    Return (success, N_used)
    """
    # --------- Pre-compute Kronecker and Λ arrays ---------
    chi_arr    = compute_kronecker(d, K)
    lambda_arr = compute_lambda(K)

    # ------------------------- RHS -------------------------

    # Compute the RHS constant based on the discriminant sign
    if d < 0:
        rhs_const = mp.mpf("0.5") * mp.log(abs(d) * E**2) / (4 * PI * E**EULER)
    else:
        rhs_const = mp.mpf("0.5") * mp.log(abs(d) / PI * E**EULER)

    # Add the approximation of logarithmic derivative contribution to the RHS
    rhs = rhs_const + logarithmic_derivative(-1, K, chi_arr, lambda_arr, remainder_bound=True)

    # ----------------------- LHS -----------------------

    # Initialize the LHS with the iota(eta) of the missing zeros guard
    lhs = 2 * iota(mp.mpf(eta))

    # Track the zeros and intervals used to save to file later
    zeros_acc: List[mp.mpf] = []
    intervals_acc: List[tuple[mp.mpf, mp.mpf]] = []

    # Contribution of the zeros
    N_used = 0
    start  = 0

    try:
        while True:
            # Load a chunk of new intervals 
            need = start + chunk 
            zeros = compute_zeros(d, need, lcalc_path)
            intervals = compute_intervals(d, need, eps, lcalc_path, zeros=zeros)

            if len(intervals) <= start:
                break   # Exhausted zero list

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
        success = True

    except Exception as err:
        # Log the error to a file if computation fails for this d
        with open(log_path, "a") as log:
            log.write(f"Error: d = {d}, N = {N_used}, reason = {repr(err)}\n")
        success = False
    
    # Save zeros, intervals, lambda, and kronecker values used to .txt files
    data_dir = Path(data_dir).expanduser().resolve()
    data_dir.mkdir(parents=True, exist_ok=True)

    write_zeros(d, [float(z) for z in zeros_acc], data_dir)
    write_intervals(d, np.asarray(intervals_acc, dtype=float), data_dir)
    write_lambda(K, lambda_arr, data_dir)
    write_kronecker(d, K, chi_arr, data_dir)

    return success, N_used
