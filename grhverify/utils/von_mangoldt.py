"""
von_mangoldt.py

Compute and persist the von Mangoldt function values Λ(n) for n up to a given bound

Functions:
    - compute_lambda(K): Build an array Λ[0..K] of von Mangoldt values
    - write_lambda(K, lambda_arr, data_dir): Save the Λ-array to a text file under a specified directory

"""

import math
import numpy as np
from pathlib import Path
from sage.all import prime_range    # Generate primes in [a, b)

def compute_lambda(K: int) -> np.ndarray:
    """
    Purpose:
        Compute the von Mangoldt function values Λ(k) for 1 <= k <= K
        Λ(n) = log(p) exactly when n is a prime power p^k, and zero otherwise
    Input:  
        K (int): Upper bound for k in Λ(k)
    Return: 
        Array of shape (K +1,) where lambda_arr[k] = Λ(k) and lambda_arr[0] is unused (set to 0)
    """
    # Input validation
    if not (isinstance(K, int) and K >= 1):
        raise ValueError(f"The upper bound K must be a positive integer")
    
    # Compute the von Mangoldt lambda value for k = 1..K
    lambda_arr = np.zeros(K + 1, dtype=float)
    for p in prime_range(2, K + 1):
        log_p = math.log(p)
        q = p
        while q <= K:
            lambda_arr[q] = log_p
            q *= p
    return lambda_arr


def write_lambda(K: int, lambda_arr: np.ndarray, data_dir: str | Path) -> None:
    """
    Purpose:
        Save the von Mangoldt array as a text file in the specified directory
    Input:  
        K (int): Upper bound for k in Λ(k)
        lambda_arr (np.ndarray): Array containing von Mangoldt values
        data_dir (str | Path): Path to data directory
    Return:
        None
    """
    # Ensure the storing directory exists
    base = Path(data_dir).expanduser()
    base.mkdir(parents=True, exist_ok=True)

    # Save the lambda values to the von_mangoldt.txt file
    txt_path = base / "von_mangoldt.txt"
    with open(txt_path, "w") as f:
        for k in range(1, K + 1):
            f.write(f"{k} {lambda_arr[k]}\n")
