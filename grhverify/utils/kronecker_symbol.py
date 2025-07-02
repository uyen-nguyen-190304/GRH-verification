"""
kronecker_symbol.py

Compute and persist the Kronecker symbol χ_n(d) for quadratic Dirichlet characters

Functions:
    - compute_kronecker(d, K): Build the character values χ_d(k) = (d|k) for k=1..K via Sage's kronecker_symbol function
    - write_kronecker(d, K, chi_arr, data_dir): Write the array of χ_d(k) values to a text file
"""

import numpy as np
from pathlib import Path
from sage.all import kronecker_symbol

def compute_kronecker(d: int, K: int) -> np.ndarray:
    """
    Purpose:
        Compute the Kronecker symbol for integers k from [1..K] with respect to d
    Input:  
        d (int): Discriminant of Dirichlet character
        K (int): Upper bound of k
    Return: 
        Array of shape (K + 1, ), where chi_arr[k] is the kronecker symbol of k and chi_arr[0] is unused (set to 0)
    """
    # Input validation
    if not (isinstance(K, int) and K >= 1):
        raise ValueError(f"The upper bound K must be a positive integer")
    
    # Use SageMath to compute the Kronecker symbol of k = 1..K
    chi_arr = np.zeros(K + 1, dtype=np.int8)
    for k in range(1, K + 1):
        chi_arr[k] = kronecker_symbol(d, k)
    return chi_arr


def write_kronecker(d: int, K: int, chi_arr: np.ndarray, data_dir = str | Path) -> None:
    """
    Purpose:
        Save the Kronecker symbol array as a text file in the specified directory
    Input:  
        d (int): Discriminant of Dirichlet character
        K (int): Upper bound for k
        chi_arr (np.ndarray): Array containing Kronecker symbol
        data_dir (str): Path to data directory
    Return:
        None
    """  
    # Ensure the storing directory exists
    target = Path(data_dir).expanduser() / ("positive_d" if d > 0 else "negative_d") / f"d_{d}"
    target.mkdir(parents=True, exist_ok=True)

    # Save the kronecker values to the kronecker.txt file
    txt_path = target / "kronecker.txt"
    with open(txt_path, "w") as f:
        for k in range(1, K + 1):
            f.write(f"{k} {chi_arr[k]}\n")  
            