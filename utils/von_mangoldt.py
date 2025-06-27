import math
import numpy as np
from pathlib import Path
from sage.all import prime_range

def compute_lambda(K: int) -> np.ndarray:
    """
    Compute the von Mangoldt function values Λ(k) for 1 <= k <= K

    Input:  K (int): Upper bound for k in Λ(k)

    Output: Array of shape (K +1,) where lambda_arr[k] = Λ(k) and 
            lambda_arr[0] is unused (set to 0)
    """
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
    Save the von Mangoldt array as a text file in the specified directory

    Input:  K (int): Upper bound for k in Λ(k)
            lambda_arr (np.ndarray): Array containing von Mangoldt values
            data_dir (str | Path): Path to data directory
    
    Output: None
    """
    # Ensure the storing directory exists
    base = Path(data_dir).expanduser()
    base.mkdir(parents=True, exist_ok=True)

    # Save the lambda values to the von_mangoldt.txt file
    txt_path = base / "von_mangoldt.txt"
    with open(txt_path, "w") as f:
        for k in range(1, K + 1):
            f.write(f"{k} {lambda_arr[k]}\n")
