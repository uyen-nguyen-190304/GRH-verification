import os
import math
import argparse
import numpy as np
from pathlib import Path
from sage.all import prime_range

def compute_lambda(K: int) -> np.ndarray:
    """
    Compute the von Mangoldt function values Λ(k) for 1 <= k <= K

    Input:  K (int): Upper bound for k in Λ(k)

    Output: Array of shape (K +1, ) where lambda_arr[k] = Λ(k) and 
            lambda_arr[0] is unused (set to 0)
    """
    lamba_arr = np.zeros(K + 1, dtype=float)
    for p in prime_range(2, K + 1):
        log_p = math.log(p)
        q = p
        while q <= K:
            lamba_arr[q] = log_p
            q *= p
    return lamba_arr

def write_lambda(lambda_arr: np.ndarray, K: int, data_dir: Path) -> None:
    """
    Save the von Mangoldt array in both .npy and .txt formats

    Input:  lambda_arr (np.ndarray): Array containing von Mangoldt values
            K (int): Upper bound for k in Λ(k)
            data_dir (Path): Path to data directory
    
    Output: None
    """
    data_dir.mkdir(parents=True, exist_ok=True)

    # Determine cache folder and ensure it exists
    cache_dir = data_dir / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Save binary .npy file
    npy_path = cache_dir / f"von_mangoldt_{K}.npy"
    np.save(npy_path, lambda_arr)

    # Save readable text file
    txt_path = data_dir / "von_mangoldt.txt"
    with open(txt_path, "w") as f:
        for k in range(1, K + 1):
            f.write(f"{k} {lambda_arr[k]}\n")

def main():
    """
    Command-line interfact for computing and caching von Mangoldt values
    """
    parser = argparse.ArgumentParser(description="Compute von Mangoldt function values")
    parser.add_argument("-K", "--upper-limit", type=int, required=True, help="Upper limit for k in von Mangoldt function")
    parser.add_argument("-data", "--data-dir", type=str, default="data", help="Output directory to save Λ(k) values (default: data)")
    args = parser.parse_args()

    # Retrieve command line parameters
    K = args.upper_limit
    data_dir = Path(args.data_dir)

    # Validate input
    if K <= 0:
        raise ValueError("Upper limit K must be a positive integer")
    
    # Skip work if cache already present
    npy_path = data_dir / "cache" / f"von_mangoldt_{K}.npy"
    if npy_path.exists(): 
        return
    
    # Else, compute Λ(k) array and save to file
    lambda_arr = compute_lambda(K)
    write_lambda(lambda_arr, K, data_dir)

if __name__ == "__main__":
    main()
