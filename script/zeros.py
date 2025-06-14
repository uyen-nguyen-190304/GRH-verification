import os
import re
import sys
import subprocess
import numpy as np
from typing import List

def compute_zeros(d: int, N: int, lcalc_path: str):
    """
    Compute the first N postive ordinates of the Dirichlet L-function associated
    with discriminant d using lcalc

    Input:  d (int): Fundamental discriminant of the Dirichlet character χ_d
            N (int): Number of non-trivial zeros to compute
            lcalc_path (str): Path to the lcalc executable

    Output: List of the first N non-trivial zeros ordinates (float)
    """
    # Validate inputs
    if not isinstance(d, int):
        sys.exit("Error: discriminant d must of type integer")
    if not isinstance(N, int) or N <= 0:
        sys.exit("Error: Number of zeros N must be a positive integer")

    # Get the path to lcalc executable
    lcalc_path = os.path.expanduser(lcalc_path)
    if not os.path.isfile(lcalc_path):
        sys.exit(f"Error: lcalc executable not found at {lcalc_path}")

    # Construct lcalc command line
    # TODO: This still doesn't work with negative discriminant
    cmd = [
        lcalc_path,
        "-z", str(N),           # Number of zeros to compute
        "--twist-quadratic",    # Twist the L-function by a quadratic character
        "--start", str(d),      # Discriminant d
        "--finish", str(d)
    ]

    # Run the command and capture output
    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
        output = result.stdout
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error running lcalc: \nSTDOUT: {e.stdout}\nSTDERR: {e.stderr}")
        sys.exit(1)

    # Extract zeros from the output 
    zeros = []
    for line in output.splitlines():
        line = line.strip()
        if re.match(r'^\d+\s+[\d\.]+', line):
            parts = line.split()
            if len(parts) >= 2:
                gamma = float(parts[1])     # Extract the imaginary part of the zero
                zeros.append(gamma)

    return zeros

def fetch_zeros(d: int, start: int, count: int, lcalc_path: str) -> List[float]:
    """
    Fetch new zeros 

    Input:  d (int): Fundamental discriminant of interest
            start (int): Number of zeros already cached (zero-based)
            count (int): Number of additional zeros needed
            lcalc_path (str): Path to the lcalc executable

    Output: 
    """
    # If count is negative or we don't need new zeros, return empty list
    if count <= 0:
        return []
    
    # Fetch ful list of zeros up to start + count, then keep only the tail
    full = compute_zeros(d, start + count, lcalc_path)
    return full[start:]         

def cached_zeros(d: int, N: int, lcalc_path: str, cache_root: str="data") -> np.ndarray:
    """
    Return the first N zeros of χ_d from cache if available, otherwise compute and cache them

    Input:  d (int): Discriminant of the Dirichlet character
            N (int): Number of zeros (intervals) to return
            lcalc_path (str): Path to the lcalc executable
            cache_root (str): Root folder for cached data (default: data)

    Output: Array of shape (N, ) containing the first N zeros (np.ndarray)
    """
    # Construct the path to zeros.npy file
    parent_dir = "positive_d" if d > 0 else "negative_d"
    dir_path   = os.path.join(cache_root, parent_dir, f"d_{d}")
    os.makedirs(dir_path, exist_ok=True)
    zfile = os.path.join(dir_path, "zeros.npy")

    # Load cached zeros if file exists, otherwise start with empty
    zeros = np.load(zfile) if os.path.exists(zfile) else np.empty(0, dtype=float)
    
    # Compute only if we need more zeros than are cached
    need = N - len(zeros)
    if need > 0:
        # Fetch the missing tail and append
        tail = np.asarray(fetch_zeros(d, len(zeros), need, lcalc_path), dtype=float)
        zeros = np.concatenate([zeros, tail])
        np.save(zfile, zeros)
    
    return zeros[:N]    

def cached_intervals(d: int, N: int, eps: float, lcalc_path: str, cache_root: str="data") -> np.ndarray:
    """
    Return an (N, 2) array of symmetric intervals [γ - eps, γ + eps] around zeros of χ_d
    Use a cache for efficiency. Recomputes if eps or N has changed 

    Input:  d (int): Discriminant of the Dirichlet character
            N (int): Number of zeros (intervals) to return
            eps (float): Half-width of the interval around each zero
            lcalc_path (str): Path to the lcalc executable
            cache_root (str): Root folder for cached data (default: data)

    Output: Array of shape (N, 2), where each row is [γ - eps, γ + eps] (np.ndarray)
    """
    # Construct the path to intervals.npy file
    parent_dir = "positive_d" if d > 0 else "negative_d"
    dir_path   = os.path.join(cache_root, parent_dir, f"d_{d}")
    os.makedirs(dir_path, exist_ok=True)
    ifile = os.path.join(dir_path, "intervals.npy")

    # Get N zeros
    zeros = cached_zeros(d, N, lcalc_path, cache_root)

    # Try to use cached intervals if valid
    if os.path.exists(ifile):
        intervals = np.load(ifile)
        # Check if cache is consistent with N and eps
        if (len(intervals) >= N) and np.allclose(intervals[0, 1] - intervals[0,0], 2 * eps):
            return intervals[:N]
    
    # Otherwise, regenerate intervals
    intervals = np.column_stack((zeros - eps, zeros + eps))
    np.save(ifile, intervals)
    return intervals[:N]
    