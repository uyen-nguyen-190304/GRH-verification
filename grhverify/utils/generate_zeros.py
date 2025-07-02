"""
generate_zeros.py

Functions to fetch and store the nontrivial zeros of quadratic Dirichlet L-functions
via the lcalc command-line tool, and for building small symmetric intervals around each zeros

Functions:
    - compute_zeros(d, N, lcalc_path): Use lcalc to compute the first N positive ordinates (imaginary parts) of the nontrivial zeros of L(s, χ_d)
    - write_zeros(d, zeros, data_dir): Write the list of zero ordinates to data_dir/{positive_d|negative_d}/d_{d}/zeros.txt
    - compute_intervals(d, N, eps, lcalc_path, zeros): Given either a precalculated list of zeros or by invoking compute_zeros, build an (N,2) numpy array of [gamma - eps, gamma + eps] rows
    - write_intervals(d, intervals, data_dir): Write the array of intervals to data_dir/{positive_d|negative_d}/d_{d}/intervals.txt

Notes
-----
- Path to lcalc executable must be provided inside the config.json file
- Work for any discriminant d; however, if d is not fundamental, lcalc will return nothing
"""

import subprocess
import numpy as np
from pathlib import Path
from typing import List, Optional


def compute_zeros(d: int, N: int, lcalc_path: str | Path) -> List[float]:
    """
    Purpose:
        Retrieve the first N positive ordinates of the Dirichlet L-function for χ_d
    Input:  
        d (int): Fundamental discriminant of the Dirichlet character χ_d
        N (int): Number of non-trivial zeros to compute
        lcalc_path (str | Path): Path to the lcalc executable
    Return:
        List[float] of the first N non-trivial zero ordinates
    """
    # Input validation
    if not isinstance(d, int):
        raise TypeError("Discriminant d must be an integer")
    if not isinstance(N, int) or N <= 0:
        raise ValueError("N must be a positive integer")

    # Locate the lcalc executable
    exe = Path(lcalc_path).expanduser()
    if not exe.is_file():
        raise FileNotFoundError(f"lcalc executable not found at {exe}")

    # Build the command to compute zeros
    cmd = [
        lcalc_path,
        "-z", str(N),                   # Number of zeros to compute
    ]

    # For nontrivial discriminants (|d| != 1), twist by the quadratic character
    if abs(d) != 1:
        cmd.append("--twist-quadratic")
    
    # Restrict to the single discriminant d
    cmd += [
        "--start", str(d),              # Discriminant d
        "--finish", str(d)
    ]

    # Execute lcalc and capture stdout
    try:
        res = subprocess.run(
            cmd, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE, 
            check=True, 
            text=True
        )
    except subprocess.CalledProcessError as err:        # Wrap and propagate any lcalc error
        # Wrap and propagate any lcalc error
        raise RuntimeError(
            f"lcalc failed with status {err.returncode}\n"
            f"stdout: {err.stdout}\nstderr: {err.stderr}"
        ) from err

    # Parse the output lines for floats
    zeros: list[float] = []
    for line in res.stdout.splitlines():
        try:
            # Last token is expected to be the ordinate
            zeros.append(float(line.split()[-1]))
        except ValueError:
            # Skip any header or footer lines that don't parse
            continue                   

    # Ensure we got N zeros back
    if len(zeros) < N:
        raise RuntimeError(
            f"Excepted {N} zeros but lcalc returned {len(zeros)}"
        )   

    return zeros[:N]  # Return only the first N zeros


def write_zeros(d: int, zeros: List[float], data_dir: str | Path) -> None:
    """
    Purpose:
        Save the zeros ordinates to the zeros.txt file in the corresponding directory
    Input:  
        d (int): Discriminant of the Dirichlet character 
        zeros (List[float]): List of zeros ordinates
        data_dir (str): Path to the data directory
    Return:
        None
    """
    # Ensure the storing directory exists
    target = Path(data_dir).expanduser() / ("positive_d" if d > 0 else "negative_d") / f"d_{d}"
    target.mkdir(parents=True, exist_ok=True)
    
    # Save the zeros to the zeros.txt file 
    txt_path = target / "zeros.txt"
    with open(txt_path, "w") as f:
        f.writelines(f"{gamma}\n" for gamma in zeros)


def compute_intervals(d: int, N: int, eps: float, lcalc_path: str, zeros: Optional[List[float]] = None) -> np.ndarray:
    """
    Purpose:
        Return an (N, 2) array of symmetric intervals [gamma - eps, gamma + eps] around zeros 
        of the Dirichlet L-function associated with discriminant d
    Input:  
        d (int): Discriminant of the Dirichlet character
        N (int): Number of zeros (intervals) to return
        eps (float): Half-width of the interval around each zero
        lcalc_path (str): Path to the lcalc executable
        zeros (Optional[List[float]]): Precomputed zeros, if available

    Output: Array of shape (N, 2), where each row is [gamma - eps, gamma + eps] (np.ndarray)
    """
    # If zeros are not provided or there are not enough, compute them
    if zeros is None or len(zeros) < N:
        zeros = compute_zeros(d, N, lcalc_path)

    # Convert zeros to a numpy array and ensure it has the correct type
    zeros = np.asarray(zeros[:N], dtype=float)

    # Construct the intervals around each zero and return them
    return np.column_stack((zeros - eps, zeros + eps))


def write_intervals(d: int, intervals: np.ndarray, data_dir: str) -> None:
    """
    Purpose:
        Save the intervals [gamma - eps, gamma + eps] to the intervals.txt file in the corresponding directory
    Input:
        d (int): Discriminant of the Dirichlet character
        intervals (np.ndarray): Array of shape (N, 2) containing the intervals
        data_dir (str | Path): Path to the data directory
    Output: 
        None
    """
    # Ensure the storing directory exists
    target = Path(data_dir).expanduser() / ("positive_d" if d > 0 else "negative_d") / f"d_{d}"
    target.mkdir(parents=True, exist_ok=True)

    # Save the intervals to the intervals.txt file
    txt_path = target / "intervals.txt"
    with open(txt_path, "w") as f:
        for gamma_minus, gamma_plus in intervals:
            f.write(f"{gamma_minus} {gamma_plus}\n")
