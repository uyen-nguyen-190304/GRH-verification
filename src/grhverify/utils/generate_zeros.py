import subprocess
import numpy as np
from pathlib import Path
from typing import List, Optional

def compute_zeros(d: int, N: int, lcalc_path: str | Path) -> List[float]:
    """
    Retrieve the first N positive ordinates of the Dirichlet L-function associated
    with discriminant d using lcalc

    Input:  d (int): Fundamental discriminant of the Dirichlet character χ_d
            N (int): Number of non-trivial zeros to compute
            lcalc_path (str | Path): Path to the lcalc executable

    Output: List of the first N non-trivial zeros ordinates (float)
    """
    # Input validation
    if not isinstance(d, int):
        raise TypeError("Discriminant d must be an integer")
    if not isinstance(N, int) or N <= 0:
        raise ValueError("N must be a positive integer")

    # Get the path to lcalc executable
    exe = Path(lcalc_path).expanduser()
    if not exe.is_file():
        raise FileNotFoundError(f"lcalc executable not found at {exe}")

    # Prepare the command to run lcalc
    cmd = [
        lcalc_path,
        "-z", str(N),                   # Number of zeros to compute
    ]

    # Add twist quadratic character if d is not primitive (d != -1, 1)
    if abs(d) != 1:
        cmd.append("--twist-quadratic")
    
    # Now, add the start and finish flags
    cmd += [
        "--start", str(d),              # Discriminant d
        "--finish", str(d)
    ]

    # Run the command and capture the output
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
    except subprocess.CalledProcessError as err:
        raise RuntimeError(
            f"lcalc failed with status {err.returncode}\n"
            f"stdout: {err.stdout}\n"
            f"stderr: {err.stderr}"
        ) from err

    # Extract the ordinates from the output
    zeros: list[float] = []
    for line in res.stdout.splitlines():
        try:
            zeros.append(float(line.split()[-1]))
        except ValueError:
            continue                    # Skip lines that do not contain valid floats

    # Check if we have enough zeros
    if len(zeros) < N:
        raise RuntimeError(
            f"Excepted {N} zeros but lcalc returned {len(zeros)}"
        )   

    return zeros[:N]  # Return only the first N zeros

def write_zeros(d: int, zeros: List[float], data_dir: str | Path) -> None:
    """
    Save the zeros ordinates to the zeros.txt file in the corresponding directory

    Input:  d (int): Discriminant of the Dirichlet character 
            zeros (List[float]): List of zeros ordinates
            data_dir (str): Path to the data directory

    Output: None
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
    Return an (N, 2) array of symmetric intervals [gamma - eps, gamma + eps] around zeros 
    of the Dirichlet L-function associated with discriminant d

    Input:  d (int): Discriminant of the Dirichlet character
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
    Save the intervals [gamma - eps, gamma + eps] to the intervals.txt file in the corresponding directory

    Input:  d (int): Discriminant of the Dirichlet character
            intervals (np.ndarray): Array of shape (N, 2) containing the intervals
            data_dir (str | Path): Path to the data directory

    Output: None
    """
    # Ensure the storing directory exists
    target = Path(data_dir).expanduser() / ("positive_d" if d > 0 else "negative_d") / f"d_{d}"
    target.mkdir(parents=True, exist_ok=True)

    # Save the intervals to the intervals.txt file
    txt_path = target / "intervals.txt"
    with open(txt_path, "w") as f:
        for gamma_minus, gamma_plus in intervals:
            f.write(f"{gamma_minus} {gamma_plus}\n")
