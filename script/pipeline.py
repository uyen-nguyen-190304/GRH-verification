import csv
import sys
import json
import math
import argparse
import subprocess
import numpy as np
from typing import Tuple
from pathlib import Path

import grhverify
from .zeros import cached_zeros, cached_intervals

def is_square_free(n: int) -> bool:
    """
    Helper function to determine whether a number is square-free 

    Input:  n (int): Number to check

    Output: True if square-free, False otherwise
    """
    n = abs(n)

    if n in (0, 1):
        return n == 1   # 0 -> False, 1 -> True
    
    if n % 4 == 0:
        return False    # Quick early rejection: divisibly be 4 can't be square-free
    
    p = 2
    while p * p <= n:
        if n % (p * p) == 0:
            return False    # Found a square divisor
        p +=1 if p == 2 else 2  # Skip event number > 2
    return True

def is_fundamental_discriminant(d: int) -> bool:
    """
    Check if an integer d is a fundamental discriminant 

    Input:  d (int): The discriminant to check

    Output: True if d is a fundamental discriminant, False otherwise
    """
    # Early check: 0 is not a fundamental discriminant
    if d == 0:
        return False
    
    # Case where d ≡ 1 mod 4: check if abs(d) is square-free 
    if d % 4 == 1:
        return is_square_free(d)

    # Case where d ≡ 0 mod 4: check if d/4 is square-free and ≡ 2 or 3 mod 4
    if d % 4 == 0:
        q = d // 4
        return is_square_free(q) and q % 4 in (2, 3)
    return False

def load_lambda(K: int, data_dir: Path) -> np.ndarray:
    """
    Load von Mangoldt function value Λ(k) up to k = K
    If not cached, trigger compute_lambda.py script to compute it

    Input:  K (int): Upper limit for k
            data_dir(Path): Base data directory

    Output: Array of Λ(k) values (np.ndarray)
    """
    # Build path to cached Λ(k) values
    f = data_dir / "cache" / f"von_mangoldt_{K}.npy"

    # If Λ(k) not cached, run compute script to generate it
    if not f.exists():
        subprocess.check_call([
            "python",
            "script/compute_lambda.py",
            "-K", str(K),
            "--data-dir", str(data_dir)
        ])

    # Load the numpy array
    return np.load(f)

def load_kronecker(d: int, K: int, data_dir: Path) -> np.ndarray:
    """
    Load Kronecker symbol values χ_d(k) up to k = K
    If not cached, trigger compute_kronecker.py script to compute it 

    Input:  d (int): Discriminant
            K (int): Upper limit for k
            data_dir (Path): Base data directory

    Output: Array of χ_d(k) values (np.ndarray)
    """
    # Build path to cached χ_d(k) values
    f = data_dir / "cache" / ("positive_d" if d > 0 else "negative_d") / f"d_{d}" / f"chi_{d}_K{K}.npy"
    
    # If χ_d(k) values don't exist, run script to compute and cache them
    if not f.exists():
        subprocess.check_call([
            "python",
            "script/compute_kronecker.py",
            "-d", str(d),
            "-K", str(K),
            "--data-dir", str(data_dir)
        ])

    # Load the numpy array
    return np.load(f)

def verify_until_success(d: int, K: int, eps: float, eta: float, lcalc_path: str, data_dir: Path, log_path: Path, chunk: int = 20) -> Tuple[bool, int]:
    """
    Verify RH for the Dirichlet L-function associated with discriminant d
    Keep increasing N (number of zeros used) until success or error
    Incrementally add zeros in blocks of size chunks until the equality is statisfied (for cached purposes)
    
    Input:  d (int): Discriminant
            K (int): Upper bound for lambda and Kronecker array
            eps (float): Half-width of interval used for zero approximation
            eta (float): Verification height
            lcalc_path (str): Path to lcalc executable
            data_dir (Path): Directory containing cached files
            log_path (Path): Path to error log
            chunk (int): Number of zeros/intervals to load at once
    
    Output: A tupe of (success (bool), number of zeros used (int))
    """
    # Construct the cache root
    cache_root = str(data_dir / "cache")

    # Load the cached values of lambda and kronecker symbol
    lambda_arr = load_lambda(K, data_dir)
    chi_arr = load_kronecker(d, K, data_dir)
    
    # ----------------------- RHS -----------------------
    
    # Compute RHS constant based on discriminant sign 
    if d < 0:
        rhs_const = 0.5 * math.log(abs(d)*math.e**2 /
                                   (4*math.pi*math.e**grhverify.EULER_CONSTANT))
    else:
        rhs_const = 0.5 * math.log(abs(d) /
                                   (math.pi*math.e**grhverify.EULER_CONSTANT))

    # Add logarithmic derivative contribution to RHS
    rhs = rhs_const + grhverify.log_derivative(chi_arr, lambda_arr, K)

    # ----------------------- LHS -----------------------
    
    # Initialize LHS with constant iota term
    lhs    = 2.0 * grhverify.iota(eta)

    # Start loading zeros incrementally
    N_used = 0 
    start  = 0

    try:
        while True:
            # Load chunk number of new intervals
            need = start + chunk
            intervals_np = cached_intervals(d, need, eps, lcalc_path, cache_root)

            if len(intervals_np) <= start:
                break
            
            intervals = intervals_np[start:]    # Only use new batch
            for gamma_minus, gamma_plus in map(tuple, intervals):
                lhs += grhverify.zero_contribution(gamma_minus, gamma_plus)
                N_used += 1
                if lhs > rhs:
                    raise StopIteration
            start += chunk
        
        # Loop exhaust without RH verified
        success = False

    except StopIteration:
        success = True

    except Exception as e:
        # Log error to file if computation fails for this d
        with open(log_path, "a") as log:
            log.write(f"Error: d = {d}, N = {N_used}, reason = {repr(e)}\n")
        success = False

    # Write zeros and intervals to .txt file
    try: 
        # Grab *exactly* the zeros/intervals we actually used
        zeros_np = cached_zeros(d, N_used, lcalc_path, cache_root)
        intervals_np  = cached_intervals(d, N_used, eps, lcalc_path, cache_root)

        # Path to correct directory
        zeros_txt = data_dir / ("positive_d" if d > 0 else "negative_d") / f"d_{d}" / f"zeros.txt"
        intervals_txt  = data_dir / ("positive_d" if d > 0 else "negative_d") / f"d_{d}" / f"intervals.txt"

        # Save the zeros and intervals used
        np.savetxt(zeros_txt, zeros_np, fmt="%.23f")
        np.savetxt(intervals_txt, intervals_np, fmt="%.23f")

    except Exception as save_err:
        with open(log_path, "a") as log:
            log.write(f"Warning: d = {d}, could not save zeros.txt and intervals.txt: {repr(save_err)}\n")

    return success, N_used

def main():
    """
    NOTE: I dont like how I still put config file right now since the /builds directory is located at different place
    TODO: A /tools directory? Or some better way for lcalc_path here

    Main entry point. Parse command-line arguments, load config, and verify RH
    for all target discrimants, and write summary output
    """
    parser = argparse.ArgumentParser(description="RH Verification")
    parser.add_argument("-d", "--discriminant", type=int, help="Discriminant of interest")
    parser.add_argument("--d-min", type=int, help="Minimum discriminant of interest")
    parser.add_argument("--d-max", type=int, help="Maximum discriminant of interest")    
    parser.add_argument("-eta", type=float, help="Width of the window to verify the RH")
    parser.add_argument("-K", "--upper-limit", type=int, default=10**5, help="Upper limit for lambda/kronecker array")
    parser.add_argument("-eps", "--epsilon", type=float, default=1e-6, help="Half-width for the zero intervals")
    # parser.add_argument("-config", "--config-file", type=str, default="example_config.json", help="Path to config.json with lcalc_path")
    parser.add_argument("-config", "--config-file", type=str, default="local_config.json", help="Path to config.json with lcalc_path")
    parser.add_argument("-data", "--data-dir", type=str, default="data", help="Directory for data/input/cache files")
    parser.add_argument("-output", "--output-dir", type=str, default="results", help="Directory for output file")
    
    args = parser.parse_args()

    # Validate discriminant input
    if args.discriminant is not None:
        if args.d_min or args.d_max:
            sys.exit("Error: Cannot provide both --discriminant and --d-min/d-max")
        d_list = [args.discriminant]
    elif (args.d_min is not None) and (args.d_max is not None):
        if args.d_min > args.d_max:
            sys.exit("Error: --d-min must be less than or equal to --d-max")
        d_list = list(range(args.d_min, args.d_max + 1))
    else:
        sys.exit("Error: Must provide either --discriminant or both --d-min and --d-max")

    # Load lcalc path from JSON config
    try:
        with open(args.config_file, "r") as f:
            config = json.load(f)
            lcalc_path = config.get("lcalc_path")
            if lcalc_path is None:
                raise ValueError(f"Error: Missing lcalc path in {args.config_file}.json")
    except Exception as e:
        sys.exit(f"Error: Failed to load config file: {e}")

    # Create output directory if missing
    data_dir = Path(args.data_dir)
    cache_root = str(data_dir / "cache")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    result_path = output_dir / "summary.csv"
    log_path = output_dir / "errors.log"

    with open(result_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["d", "N_needed"])
        
        for d in d_list:
            if not is_fundamental_discriminant(d):
                continue    # Skip if not a fundamental discriminant

            # Choose eta: use input if provided, otherwise use first zero height
            eta = args.eta
            if eta is None:
                try:
                    zero_list = cached_zeros(d, 1, lcalc_path, cache_root)
                    eta = zero_list[0]      # First non-trivial zero ordinate
                except Exception as e:
                    eta = 50.0     # Fallback value

            # Attempt to verify RH for this discriminant
            success, n_required = verify_until_success(d, args.upper_limit, args.epsilon, eta, lcalc_path, data_dir, log_path)
        
            # Log result to CSV
            writer.writerow([d, n_required if success else float('nan')])
            f.flush()

            # Print to console
            if success:
                print(f"For d = {d}, {n_required} zeros are needed to verify the RH up to height {eta}")
            else:
                print(f"Fail to verify the RH with discriminant d = {d}")
                      
if __name__ == "__main__":
    main()
