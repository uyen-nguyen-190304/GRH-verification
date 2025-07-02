"""
driver.py

Command-line entry point for verifying the Generalized Riemann Hypothesis (GRH) for 
quadratic Dirichlet L-functions using first logarithmic derivative 

This script parses user arguments, filters to fundamental discriminants, chooses a 
verification height η (either user-provided or based on the first zero), and then
calls base_case_verify for each d. Results are logged to stdout and appended to a CSV summary

Arguments:
    -d, --discriminant      Single discriminant to test (mutually exclusive with d-min/d-max)
    --d-min, --d-max        Inclusive range of discriminants to test
    -eta, --height          Optional window half-width η
                            If not provided, set to (first zero ordinate + 2*ε)
    -k, --power             Which logarithmic derivative to use (base case k = 1)
    -K, --upper-limit       Truncation parameter K for χ/Λ arrays (default 1e5)
    -eps, --epsilon         Half-width epsilon for zero intervals (default 1e-6)
    -config, --config-file  Path to JSON config containing {"lcalc_path": "..."}
    -data, --data-dir       Directory for caching zeros, intervals, χ, Λ (default "data")
    -output, --output-dir   Directory for output and logs (default "results")

Usage:
    python driver.py --d-min -1000 --d-max 1000

Outputs:
    results/summary.csv     CSV with columns [d, eta, N_needed]
    results/errors.log      Any runtime errors per discriminant
    data/...                Computed values of zeros, intervals, χ, Λ for each d
"""

import csv
import json
import argparse
from typing import List
from pathlib import Path

from grhverify.utils.discriminant import is_fundamental_discriminant
from grhverify.utils.generate_zeros import compute_zeros
from grhverify.base_case import base_case_verify

# =========================== BASE CASE ENTRY ===========================

def main() -> None:
    """
    Parse command-line arguments, prepare I/O, and run base_case_verify
    for each requested fundamental discriminant
    """
    parser = argparse.ArgumentParser(description="RH Verification")

    # Discriminant selection
    parser.add_argument("-d", "--discriminant", type=int, help="Discriminant of interest")
    parser.add_argument("--d-min", type=int, help="Minimum discriminant of interest")
    parser.add_argument("--d-max", type=int, help="Maximum discriminant of interest")    
    
    # Verification parameters
    parser.add_argument("-eta", "--height", type=float, help="Width of the window to verify the RH")
    parser.add_argument("-k", "--power", type=int, default=1, help="K-th logarithmic derivative L'/L used in verification")
    parser.add_argument("-K", "--upper-limit", type=int, default=10**5, help="Upper limit for lambda/kronecker array")
    parser.add_argument("-eps", "--epsilon", type=float, default=1e-6, help="Half-width for the zero intervals")
    
    # Paths for tools & I/O
    # parser.add_argument("-config", "--config-file", type=str, default="example_config.json", help="Path to config.json with lcalc_path")
    parser.add_argument("-config", "--config-file", type=str, default="local_config.json", help="Path to config.json with lcalc_path")
    parser.add_argument("-data", "--data-dir", type=str, default="data", help="Directory for data/input/cache files")
    parser.add_argument("-output", "--output-dir", type=str, default="results", help="Directory for output file")
    
    args = parser.parse_args()

    # Validate discriminant input
    if args.discriminant is not None:
        # Single-d mode cannot mix with range
        if args.d_min or args.d_max:
            raise ValueError("Cannot provide both --discriminant and --d-min/d-max")
        d_list: List[int] = [args.discriminant]
    elif (args.d_min is not None) and (args.d_max is not None):
        if args.d_min > args.d_max:
            raise ValueError("--d-min must be less than or equal to --d-max")
        d_list = list(range(args.d_min, args.d_max + 1))
    else:
        raise ValueError("Must provide either --discriminant or both --d-min and --d-max")

    # Load the config file (lcalc path)
    config_path = Path(args.config_file).expanduser()
    try:
        lcalc_path = Path(json.loads(config_path.read_text())["lcalc_path"]).expanduser()
    except Exception as err:
        raise RuntimeError(f"Failed to read lcalc_path from {config_path}") from err

    # Prepare I/O directories
    output_dir = Path(args.output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path   = output_dir / "errors.log"
    data_dir   = Path(args.data_dir).expanduser()
    data_dir.mkdir(parents=True, exist_ok=True)

    # Open summary CSV for append
    summary_path = output_dir / "summary.csv"
    new_file = not summary_path.exists()
    summary_fh = summary_path.open("a", newline="")
    writer = csv.writer(summary_fh)
    if new_file:
        # write header on first run
        writer.writerow(["d", "eta", "N_needed"])

    # Direct the RH verification based on the power k 
    k = args.power
    if k == 1:
        # Direct to base case verification
        for d in d_list:
            # Skip non-fundamental discriminants silently
            if not is_fundamental_discriminant(d):
                continue

            # Determine height η: user input or based on first zero + padding
            eta: float
            if args.height is not None:
                eta = args.height
            else:
                try:
                    first_zero = float(compute_zeros(d, 1, lcalc_path)[0])
                    eta = first_zero + 2 * args.epsilon     # Padding so that η is larger than the upper interval bound
                except Exception as err:
                    raise RuntimeError(f"Fail to compute the first zero ordinate to use at height eta for GRH verification")

            # Call the base case verification function
            success, eta, N_used = base_case_verify(
                d=d,
                K=args.upper_limit,
                eta=eta,
                eps=args.epsilon,
                lcalc_path=lcalc_path,
                data_dir=data_dir,
                log_path=log_path
            )

            # Human‐readable console output
            writer.writerow([d, eta, N_used])
            summary_fh.flush()

            # Print to console
            if success:
                print(f"For d = {d}, {N_used} zeros are needed to verify the RH up to height eta = {eta}")
            else:
                print(f"Fail to verify the RH with discriminant d = {d} up to height eta = {eta}")

    elif (k % 2 == 0 and k >= 2): 
        # Later implementation for higher power
        pass

    else:
        raise ValueError("k must be either 1 or an even integer greater than or equal to 2")


if __name__ == "__main__":
    main()
    