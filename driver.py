import csv
import json
import argparse
from typing import List
from pathlib import Path

from grhverify.base_case import base_case_verify
from grhverify.utils.generate_zeros import compute_zeros

# =========================== HELPER FUNCTIONS ===========================

def is_square_free(n: int) -> bool:
    """
    Return True if n is square-free, else False
    """
    n = abs(n)

    if n in (0, 1):
        return n == 1   # 0 -> False, 1 -> True
    
    if n % 4 == 0:
        return False    # Quick early rejection: divisibility be 4 can't be square-free
    
    p = 2
    while p * p <= n:
        if n % (p * p) == 0:
            return False    # Found a square divisor
        p +=1 if p == 2 else 2  # Skip even number > 2
    return True

def is_fundamental_discriminant(d: int) -> bool:
    """
    Check if an integer d is a fundamental discriminant 
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

# =========================== PIPELINE ENTRY ===========================

def main() -> None:
    """
    Main function to run the GRH verification process
    """
    parser = argparse.ArgumentParser(description="RH Verification")
    parser.add_argument("-d", "--discriminant", type=int, help="Discriminant of interest")
    parser.add_argument("--d-min", type=int, help="Minimum discriminant of interest")
    parser.add_argument("--d-max", type=int, help="Maximum discriminant of interest")    
    
    parser.add_argument("-eta", "--height", type=float, help="Width of the window to verify the RH")
    parser.add_argument("-k", "--power", type=int, default=1, help="K-th logarithmic derivative L'/L used in verification")
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

    # Prepare I/O dirs
    output_dir = Path(args.output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "errors.log"

    data_dir = Path(args.data_dir).expanduser()
    data_dir.mkdir(parents=True, exist_ok=True)

    summary_path = output_dir / "summary.csv"
    new_file = not summary_path.exists()
    summary_fh = summary_path.open("a", newline="")
    writer = csv.writer(summary_fh)
    if new_file:
        writer.writerow(["d", "N_needed"])

    # Direct the RH verification based on the power k 
    k = args.power
    if k == 1:
        # Direct to base case verification
        for d in d_list:
            # First, check if d is a fundamental discriminant
            if not is_fundamental_discriminant(d):
                continue

            # Compute the height eta if user not provided
            eta: float
            if args.height is not None:
                eta = args.height
            else:
                try:
                    eta = float(compute_zeros(d, 1, lcalc_path)[0])
                except Exception as err:
                    raise RuntimeError(f"Fail to compute the first zero ordinate to use at height eta for GRH verification")

            # Call the base case verification function
            success, N_used = base_case_verify(
                d=d,
                K=args.upper_limit,
                eta=eta,
                eps=args.epsilon,
                lcalc_path=lcalc_path,
                data_dir=data_dir,
                log_path=log_path
            )

            # Log result to CSV
            writer.writerow([d, N_used])
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