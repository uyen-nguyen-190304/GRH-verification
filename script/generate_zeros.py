import os
import re
import sys
import json
import argparse
import subprocess

def compute_zeros(d: int, N: int, lcalc_path: str):
    """
    Compute the first N zeros of the Dirichlet L-function with discriminant d

    Input:  d (int): Discriminant of the Dirichlet L-function
            N (int): Number of zeros to compute

    Output: List of zeros (float)
    """
    # Validate inputs
    if not isinstance(d, int) or d <= 0:
        print("Error: Discriminant d must be a positive integer.")
        sys.exit(1)
    if not isinstance(N, int) or N <= 0:
        print("Error: Number of zeros N must be a positive integer.")
        sys.exit(1)

    # Get the path to the lcalc executable
    lcalc_path = os.path.expanduser(lcalc_path)

    if not os.path.isfile(lcalc_path):
        print(f"Error: lcalc executable not found at: {lcalc_path}")
        sys.exit(1)

    # TODO: Double check this command line
    # Command to run lcalc 
    cmd = [
        lcalc_path, 
        "-z", str(N),           # Number of zeros to compute
        "--twist-quadratic",    # Twist the L-function by a quadratic character
        "--start", str(d),      # Discriminant d
        "--finish", str(d)  
    ]

    # Calling the lcalc executable -> compute zeros
    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
        output = result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running lcalc: \nSTDOUT: {e.stdout}\nSTDERR: {e.stderr}")
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

def write_zeros(zeros: list, filename: str):
    """
    Write the computed zeros to a file
    
    Input:  zeros (list): List of zeros to write
            filename (str): Name of the file to write the zeros to

    Output: None
    """
    with open(filename, 'w') as f:
        for zero in zeros:
            f.write(f"{zero}\n")

def compute_intervals(zeros: list, err: float = "1e-8"):
    """
    Compute intervals for the imaginary part of the zeros based on the error window

    Input:  zeros (list): List of computed zeros
            err (float): Error window for the intervals, default of 1e-8

    Output: List of tuples representing intervals for each zero
    """
    intervals = []
    for zero in zeros:
        gamma_minus = zero - err
        gamma_plus  = zero + err
        intervals.append((gamma_minus, gamma_plus))
    return intervals

def write_intervals(intervals: list, filename: str):
    """
    Write the computed intervals to a file

    Input:  intervals (list): List of intervals to write
            filename (str): Name of the file to write the intervals to

    Output: None
    """
    with open(filename, 'w') as f:
        for gamma_minus, gamma_plus in intervals:
            f.write(f"{gamma_minus} {gamma_plus}\n")

def main():    
    """
    Pipeline to compute zeros of Dirichlet L-functions and write them to files.
    1. Compute zeros using lcalc
    2. Write zeros to a file
    3. Compute intervals for the zeros
    4. Write intervals to a file
    """
    parser = argparse.ArgumentParser(description="Compute zeros of Dirichlet L-functions and intervals")
    parser.add_argument('-d', '--discriminant', type=int, required=True, help='Discriminant of the Dirichlet L-function')
    parser.add_argument('-N', '--num-zeros', type=int, default=50, help='Number of zeros to compute (default: 50)')
    parser.add_argument('-err', '--error-window', type=float, default=1e-8, help='Error window for intervals (default: 1e-8)')
    parser.add_argument('-config', '--config-file', type=str, default='../config.json', help='Path to config file (default: ../config.json)')
    args = parser.parse_args()

    # Retrieve command line arguments
    d = args.discriminant
    N = args.num_zeros
    err = args.error_window
    lcalc_path = args.config_file  

    # Load lcalc_path from JSON config file
    try:
        with open(args.config_file, 'r') as f:
            config = json.load(f)
            lcalc_path = config.get('lcalc_path')
            if not lcalc_path:
                print("Error: 'lcalc_path' not found in config file.")
                exit(1)
    except FileNotFoundError:
        print(f"Error: Config file {args.config_file} not found.")
        exit(1)

    # Compute the zeros of the Dirichlet L-function
    zeros = compute_zeros(d, N, lcalc_path)
    write_zeros(zeros, filename="../data/zeros.txt")

    # Compute intervals for the zeros
    intervals = compute_intervals(zeros, err)
    write_intervals(intervals, filename="../data/intervals.txt")

if __name__ == "__main__":
    main()
