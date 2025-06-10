#!usr/bin/env python3

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
            lcalc_path (str): Path to the lcalc executable

    Output: List of zeros (float)
    """
    # Validate inputs
    if not isinstance(N, int) or N <= 0:
        print("Error: Number of zeros N must be a positive integer")
        sys.exit(1)

    # Expand path to lcalc executable
    lcalc_path = os.path.expanduser(lcalc_path)
    if not os.path.isfile(lcalc_path):
        print(f"Error: lcalc executable not found at: {lcalc_path}")
        sys.exit(1)

    # Prepare the lcalc command to compute zeros
    cmd = [
        lcalc_path, 
        "-z", str(N),           # Number of zeros to compute
        "--twist-quadratic",    # Twist the L-function by a quadratic character
        "--start", str(d),      # Discriminant d
        "--finish", str(d)  
    ]

    # Calling the lcalc executable 
    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
        output = result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running lcalc: \nSTDOUT: {e.stdout}\nSTDERR: {e.stderr}")
        sys.exit(1)
    
    # Parse zeros from the output
    zeros = []
    for line in output.splitlines():
        line = line.strip()
        if re.match(r'^\d+\s+[\d\.]+', line):
            parts = line.split()
            if len(parts) >= 2:
                gamma = float(parts[1])     # Imaginary part of the zero
                zeros.append(gamma)

    return zeros


def write_zeros(zeros: list, filename: str):
    """
    Write the computed zeros to a file
    
    Input:  zeros (list): List of zeros to write
            filename (str): Name of the file to write the zeros to

    Output: None
    """
    dir_name = os.path.dirname(filename)
    os.makedirs(dir_name, exist_ok=True)  # Create directory if it doesn't exist

    with open(filename, 'w') as f:
        for zero in zeros:
            f.write(f"{zero}\n")


def compute_intervals(zeros: list, err: float = 1e-8):
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
    dir_name = os.path.dirname(filename)
    os.makedirs(dir_name, exist_ok=True)  # Create directory if it doesn't exist

    with open(filename, 'w') as f:
        for gamma_minus, gamma_plus in intervals:
            f.write(f"{gamma_minus} {gamma_plus}\n")


def main():    
    """
    Pipeline to compute zeros of Dirichlet L-functions and write them to files
    1. Compute zeros using lcalc
    2. Write zeros to a file
    3. Compute intervals for the zeros
    4. Write intervals to a file
    """
    parser = argparse.ArgumentParser(description="Compute zeros of Dirichlet L-functions and intervals")
    parser.add_argument('-d', '--discriminant', type=int, required=True, help='Discriminant of the Dirichlet L-function')
    parser.add_argument('-N', '--num-zeros', type=int, default=20, help='Number of zeros to compute (default: 20)')
    parser.add_argument('-err', '--error-window', type=float, default=1e-8, help='Error window for intervals (default: 1e-8)')
    parser.add_argument('-config', '--config-file', type=str, default='config.json', help='Path to config file (default: config.json)')
    parser.add_argument('-output', '--output-dir', type=str, default=None, help='Output directory to save results')
    args = parser.parse_args()

    # Retrieve command line arguments
    d = args.discriminant
    N = args.num_zeros
    err = args.error_window
    lcalc_path = args.config_file  
    output_dir = args.output_dir

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

    # Determine output directory based on sign of discriminant
    if output_dir is None:
        parent_dir = "data/positive_d" if d > 0 else "data/negative_d"
        output_dir = os.path.join(parent_dir, f"d_{d}")
    os.makedirs(output_dir, exist_ok=True)  # Create directory if it doesn't exist

    # Compute the zeros of the Dirichlet L-function
    zeros = compute_zeros(d, N, lcalc_path)
    zeros_file = os.path.join(output_dir, "zeros.txt")
    write_zeros(zeros, zeros_file)

    # Compute intervals for the zeros
    intervals = compute_intervals(zeros, err)
    intervals_file = os.path.join(output_dir, "intervals.txt")
    write_intervals(intervals, intervals_file)


if __name__ == "__main__":
    main()
