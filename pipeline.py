import sys
import argparse
import subprocess

def compute_default_eta(zeros_file):
    """
    Compute the default eta value as the midpoint of the first two zeros

    Input:  zeros (list): List of computed zeros

    Output: float: Default eta value, or 50.0 if any errors occur/not enough zeros are available
    """
    try:
        with open(zeros_file, 'r') as f:
            lines = f.readlines()
            if len(lines) >= 2:
                gamma1 = float(lines[0].strip())
                gamma2 = float(lines[1].strip())
                return (gamma1 + gamma2) / 2.0
            else:
                # Not enough zeros to compute eta
                return 50.0

    except Exception as e:
        # If there's an error reading the file or parsing the zeros, return a default value
        return 50.0     


def main():
    """
    Main function to run the GRH verification pipeline
    """
    parser = argparse.ArgumentParser(description="GRH Verification Pipeline")
    parser.add_argument("-d", "--discriminant", dest="d", type=int, required=True, help="Discriminant of the Dirichlet L-function")
    parser.add_argument("-N", "--num-zeros", dest="N", type=int, default=100, help="Number of zeros to compute")
    parser.add_argument("-K", "--upper-limit", dest="K", type=int, default=20, help="Upper limit for von Mangoldt and Kronecker computation")
    parser.add_argument("-eta", type=float, help="Height parameter of interest eta")
    parser.add_argument("-err", "--error-window", dest="err", type=float, default=1e-8, help="Error window for the intervals")
    parser.add_argument("-config", "--config-file", dest="config", type=str, default="example_config.json", help="Path to config.json with lcalc_path")
    args = parser.parse_args()
    
    # Retrieve command line arguments
    d = args.d
    N = args.N
    K = args.K
    eta = args.eta
    err = args.err
    config_file = args.config

    try: 
        # 1. Compute zeros and intervals of Dirichlet L-functions using lcalc
        subprocess.run([
            "python3", "script/generate_zeros.py", 
            "-d", str(d),           # Discriminant
            "-N", str(N),           # Number of zeros
            "-err", str(err),       # Error window
            "-config", config_file  # Path to config.json
        ], check=True)

        # 2. Compute von Mangoldt values
        subprocess.run([
            "python3", "script/compute_von_mangoldt.py",
            "-n", str(K)            # Upper limit for von Mangoldt
        ], check=True)

        # 3. Compute Kronecker symbols
        subprocess.run([
            "python3", "script/compute_kronecker.py",
            "-d", str(d),           # Discriminant
            "-n", str(K)            # Upper limit for Kronecker
        ], check=True)

        # Check if eta is provided
        if eta is None:
            # If eta is not provided, set it to be the middle point of the first two zeros
            eta = compute_default_eta("data/zeros.txt")

        # 4. Verify RH inequality (Corollary 1)
        subprocess.run([
            "src/rh_verify",
            str(d),                 # Discriminant
            str(eta),               # Height parameter of interest
            str(K),                 # Upper limit for von Mangoldt and Kronecker
        ], check=True)

    # If any subprocess fails, catch the error and print a message
    except subprocess.CalledProcessError as e:
        print(f"Error running the pipeline: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
