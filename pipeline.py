import json
import argparse
import subprocess

def main():
    parser = argparse.ArgumentParser(description="GRH Verification Pipeline")    
    parser.add_argument("-d", "--discriminant", type=int, required=True, help="Discriminant of the Dirichlet L-function")
    parser.add_argument("-N", type=int, default=100, help="Number of zeros to compute")
    parser.add_argument("-K", type=int, default=20, help="")
    parser.add_argument("-eta", type=float, default= , help="Height parameter of interest eta")
    parser.add_argument("-err","--error_window", type=float, default=1e-8, help="Error window for the intervals")
    
    args = parser.parse_args()

    # Retrieve command line arguments
    d = args.d
    N = args.N
    K = args.K
    eta = args.eta
    error_window = args.error_window

    try: 
        # 1. Compute zeros and intervals of Dirichlet L-functions using lcalc
        subprocess.run([
            "python3", "script/generate_zeros.py", 
            str(d),         # Discriminant
            str(N),         # Number of zeros
            "config.json",  # Configuration file -> lcalc path
            error_window    # Error window
            ]
        )

        # 



    except subprocess.CalledProcessError as e:
        print(f"Error running the pipeline: {e}")
        exit(1)

if __name__ == "__main__":
    main()