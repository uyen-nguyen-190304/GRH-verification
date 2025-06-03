import os
import argparse
from sage.all import *

def compute_kronecker(d: int, n: int) -> list:
    """
    Compute the Kronecker symbol for integers k from [1..n] with respect to d

    Input:  d (int): Discriminant for the Kronecker symbol
            n (int): Upper limit for k

    Output: list of ints: Kronecker symbol values
    """
    kronecker_list = []
    for k in range(1, n + 1):
        kronecker_list.append(kronecker(d, k))
    return kronecker_list

def write_kronecker(kronecker_list: list, filename: str):
    """
    Write the computed Kronecker symbols to a file
    
    Input:  kronecker_list (list): List of Kronecker symbols to write
            filename (str): Name of the file to write the symbols to

    Output: None
    """
    dir_name = os.path.dirname(filename)
    if dir_name and not os.path.exists(dir_name):
        os.makedirs(dir_name, exist_ok=True)

    with open(filename, 'w') as f:
        for k, kronecker_value in enumerate(kronecker_list, start=1):
            f.write(f"{k} {kronecker_value}\n")

def main():
    """
    Pipeline to compute Kronecker symbols from [1..n] and write them to a file
    """
    parser = argparse.ArgumentParser(description="Compute Kronecker symbols for integers [1..n] with respect to discriminant d")
    parser.add_argument('-d', '--discriminant', dest="d", type=int, required=True, help='Discriminant for the Kronecker symbol')
    parser.add_argument('-n', '--upper_limit', dest="n", type=int, required=True, help='Upper limit for k in the Kronecker symbol')
    args = parser.parse_args()

    # Retrieve the discriminant and upper limit for k
    d = args.d
    n = args.n

    # Compute the Kronecker symbols
    knonecker_list = compute_kronecker(d, n)
    write_kronecker(knonecker_list, filename="data/kronecker.txt")

if __name__ == "__main__":
    main()