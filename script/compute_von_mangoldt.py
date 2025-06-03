import os
import math
import argparse
from sage.all import *     

def mangoldt_lambda(k: int) -> float:
    """
    Compute the von Mangoldt function value for a given integer k
    
    Input:  k (int): Integer for which to compute the von Mangoldt function

    Output: float: Value of the von Mangoldt function at k 
    """
    f = factor(k)
    if len(f) == 1:
        p, m = f[0]
        if p.is_prime():
            return math.log(p)
    return 0.0

def compute_von_mangoldt(n: int) -> list:
    """
    Compute the von Mangoldt function for integers k from [1..n]

    Input:  n (int): Upper limit for k

    Output: list of floats: von Mangoldt values 
    """
    von_mangoldt_list = []
    for k in range(1, n + 1):
        von_mangoldt_list.append(mangoldt_lambda(k))
    return von_mangoldt_list

def write_von_mangoldt(von_mangoldt_list: list, filename: str):
    """
    Write the computed von Mangoldt values to a file
    
    Input:  von_mangoldt_list (list): List of von Mangoldt values to write
            filename (str): Name of the file to write the values to

    Output: None
    """
    dir_name = os.path.dirname(filename)
    if dir_name and not os.path.exists(dir_name):
        os.makedirs(dir_name, exist_ok=True)

    with open(filename, 'w') as f:
        for k, von_mangoldt_value in enumerate(von_mangoldt_list, start=1):
            f.write(f"{k} {von_mangoldt_value}\n")
            
def main():
    """
    Pipeline to compute von Mangoldt function values from [1..n] and write them to a file
    """
    parser = argparse.ArgumentParser(description="Compute von Mangoldt function values")
    parser.add_argument("-n", "--upper-limit", dest="n", type=int, required=True, help="Upper limit for k in von Mangoldt function") 
    args = parser.parse_args()

    # Retrieve the upper limit for k
    n = args.n

    # Compute the von Mangoldt function values
    von_mangoldt_list = compute_von_mangoldt(n)
    write_von_mangoldt(von_mangoldt_list, filename="data/von_mangoldt.txt")

if __name__ == "__main__":
    main()