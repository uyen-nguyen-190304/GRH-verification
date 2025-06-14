import argparse
import numpy as np
from pathlib import Path
from sage.all import kronecker_symbol

def compute_kronecker(d: int, K: int) -> np.ndarray:
    """
    Compute the Kronecker symbol for integers k from [1..K] with respect to d

    Input:  d (int): Discriminant of Dirichlet character
            K (int): Upper bound of k
    
    Output: Array of shape (K + 1, ), where chi_arr[k] is the kronecker symbol of k
            and chi_arr[0] is unused (set to 0)
    """
    chi_arr = np.zeros(K + 1, dtype=np.int8)
    for k in range(1, K + 1):
        chi_arr[k] = kronecker_symbol(d, k)
    return chi_arr

def write_kronecker(d: int, K: int, chi_arr: np.ndarray, data_dir = Path) -> None:
    """
    Save the Kronecker symbol array in both .npy and .txt formats

    Input:  d (int): Discriminant of Dirichlet character
            K (int): Upper bound for k
            chi_arr (np.ndarray): Array containing Kronecker symbol
            data_dir (Path): Path to data directory

    Output: None
    """  
    data_dir.mkdir(parents=True, exist_ok=True)

    # Determine cache folder and ensure it exists
    cache_dir = data_dir / "cache" / ("positive_d" if d > 0 else "negative_d") / f"d_{d}"
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Save binary .npy file
    npy_path = cache_dir / f"chi_{d}_K{K}.npy"
    np.save(npy_path, chi_arr)

    # Save readable text file
    txt_dir = data_dir / ("positive_d" if d > 0 else "negative_d") / f"d_{d}"
    txt_dir.mkdir(parents=True, exist_ok=True)
    txt_path = txt_dir / "kronecker.txt"
    with open(txt_path, "w") as f:
        for k in range(1, K + 1):
            f.write(f"{k} {chi_arr[k]}\n")  # only write from k = 1 to K 

def main():
    """
    Command-line interface for computing and caching Kronecker symbols
    """
    parser = argparse.ArgumentParser(description="Compute Kronecker symbols for integers [1..K] with respect to discriminant d")
    parser.add_argument("-d", "--discriminant", type=int, required=True, help='Discriminant for the Kronecker symbol')
    parser.add_argument('-K', '--upper_limit', type=int, required=True, help='Upper limit for k in the Kronecker symbol')
    parser.add_argument("-data", "--data-dir", type=str, default="data", help="Output directory to save Kronecker symbols (default: data)")
    args = parser.parse_args()

    # Retrieve command line parameters
    d = args.discriminant
    K = args.upper_limit
    data_dir = Path(args.data_dir)

    # Skip if binary cache already exists
    cache_npy = data_dir / "cache" / ("positive_d" if d > 0 else "negative_d") / f"d_{d}" / f"chi_{d}_K{K}.npy"
    if cache_npy.exists():
        return

    # Compute and write Kronecker symbols
    chi_arr = compute_kronecker(d, K)
    write_kronecker(d, K, chi_arr, data_dir)

if __name__ == "__main__":
    main()
