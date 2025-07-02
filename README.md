# Generalized Riemann Hypothesis Verification for Quadratic Dirichlet L-Functions

A Python toolkit for verifying instances of the Generalized Riemann Hypothesis (GRH) for quadratic Dirichlet L-functions. It calculates Kronecker symbols and von Mangoldt values, fetches zeros using `lcalc`, and verifies the absence of zeros up to a specified height $\eta$.


## Repository Structure

```
GRH-verification/
├─ grhverify/
│   ├─ __init__.py
│   ├─ base_case.py
│   ├─ utils/
│   │   ├─ discriminants.py
│   │   ├─ generate_zeros.py
│   │   ├─ kronecker_symbol.py
│   │   └─ von_mangoldt.py
│   └─ … (planned future modules)
│
├─ driver.py               # Command-line interface entry point
├─ local_config.json       # e.g., {"lcalc_path": "/path/to/lcalc"}
├─ README.md               # Project description file
├─ pyproject.toml          # Package metadata for editable installs
└─ tests/
    └─ (unit tests in development)
```


## Installation

```bash
# Clone the repository
git clone https://github.com/uyen-nguyen-190304/GRH-verification.git

# Navigate into the repository directory
cd GRH-verification

# Install the package in editable mode
pip install -e .                                  

# Activate SageMath environment
conda activate sage                                  
```


## Configuration

Create a `local_config.json` in the project root with similar structure to the provided `example_config.json`:

```json
{
  "lcalc_path": "/full/path/to/lcalc"
}
```

This specifies the location of the `lcalc` executable.


## Command-Line Parameters

Run `python driver.py --help` for detailed usage. Key parameters include:

| Option                     | Arg     | Default             | Description                                                         |
| -------------------------- | ------- | ------------------- | ------------------------------------------------------------------- |
| `-d`, `--discriminant`     | *int*   | —                   | Single discriminant to verify. Exclusive with `--d-min/max`        |
| `--d-min`                  | *int*   | —                   | Minimum discriminant (inclusive)                                   |
| `--d-max`                  | *int*   | —                   | Maximum discriminant (inclusive)                                   |
| `-eta`, `--height`         | *float* | auto                | Height η; defaults to (first positive zero + 2 $\varepsilon$) if unspecified    |
| `-k`, `--power`            | *int*   | `1`                 | Order of logarithmic derivative $L^{(k)}/L$ (currently only `1`) |
| `-K`, `--upper-limit`      | *int*   | `100000`            | Truncation limit $K$ for $\chi$ and $\Lambda$ arrays                          |
| `-eps`, `--epsilon`        | *float* | `1e-6`              | Half-width $\varepsilon$ for zero intervals $[\gamma - \varepsilon, \gamma + \varepsilon]$                    |
| `-config`, `--config-file` | *path*  | `local_config.json` | JSON config file path                                              |
| `-data`, `--data-dir`      | *path*  | `data`              | Directory for cached zeros, intervals, and $\chi, \Lambda$ data                |
| `-output`, `--output-dir`  | *path*  | `results`           | Output directory for results and error logs                              |


## Quickstart Examples

### Verify a range of discriminants

```bash
python driver.py --d-min -500 --d-max 500 \
  --eps 1e-6 --upper-limit 50000
```

### Verify a single discriminant with explicit height

```bash
python driver.py -d 13 -eta 3.5
```


## Output Files

Output directories:

* **`results/summary.csv`**
  * CSV summary: `d, eta, N_needed`

* **`results/errors.log`**
  * Logs runtime errors or failures

* **`data/von_mangoldt.txt`**
  * von Mangoldt values $\Lambda$

* **`data/positive_d/d_<d>/`**, **`data/negative_d/d_<d>/`**
  * `zeros.txt`: zeros $\gamma$
  * `intervals.txt`: intervals $[\gamma - \varepsilon, \gamma + \varepsilon]$
  * `kronecker.txt`: Kronecker symbols $\chi$


## Cleaning Up

Remove generated files and cached data:

```bash
# Remove results and logs
rm -rf results/

# Clear cached data
rm -rf data/

# Remove Python caches
find . \
  -type d \( -name "__pycache__" -o -name ".pytest_cache" -o -name ".mypy_cache" \) \
  -prune \
  -exec rm -rf {} +

# Remove build artifacts
rm -rf build/ dist/ *.egg-info

# Uninstall if previously installed in editable mode
pip uninstall grhverify
```


## Future Development

* Implement verification for higher-order logarithmic derivatives
* Introduce unit testing for arithmetic functions
