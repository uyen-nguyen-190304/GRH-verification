# Generalized Riemann Hypothesis Verification

A lightweight numerical tool-chain for verifying the Generalized Riemann Hypothesis for quadratic Dirichlet L-functions


## 1. Building the C++ extension

```bash
mkdir build ** cd build
cmake ..                           # need CMake >= 3.16
make                               # build grhverify*.so
cd ..
```

## 2. Telling the pipeline where **lcalc** lives

```jsonc
// local_config.json   (repo root)
{
  "lcalc_path": "/usr/bin/lcalc"   // edit if lcalc lives elsewhere
}
```

## 3. Running the verification pipeline

The entry point to the verification pipeline is `script/pipeline.py`. To run the verification pipeline, you have to provide either:

* exactly one discriminant (`-d/--discriminant`), or
* an inclusive range (`--d-min` and `--d-max`)

The two options are mutually exclusive. If you are trying to pass both flag/flags combination, the script will abort with:

```bash
Error: Error: Cannot provide both --discriminant and --d-min/d-max
```

### A. Single discriminant

```bash
python -m script.pipeline -d 1299721
```

### B. Custom range of discriminants

```bash
python -m script.pipeline --d-min -1000 --d-max 1000
```

One of the two options above must be provided. Then, for each fundamental $d$ in the set, the pipeline keeps adding zeros until the RH is verified and records the smallest number of zeros $N$ needed. Besides, there are some other optional flags that could be supplied to the pipeline if wanted:

| Flag                       | Default                 | Meaning                                                                           |
| -------------------------- | ----------------------- | --------------------------------------------------------------------------------- |
| `-eta`                     | *first zero*            | Verification height $η$.  If omitted, the ordinate of the first zero is used |
| `-K`, `--upper-limit`      | **$100000 = 10^5$**              | Length of the Λ(k) and χ<sub>d</sub>(k) arrays                                  |
| `-eps`, `--epsilon`        | **$1 × 10^{-6}$**            | Half-width $ε$ of each symmetric interval $[γ-ε, γ+ε]$                         |
| `-config`, `--config-file` | **`local_config.json`** | JSON containing `"lcalc_path": "/path/to/lcalc"`                                 |
| `-data`, `--data-dir`      | **`data`**              | Root directory for caches (Λ, χ, zeros)                                          |
| `-output`, `--output-dir`  | **`results`**           | Destination for `summary.csv` & `errors.log`                                     |

---

### Complete example

```bash
python -m script.pipeline --d-min -1000 --d-max 1000    \
        -K 200000 --epsilon 5e-6 --eta 0.35             \
        --config local_config.json                      \
        --data-dir data --output-dir results
```
 
## 4. Cleaning Up

Remove build artefacts, Python byte-code, and generated caches:

```bash
# delete CMake build dir and compiled .so
rm -rf build grhverify/*.so

# delete __pycache__ folders
find . -name '__pycache__' -type d -exec rm -r {} +

# delete cached zeros / Λ / χ and previous results
rm -rf data/ results/
```
