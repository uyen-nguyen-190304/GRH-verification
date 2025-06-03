# GRH Verification

## Usage

1. Compile the C++ executable

```
 g++ -o src/rh_verify src/rh_verify.cpp -std=c++11
```

2. Run the pipeline

```
python3 pipeline.py -d <discriminant> [options]
```

Required:
* `-d`: discriminant (integer)

Optional:
* `-N`: number of zeros to compute (default: 100)
* `-K`: upper limit for von Mangoldt and Kronecker computation (default: 20)
* `-eta`: height parameter (default: computed from zeros or 50.0 if fail)
* `-err`: error window for intervals (default: 1e-8)
* `-config`: path to config.json file (default: `example_config.json`)

Example

```
python3 pipeline.py -d 5
```

## Output Files

* `data/zeros.txt`: computed zeros of the Dirichlet L-function
* `data/intervals.txt`: intervals around zeros
* `data/von_mangoldt.txt`: von Mangoldt values
* `data/kronecker.txt`: Kronecker symbols
* Results are currently printed to the console