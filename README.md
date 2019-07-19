# convex-qsr

Reference-guided viral quasispecies reconstruction using convex optimization. **Under development**.

## Install

```
git clone https://github.com/stephenshank/convex-qsr
cd convex-qsr
python setup.py install
```

`pip` and `conda` installable packages coming soon.

## Run

```
cqsr_ec -b sorted.bam -o output
cqsr_rg -i output -o output
cqsr_qr -i output -o output
cqsr_viz -o output
```

## Test
```
python -m unittest discover
```