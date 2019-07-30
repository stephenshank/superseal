# convex-qsr

Reference-guided viral quasispecies reconstruction using convex optimization. **Under development**.

## Install

Available via [conda](https://anaconda.org/) from the [stephenshank](https://anaconda.org/stephenshank/) channel (coming soon to [bioconda](https://bioconda.github.io)):

```
conda config --add channels stephenshank
conda install convex-qsr
```

Also available via `pip` and the [PyPI](https://pypi.org/project/convex-qsr/):

```
pip install convex-qsr
```

## Run

```
cqsr_ec -b sorted.bam -o output
cqsr_rg -i output -o output
cqsr_qr -i output -o output
```

## Test
```
python -m unittest discover
```