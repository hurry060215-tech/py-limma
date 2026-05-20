# limmapy

A **pure-Python port of [Bioconductor limma](https://bioconductor.org/packages/release/bioc/html/limma.html)** (Ritchie et al., *Nucleic Acids Research* 2015) — linear models and empirical-Bayes moderated statistics for differential expression.

- **No `rpy2`**, no R install — the core limma linear-model workflow reimplemented in NumPy / SciPy
- The canonical pipeline: `voom → lmFit → contrasts.fit → eBayes → topTable`
- Smyth 2004 empirical-Bayes variance moderation; `treat`, `decideTests`, `duplicateCorrelation`, `removeBatchEffect`
- Both Python-style (`lm_fit`, `ebayes`, `top_table`) and R-style (`lmFit`, `eBayes`, `topTable`) names exported
- `pandas` / `numpy`-friendly

> The import name is **`pylimma`**; the PyPI distribution name is **`omicverse-limma`** (`pip install omicverse-limma`). The names `pylimma` / `limmapy` were unavailable on PyPI, so the distribution carries the omicverse-ecosystem name.

> This is a **standalone mirror** of the implementation developed in [`omicverse`](https://github.com/Starlitnightly/omicverse). It powers `omicverse`'s edgeR / limma-voom differential-expression backends and the `pydeqms` proteomics workflow.

## Install

```bash
pip install omicverse-limma
```

## Quick start

```python
import numpy as np
import pylimma

# expr: genes x samples log-expression matrix; design: samples x coefficients
fit = pylimma.lmFit(expr, design)
fit = pylimma.eBayes(fit)
res = pylimma.topTable(fit, coef=1, number=np.inf)
res.head()
```

### RNA-seq with voom

```python
import pylimma
v = pylimma.voom(counts, design)            # EList: mean-variance weights
fit = pylimma.eBayes(pylimma.lmFit(v.E, design, weights=v.weights))
pylimma.topTable(fit, coef=1)
```

## API

| Python | R counterpart |
|---|---|
| `lm_fit` / `lmFit` | `lmFit` |
| `ebayes` / `eBayes` | `eBayes` |
| `contrasts_fit` | `contrasts.fit` |
| `top_table` / `topTable` | `topTable` |
| `voom` | `voom` |
| `treat` | `treat` |
| `decide_tests` / `decideTests` | `decideTests` |
| `duplicate_correlation` / `duplicateCorrelation` | `duplicateCorrelation` |
| `remove_batch_effect` / `removeBatchEffect` | `removeBatchEffect` |
| `MArrayLM`, `EList`, `TestResults` | the corresponding S4 classes |

## Citation

> Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., Smyth, G.K. **limma powers differential expression analyses for RNA-sequencing and microarray studies.** *Nucleic Acids Research* 43(7), e47 (2015).

…and acknowledge omicverse / this repo for the Python port.

## License

LGPL-3.0-or-later — matches the upstream Bioconductor package.
