# sigscore

A fast, lightweight, and parallelizable Python library for gene signature scoring in single-cell RNA-seq data.

## Installation

```bash
pip install sigscore
```

Or install from source:

```bash
git clone https://github.com/username/sigscore.git
cd sigscore
pip install -e .
```

## Quick Start

```python
import scanpy as sc
import sigscore as ss

# Load your AnnData object
adata = sc.read_h5ad("your_data.h5ad")

# Define your gene signatures
signatures = {
    "Signature1": ["GENE1", "GENE2", "GENE3"],
    "Signature2": ["GENE4", "GENE5", "GENE6"]
}

# Compute signature scores
scores = ss.compute_signature_scores(
    adata=adata,
    signatures=signatures,
    conserved=0.7,  # Minimum fraction of genes to be conserved
    key_added="sig_scores"  # Where to store the results in adata.obsm
)

# Access the scores
print(scores.head())
# Or access from the AnnData object
print(adata.obsm["sig_scores"].head())
```

## Using Built-in Signatures

sigscore comes with built-in gene signatures that you can use directly:

```python
# List available signatures
print(ss.list_signatures())

# Use tumor meta-programs
from sigscore import tumor_MPs

scores = ss.compute_signature_scores(
    adata=adata,
    signatures=tumor_MPs,
    key_added="tumor_mp_scores"
)
```

## Advanced Usage

### Cell Annotation Based on Signatures

You can annotate cells based on their highest-scoring signatures:

```python
# Annotate cells across different samples
annotations = ss.annotate_signatures_per_sample(
    adata=adata,
    signatures=signatures,
    sample_key="sample",  # Column in adata.obs that identifies samples
    added_key="signature_annotation",  # Where to store annotations in adata.obs
    min_conserved_genes=0.4,  # Minimum fraction of genes to be conserved
    min_MP_score=0.0,  # Minimum score threshold
    min_cell_fraction=0.05  # Minimum fraction of cells per signature
)

# Access annotations
print(adata.obs["signature_annotation"].value_counts())
```

### Parallel Processing

For large datasets, you can use parallel processing to speed up the annotation:

```python
annotations = ss.annotate_signatures_per_sample(
    adata=adata,
    signatures=signatures,
    sample_key="sample",
    parallel=True,  # Enable parallel processing
    n_cpus=4,  # Number of CPU cores to use
    verbose=10  # Verbosity level
)
```

## API Reference

### Core Functions

- `compute_signature_scores(adata, signatures, ...)`: Compute signature scores for cells in an AnnData object
- `annotate_signatures_per_sample(adata, signatures, sample_key, ...)`: Annotate cells based on signature scores across samples
- `annotate_signatures_to_cells(adata, signatures, ...)`: Annotate cells based on signature scores in a single sample

### Parameters

- `conserved`: Minimum fraction of genes to be conserved in the signature (default: 0.7)
- `center_genes`: Whether to center gene expression (default: True)
- `center`: Whether to center signature scores (default: True)
- `expr_center`: Whether to center expression-based scores (default: True)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the terms of the MIT license.
