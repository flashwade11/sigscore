from typing import Dict, List

import hdf5plugin  # noqa: F401
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import issparse


def filter_signatures(
    signatures: Dict[str, List[str]],
    reference_genes: List[str],
    conserved: float = 0.7
) -> Dict[str, List[str]] | None:
    all_genes = []
    for gene_list in signatures.values():
        all_genes.extend(gene_list)
    
    if all(gene in reference_genes for gene in all_genes):
        return signatures
    
    num_original_sigs = len(signatures)
    num_original_genes_per_sig = {name: len(genes) for name, genes in signatures.items()}
    
    filtered_sigs = {
        name: [gene for gene in genes if gene in reference_genes] \
            for name, genes in signatures.items()
    }
    
    num_filtered_genes_per_sig = {name: len(genes) for name, genes in filtered_sigs.items()}
    
    if all(count == 0 for count in num_filtered_genes_per_sig.values()):
        return None
    
    frac_conserved = {
        name: num_filtered_genes_per_sig[name] / num_original_genes_per_sig[name] \
            for name in filtered_sigs.keys()
    }
    
    filtered_sigs = {
        name: genes for name, genes in filtered_sigs.items() \
            if frac_conserved[name] >= conserved
    }
    
    num_filtered_sigs = len(filtered_sigs)
    
    if num_filtered_sigs == 0:
        raise ValueError("No signatures left to score with.")
    
    if num_filtered_sigs < num_original_sigs:
        print(f"Filtered {num_original_sigs - num_filtered_sigs} out of {num_original_sigs} signatures with < {conserved:.2%} conserved genes.")
    
    return filtered_sigs


def compute_base_signature_scores(
    X: np.ndarray,
    signatures: Dict[str, List[str]],
    gene_names: List[str],
    conserved: float = 0.7,
) -> pd.DataFrame:
    filtered_sigs = filter_signatures(
        signatures=signatures,
        reference_genes=gene_names,
        conserved=conserved
    )
    
    scores_dict = {}
    
    for sig_name, genes in filtered_sigs.items():
        gene_indices = [gene_names.index(gene) for gene in genes if gene in gene_names]
        if gene_indices:
            scores_dict[sig_name] = np.mean(X[:, gene_indices], axis=1)
    
    return pd.DataFrame(scores_dict)


def sample_bin_matched_genes(
    X: np.ndarray,
    gene_names: List[str],
    n_bins: int = 30
) -> Dict[str, int]:
    gene_means = np.mean(X, axis=0)
    bins = pd.qcut(gene_means, q=n_bins, labels=False, duplicates="drop")
    return {gene_names[i]: bin_ for i, bin_ in enumerate(bins)}

def bin_match(
    signature_genes: List[str],
    bins: Dict[str, int],
    n: int = 100,
    replace: bool = False
) -> List[str]:
    all_genes = list(bins.keys())
    matched_genes = []
    for gene in signature_genes:
        if gene in bins:
            gene_bin = bins[gene]
            bin_genes = [g for g in all_genes if g != gene and bins[g] == gene_bin]
            if len(bin_genes) > 0:
                sample_size = min(len(bin_genes), n)
                sampled = np.random.choice(bin_genes, size=sample_size, replace=replace)
                matched_genes.extend(sampled)
    return matched_genes

def compute_signature_scores(
    adata: AnnData,
    signatures: Dict[str, List[str]],
    conserved: float = 0.7,
    sample_key: str | None = None,
    layer: str | None = None,
    center_genes: bool = True,
    center: bool = True,
    expr_center: bool = True,
    expr_bin_adata: AnnData | None = None,
    expr_bins: Dict[str, int] | None = None,
    expr_sigs: Dict[str, List[str]] | None = None,
    expr_nbin: int = 30,
    expr_binsize: int = 100,
    replace: bool = False,
    key_added: str = "sig_scores",
) -> pd.DataFrame:
    filtered_sigs = filter_signatures(signatures=signatures, reference_genes=adata.var_names, conserved=conserved)
    if sample_key is not None:
        all_scores = []
        all_samples = adata.obs[sample_key].unique()
        for sample in all_samples:
            sample_adata = adata[adata.obs[sample_key] == sample].copy()
            sample_scores = compute_signature_scores_core(
                adata=sample_adata,
                signatures=filtered_sigs,
                layer=layer,
                center_genes=center_genes,
                center=center,
                expr_center=expr_center,
                expr_bin_adata=expr_bin_adata,
                expr_bins=expr_bins,
                expr_sigs=expr_sigs,
                expr_nbin=expr_nbin,
                expr_binsize=expr_binsize,
                conserved=conserved,
                replace=replace
            )
            all_scores.append(sample_scores)
        scores_df = pd.concat(all_scores)
        scores_df = scores_df.loc[adata.obs_names]
    else:
        scores_df = compute_signature_scores_core(
            adata=adata,
            signatures=filtered_sigs,
            layer=layer,
            center_genes=center_genes,
            center=center,
            expr_center=expr_center,
            expr_bin_adata=expr_bin_adata,
            expr_bins=expr_bins,
            expr_sigs=expr_sigs,
            expr_nbin=expr_nbin,
            expr_binsize=expr_binsize,
            conserved=conserved,
            replace=replace
        )
    adata.obsm[key_added] = scores_df
    return scores_df
              
def compute_signature_scores_core(
    adata: AnnData,
    signatures: Dict[str, List[str]],
    layer: str | None = None,
    conserved: float = 0.7,
    center_genes: bool = True,
    center: bool = True,
    expr_center: bool = True,
    expr_bin_adata: AnnData | None = None,
    expr_bins: Dict[str, int] | None = None,
    expr_sigs: Dict[str, List[str]] | None = None,
    expr_nbin: int = 30,
    expr_binsize: int = 100,
    replace: bool = False
) -> pd.DataFrame:
    X = adata.layers[layer] if layer is not None else adata.X
    if issparse(X):
        X = X.toarray()
    
    X_centered = X - np.mean(X, axis=0, keepdims=True) if center_genes else X
        
    gene_names = adata.var_names.tolist()
    scores = compute_base_signature_scores(
        X=X_centered,
        signatures=signatures,
        gene_names=gene_names,
        conserved=conserved
    )
    
    if not center:
        expr_center = False
    
    if expr_center:
        if expr_sigs is None:
            if expr_bins is None:
                if expr_bin_adata is None:
                    expr_bin_X = X
                else:
                    expr_bin_X = expr_bin_adata.layers[layer] if layer is not None else expr_bin_adata.X
                    if issparse(expr_bin_X):
                        expr_bin_X = expr_bin_X.toarray()
                expr_bins = sample_bin_matched_genes(X=expr_bin_X, gene_names=gene_names, n_bins=expr_nbin)
                all_sig_genes = set()
                for gene_list in signatures.values():
                    all_sig_genes.update(gene_list)
                if not all(gene in expr_bins for gene in all_sig_genes):
                    raise ValueError("Some signature genes are not found in expression bins")
            expr_sigs = {
                name: bin_match(genes, bins=expr_bins, n=expr_binsize, replace=replace) \
                    for name, genes in signatures.items()
            }
        
        expr_scores = compute_base_signature_scores(
            X=X_centered if center_genes else X,
            signatures=expr_sigs,
            gene_names=gene_names
        )
        
        scores = scores.subtract(expr_scores)
    elif center:
        if expr_bin_adata is not None:
            expr_bin_X = expr_bin_adata.layers[layer] if layer is not None else expr_bin_adata.X
            if issparse(expr_bin_X):
                expr_bin_X = expr_bin_X.toarray()
            
            center_scores = np.mean(expr_bin_X, axis=1)
        else:
            center_scores = np.mean(X, axis=1)
        for col in scores.columns:
            scores[col] = scores[col] - center_scores
            
    scores.index = adata.obs_names
    return scores