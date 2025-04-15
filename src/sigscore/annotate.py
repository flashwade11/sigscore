from typing import Dict, List

import pandas as pd
from anndata import AnnData
from joblib.parallel import Parallel, delayed
from tqdm.rich import tqdm

from .core import compute_signature_scores


def annotate_signatures_to_cells(
    adata: AnnData,
    signatures: Dict[str, List[str]],
    min_conserved_genes: float = 0.4,
    min_MP_score: float = 0.0,
    min_cell_fraction: float = 0.05
):
    scores_df: pd.DataFrame = compute_signature_scores(
        adata=adata,
        signatures=signatures,
        conserved=min_conserved_genes,
        key_added="mp_scores"
    )
    
    max_scores = scores_df.max(axis=1)
    cells_to_keep = max_scores >= min_MP_score
    if not all(cells_to_keep):
        scores_df = scores_df.loc[cells_to_keep]
    total_cells = len(scores_df)
    if total_cells == 0:
        return None
    MP_for_cells = scores_df.idxmax(axis=1)
    MP_counts = MP_for_cells.value_counts()
    MP_to_keep = (MP_counts / total_cells) >= min_cell_fraction
    if not any(MP_to_keep):
        return None
    return MP_for_cells[MP_for_cells.map(MP_to_keep)]

def annotate_signatures_per_sample(
    adata: AnnData,
    signatures: Dict[str, List[str]],
    sample_key: str,
    added_key: str = "sigscore_annotation",
    parallel: bool = False,
    n_cpus: int = 1,
    min_conserved_genes: float = 0.4,
    min_MP_score: float = 0.0,
    min_cell_fraction: float = 0.05,
    verbose: int = 10,
) -> pd.DataFrame:
    sample_adata_list = [
        adata[adata.obs[sample_key] == sample].copy() \
            for sample in adata.obs[sample_key].unique()
    ]
    
    if parallel:
        if n_cpus < 1 and not isinstance(n_cpus, int):
            raise ValueError("n_cpus must be a positive integer")
        results = Parallel(n_jobs=n_cpus, verbose=verbose)(
            delayed(annotate_signatures_to_cells)(
                adata=sample_adata,
                signatures=signatures,
                min_conserved_genes=min_conserved_genes,
                min_MP_score=min_MP_score,
                min_cell_fraction=min_cell_fraction
            ) for sample_adata in sample_adata_list
        )
    else:
        results = []
        for sample_adata in tqdm(sample_adata_list):
            results.append(annotate_signatures_to_cells(
                adata=sample_adata,
                signatures=signatures,
                min_conserved_genes=min_conserved_genes,
                min_MP_score=min_MP_score,
                min_cell_fraction=min_cell_fraction
            ))
    results = pd.concat([result for result in results if result is not None])
    adata.obs[added_key] = results
    return results