from .core import compute_signature_scores
from .annotate import annotate_signatures_per_sample, annotate_signatures_to_cells
from .signatures import tumor_MPs

__all__ = ["compute_signature_scores", "annotate_signatures_per_sample", "annotate_signatures_to_cells", "tumor_MPs"]
__version__ = "0.1.0"


def list_signatures():
    return ["tumor_MPs"]
