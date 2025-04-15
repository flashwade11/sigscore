from .core import compute_signature_scores
from .meta_program import assgin_MP_to_samples, assign_MP_to_cells
from .signatures import tumor_MPs

__all__ = ["compute_signature_scores", "assgin_MP_to_samples", "assign_MP_to_cells", "tumor_MPs"]
__version__ = "0.1.0"


def list_signatures():
    return ["tumor_MPs"]
