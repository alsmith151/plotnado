"""
Aesthetics factory for compatibility.
"""

from pydantic import BaseModel
from .bigwig import BigwigAesthetics
from .scalebar import ScaleBarAesthetics
from .genes import GenesAesthetics


class Aesthetics(BaseModel):
    """
    Factory class for track aesthetics.
    Maintains compatibility with the old API.
    """

    @classmethod
    def bigwig(cls, **kwargs) -> BigwigAesthetics:
        """Create BigwigAesthetics instance."""
        return BigwigAesthetics(**kwargs)

    @classmethod
    def scalebar(cls, **kwargs) -> ScaleBarAesthetics:
        """Create ScaleBarAesthetics instance."""
        return ScaleBarAesthetics(**kwargs)

    @classmethod
    def genes(cls, **kwargs) -> GenesAesthetics:
        """Create GenesAesthetics instance."""
        return GenesAesthetics(**kwargs)
