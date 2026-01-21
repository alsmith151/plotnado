"""
Genomic region model.
"""

from typing import Any, Literal
from pydantic import BaseModel


class GenomicRegion(BaseModel):
    """
    Represents a genomic region with chromosome, start, end, and strand.

    Attributes:
        chromosome: Chromosome name (e.g., 'chr1')
        start: Start position (0-based)
        end: End position
        strand: Strand ('+' or '-')
    """

    chromosome: str
    start: int
    end: int
    strand: Literal["+", "-"] = "+"

    @property
    def length(self) -> int:
        """Length of the region in base pairs."""
        return self.end - self.start

    @property
    def center(self) -> int:
        """Center position of the region."""
        return (self.start + self.end) // 2

    def __str__(self) -> str:
        return f"{self.chromosome}:{self.start}-{self.end}({self.strand})"

    @classmethod
    def from_str(cls, region_str: str) -> "GenomicRegion":
        """
        Parse a region string like 'chr1:1000-2000' or 'chr1:1000-2000(+)'.

        Args:
            region_str: Region string to parse

        Returns:
            GenomicRegion instance
        """
        strand = "+"
        if "(" in region_str:
            main_part, strand_part = region_str.rsplit("(", 1)
            strand = strand_part.rstrip(")")
            region_str = main_part

        if ":" not in region_str:
            raise ValueError(
                f"Invalid region string: {region_str}. Expected 'chrom:start-end'"
            )

        chromosome, positions = region_str.split(":")
        if "-" not in positions:
            raise ValueError(f"Invalid positions: {positions}. Expected 'start-end'")

        start_str, end_str = positions.split("-")
        return cls(
            chromosome=chromosome,
            start=int(start_str.replace(",", "")),
            end=int(end_str.replace(",", "")),
            strand=strand,
        )

    @classmethod
    def from_tuple(cls, region_tuple: tuple) -> "GenomicRegion":
        return cls(
            chromosome=region_tuple[0],
            start=region_tuple[1],
            end=region_tuple[2],
            strand=region_tuple[3],
        )

    @classmethod
    def from_list(cls, region_list: list) -> "GenomicRegion":
        return cls(
            chromosome=region_list[0],
            start=region_list[1],
            end=region_list[2],
            strand=region_list[3],
        )

    @classmethod
    def from_dict(cls, region_dict: dict) -> "GenomicRegion":
        return cls(
            chromosome=region_dict["chromosome"],
            start=region_dict["start"],
            end=region_dict["end"],
            strand=region_dict["strand"],
        )

    @classmethod
    def from_named_tuple(cls, region_named_tuple: Any) -> "GenomicRegion":
        return cls(
            chromosome=region_named_tuple.chromosome,
            start=region_named_tuple.start,
            end=region_named_tuple.end,
            strand=region_named_tuple.strand,
        )

    @classmethod
    def into(cls, gr: Any) -> "GenomicRegion":
        """Convert input into a GenomicRegion."""
        if isinstance(gr, GenomicRegion):
            return gr
        if isinstance(gr, str):
            return cls.from_str(gr)
        if isinstance(gr, tuple):
            return cls.from_tuple(gr)
        if isinstance(gr, list):
            return cls.from_list(gr)
        if isinstance(gr, dict):
            return cls.from_dict(gr)
        if hasattr(gr, "chromosome") and hasattr(gr, "start"):
            return cls.from_named_tuple(gr)
        raise ValueError(f"Cannot convert {type(gr)} to GenomicRegion")
