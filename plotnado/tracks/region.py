"""
Genomic region model.
"""

from typing import Any
from pydantic import BaseModel, ConfigDict, Field

from .enums import Strand


class GenomicRegion(BaseModel):
    """
    Represents a genomic region with chromosome, start, end, and strand.

    Attributes:
        chromosome: Chromosome name (e.g., 'chr1')
        start: Start position (0-based)
        end: End position
        strand: Strand ('+' or '-')
    """

    chromosome: str = Field(description="Chromosome name (for example `chr1`).")
    start: int = Field(description="0-based inclusive start coordinate.")
    end: int = Field(description="0-based exclusive end coordinate.")
    strand: Strand = Field(default=Strand.PLUS, description="Strand orientation for directional operations.")

    model_config = ConfigDict(use_enum_values=True)

    @property
    def length(self) -> int:
        """Length of the region in base pairs."""
        return self.end - self.start

    @property
    def center(self) -> int:
        """Center position of the region."""
        return (self.start + self.end) // 2

    def extend(self, upstream: int = 0, downstream: int = 0) -> "GenomicRegion":
        """
        Extend the region upstream and/or downstream.

        Args:
            upstream: Number of base pairs to extend upstream
            downstream: Number of base pairs to extend downstream

        Returns:
            A new GenomicRegion instance
        """
        new_start = self.start
        new_end = self.end

        if self.strand == "+":
            new_start -= upstream
            new_end += downstream
        else:  # strand == "-"
            new_start -= downstream
            new_end += upstream

        # Ensure start is not negative
        new_start = max(0, new_start)

        return GenomicRegion(
            chromosome=self.chromosome,
            start=new_start,
            end=new_end,
            strand=self.strand,
        )

    def __str__(self) -> str:
        return f"{self.chromosome}:{self.start}-{self.end}({self.strand})"

    @classmethod
    def from_str(cls, region_str: str) -> "GenomicRegion":
        """
        Parse a region string like 'chr1:1000-2000' or 'chr1:1M-2M' or 'chr1:1M+20k'.

        Args:
            region_str: Region string to parse

        Returns:
            GenomicRegion instance
        """
        from .utils import parse_genomic_value

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

        if "-" in positions:
            start_str, end_str = positions.split("-")
            start = parse_genomic_value(start_str)
            end = parse_genomic_value(end_str)
        elif "+" in positions:
            center_str, size_str = positions.split("+")
            center = parse_genomic_value(center_str)
            size = parse_genomic_value(size_str)
            start = center - (size // 2)
            end = center + (size // 2)
        else:
            raise ValueError(
                f"Invalid positions: {positions}. Expected 'start-end' or 'center+size'"
            )

        return cls(
            chromosome=chromosome,
            start=max(0, start),
            end=end,
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


GenomicRange = GenomicRegion
