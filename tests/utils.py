"""
Utility functions for testing plotnado.
"""

import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig
from pybedtools import BedTool


def create_test_bigwig(chromosome="chr1", start=0, end=1000, num_points=10, values=None, noise=0.1):
    """
    Create a test bigwig file with random data.
    
    Args:
        chromosome (str): Chromosome name
        start (int): Start position
        end (int): End position
        num_points (int): Number of data points
        values (list): Optional list of values to use (overrides num_points)
        noise (float): Amount of noise to add
        
    Returns:
        str: Path to the temporary bigwig file
    """
    # Create a temporary bigwig file
    temp_bw = tempfile.NamedTemporaryFile(suffix='.bw', delete=False)
    temp_path = temp_bw.name
    temp_bw.close()
    
    # Create the bigwig file
    bw = pyBigWig.open(temp_path, "w")
    bw.addHeader([(chromosome, end)], maxZooms=0)
    
    # Create the data
    intervals = end - start
    step = intervals // num_points
    
    if values is None:
        # Create some synthetic data with a peak in the middle
        x = np.linspace(0, 1, num_points)
        # Gaussian peak
        y = 0.5 * np.exp(-((x - 0.5) ** 2) / 0.05) + noise * np.random.randn(num_points)
        y = np.maximum(0, y)  # Ensure positive values
    else:
        y = values
        num_points = len(values)
    
    # Add entries
    starts = [start + i * step for i in range(num_points)]
    ends = [start + (i + 1) * step for i in range(num_points)]
    
    bw.addEntries(chromosome, starts, ends=ends, values=y.tolist())
    bw.close()
    
    return temp_path


def create_test_bed(chromosome="chr1", regions=None, strand=None, name_prefix="feature"):
    """
    Create a test BED file.
    
    Args:
        chromosome (str): Chromosome name
        regions (list): List of (start, end) tuples
        strand (str): Strand for all features, or None for random
        name_prefix (str): Prefix for feature names
        
    Returns:
        str: Path to the temporary bed file
    """
    # Create a temporary bed file
    temp_bed = tempfile.NamedTemporaryFile(suffix='.bed', delete=False)
    temp_path = temp_bed.name
    temp_bed.close()
    
    # Default regions if none provided
    if regions is None:
        regions = [(100, 200), (300, 400), (500, 600)]
    
    # Create the bed content
    lines = []
    for i, (start, end) in enumerate(regions):
        if strand is None:
            s = "+" if np.random.rand() > 0.5 else "-"
        else:
            s = strand
        name = f"{name_prefix}{i+1}"
        score = int(10 * np.random.rand())
        lines.append(f"{chromosome}\t{start}\t{end}\t{name}\t{score}\t{s}\n")
    
    # Write to file
    with open(temp_path, 'w') as f:
        f.writelines(lines)
    
    return temp_path


def create_test_gtf(chromosome="chr1", genes=None):
    """
    Create a test GTF file with genes and exons.
    
    Args:
        chromosome (str): Chromosome name
        genes (list): List of (gene_id, start, end, strand) tuples
        
    Returns:
        str: Path to the temporary gtf file
    """
    # Create a temporary gtf file
    temp_gtf = tempfile.NamedTemporaryFile(suffix='.gtf', delete=False)
    temp_path = temp_gtf.name
    temp_gtf.close()
    
    # Default genes if none provided
    if genes is None:
        genes = [
            ("gene1", 100, 500, "+"),
            ("gene2", 700, 1200, "-")
        ]
    
    # Create the gtf content
    lines = []
    for gene_id, start, end, strand in genes:
        # Add gene entry
        gene_attr = f'gene_id "{gene_id}"; gene_name "{gene_id}";'
        lines.append(f"{chromosome}\t.\tgene\t{start}\t{end}\t.\t{strand}\t.\t{gene_attr}\n")
        
        # Add transcript entry
        transcript_id = f"{gene_id}_transcript1"
        transcript_attr = f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'
        lines.append(f"{chromosome}\t.\ttranscript\t{start}\t{end}\t.\t{strand}\t.\t{transcript_attr}\n")
        
        # Add exons - create 2-3 exons per gene
        exon_count = 2 if start + 200 <= end else 1
        
        if exon_count == 1:
            # Single exon gene
            exon_attr = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; exon_number "1";'
            lines.append(f"{chromosome}\t.\texon\t{start}\t{end}\t.\t{strand}\t.\t{exon_attr}\n")
        else:
            # Multi-exon gene
            exon_length = (end - start) // (2 * exon_count - 1)  # Leave space for introns
            for i in range(exon_count):
                exon_start = start + i * 2 * exon_length
                exon_end = exon_start + exon_length
                exon_attr = f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; exon_number "{i+1}";'
                lines.append(f"{chromosome}\t.\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\t{exon_attr}\n")
    
    # Write to file
    with open(temp_path, 'w') as f:
        f.writelines(lines)
    
    return temp_path


def tabix_index_file(file_path):
    """
    Create a tabix index for a file.
    
    Args:
        file_path (str): Path to the file to index
        
    Returns:
        str: Path to the indexed file
    """
    # Check file type and determine correct parameters
    if file_path.endswith('.bed'):
        bt = BedTool(file_path)
        return bt.tabix(force=True).fn
    elif file_path.endswith('.gtf'):
        from plotnado.tracks import tabix_gtf
        return tabix_gtf(Path(file_path))
    else:
        raise ValueError(f"Unsupported file type: {file_path}")


def create_test_genomic_region(chromosome="chr1", start=1000, end=2000, strand="+"):
    """
    Create a test GenomicRegion instance.
    
    Args:
        chromosome (str): Chromosome name
        start (int): Start position
        end (int): End position
        strand (str): Strand
        
    Returns:
        GenomicRegion: A GenomicRegion instance
    """
    from plotnado.tracks import GenomicRegion
    return GenomicRegion(
        chromosome=chromosome,
        start=start,
        end=end,
        strand=strand
    )
