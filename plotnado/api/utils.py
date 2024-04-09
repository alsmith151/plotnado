import pathlib
from typing import Iterable, Optional, Union

import coolbox.api as cb
import pyBigWig
import pyranges as pr


def get_max_value_over_region(
    bigwigs: Iterable[Union[pathlib.Path, str]],
    region: cb.GenomeRange,
    exclude: cb.GenomeRange = None,
) -> float:
    """
    Get the maximum value over a region in a list of bigwig files

    Args:
        bigwigs (Iterable[Union[pathlib.Path, str]]): List of bigwig files
        region (cb.GenomeRange): Region to get the maximum value over
        exclude (cb.GenomeRange, optional): Region to exclude from the maximum value calculation. Defaults to None.
    """

    max_values = []
    for bw in bigwigs:
        with pyBigWig.open(str(bw)) as f:
            chrom_len = f.chroms(region.chrom)
            start = int(max(0, region.start))
            end = int(min(region.end, chrom_len))

            if exclude:
                # Split the region into two parts
                if start < exclude.start:
                    start1, end1 = start, exclude.start
                    start2, end2 = exclude.end, end
                else:
                    start1, end1 = exclude.end, end
                    start2, end2 = start, exclude.start

                max_val1 = f.stats(region.chrom, start1, end1, type="max")[0]
                max_val2 = f.stats(region.chrom, start2, end2, type="max")[0]
                max_values.append(max(max_val1, max_val2))

            else:
                max_val = f.stats(region.chrom, start, end, type="max")[0]
                max_values.append(max_val)

    return round(max(max_values))


def interval_to_pyranges(interval: cb.GenomeRange) -> pr.PyRanges:
    """
    Convert a GenomeRange to a PyRanges object

    Args:
        interval (cb.GenomeRange): GenomeRange to convert
    """

    return pr.PyRanges(
        {
            "Chromosome": [interval.chrom],
            "Start": [interval.start],
            "End": [interval.end],
        }
    )


def get_human_readable_number_of_bp(bp: int) -> str:
    """Converts integer into human readable basepair number"""

    if bp < 1000:
        bp = f"{bp}bp"
    elif (bp / 1e3) < 1000:
        bp = f"{bp / 1e3}kb"
    elif (bp / 1e6) < 1000:
        bp = f"{bp / 1e6}mb"

    return bp
