"""
Data schemas for plotnado tracks.
"""

import pandas as pd
import pandera.pandas as pa
from pandera.typing import Series


class BedgraphDataFrameSchema(pa.DataFrameModel):
    """Schema for Bedgraph-like data (chrom, start, end, value)."""

    chrom: Series[str]
    start: Series[int]
    end: Series[int]
    value: Series[float]

    class Config:
        coerce = True
        strict = True


class BedgraphDataFrame(pd.DataFrame):
    """
    A pandas DataFrame subclass that validates against BedgraphDataFrameSchema.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not self.empty:
            BedgraphDataFrameSchema.validate(self, inplace=True)
            self["chrom"] = self["chrom"].astype(object)
