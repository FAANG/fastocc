from dataclasses import dataclass, replace
from typing import Iterator

import numpy as np
import pandas as pd

from .range import Range
from .options import Options, get_options


@dataclass
class Fragments:
    """
    Dataclass containing raw ATAC fragments from a single genomic interval.
    """
    range: Range
    start: np.ndarray
    size: np.ndarray
    opts: Options = None

    def filter_by_max_size(self, max_size: int):
        idx = self.size < max_size

        opts = replace(self.opts, max_insert = max_size)

        return Fragments(
            self.range,
            self.start[idx],
            self.size[idx],
            opts
        )

    def __len__(self):
        return len(self.nfr)


@dataclass
class BinnedFragments:
    """
    Dataclass containing binned ATAC fragments within a range.
 
    Fragments are sorted into fixed bins (default 10bp) based on their midpoints,
    and classified based on insertion length (insertions less than 120 are 
    nucleosome-free by default). The insertion length is capped (default 300bp),
    since these fragments are less informative.
    """

    range: Range
    nfr: np.ndarray
    nuc: np.ndarray
    opts: Options = None

    @classmethod
    def from_fragments(cls, frags: Fragments, opts: Options=None):
        opts = get_options(frags.opts, opts)

        if not frags.range.align(opts.bin_size) == frags.range:
            raise ValueError("Fragments.range must be aligned to provided bin_size")

        n_bins = len(frags.range) // opts.bin_size

        nfr = np.zeros(n_bins, dtype=np.intc)
        nuc = np.zeros(n_bins, dtype=np.intc)

        midpoint = frags.start + frags.size // 2
        bin = (midpoint - frags.range.start + 1) // opts.bin_size

        in_bounds = (bin >= 0) & (bin < n_bins)

        is_nfr = frags.size < opts.nfr_cutoff

        is_nuc = frags.size >= opts.nfr_cutoff

        np.add.at(nfr, bin[in_bounds & is_nfr], 1)
        np.add.at(nuc, bin[in_bounds & is_nuc], 1)

        return cls(
            frags.range,
            nfr,
            nuc,
            opts
        )

    def to_data_frame(self) -> pd.DataFrame:
        bins = self.range.make_bins(self.opts.bin_size)

        bed = pd.DataFrame(
            [(r.chr, r.start, r.end) for r in bins],
            columns=["chr", "start", "end"]
        )

        bed["nfr"] = self.nfr
        bed["nuc"] = self.nuc

        return bed


class FragmentLoader:
    """
    Just needs iter -> Range method.
    """
    def __iter__(self) -> Iterator[Fragments]:
        pass
