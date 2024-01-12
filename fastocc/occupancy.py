from dataclasses import dataclass, replace

import numpy as np
import pandas as pd

from .options import Options, get_options
from .fragments import Fragments, BinnedFragments
from .range import Range
from .np_utils import rolling_window, smooth


def occupancy(nfr: np.ndarray, nuc: np.ndarray, opts: Options):
    """
    Calculates nucleosome occupancy from a vector of nucleosome-free and nucleosomal 
    reads, based on the formula (NFR*1.22) / (NFR + NUC).

    This function does not add a pseudocount, so if used outside the main script
    it can cointain NaNs in the output.
    """
    nfr_w = rolling_window(nfr, opts.window_size)
    nuc_w = rolling_window(nuc, opts.window_size)
    
    occ = (nfr_w * opts.nfr_ratio + opts.pseudocount_nuc) / (nfr_w + nuc_w + opts.pseudocount)

    # corrected values can't exceed 1
    return occ.clip(max=1)


@dataclass
class Occupancy:
    """
    Dataclass containing occupancy caclculated in a given genomic interval. Binned
    fragment counts are also stored for convenience.
    """
    range: Range
    nfr: np.ndarray
    nuc: np.ndarray
    occ: np.ndarray
    opts: Options

    def to_data_frame(self) -> pd.DataFrame:
        bins = self.range.make_bins(self.opts.bin_size)

        bed = pd.DataFrame(
            [(r.chr, r.start, r.end) for r in bins],
            columns=["chr", "start", "end"]
        )

        bed["nfr"] = self.nfr
        bed["nuc"] = self.nuc
        bed["occ"] = self.occ

        return bed

    @classmethod
    def from_fragments(cls, frags: Fragments, opts: Options=None):
        return cls.from_binned_fragments(BinnedFragments.from_fragments(frags, opts))

    @classmethod
    def from_binned_fragments(cls, fragments: BinnedFragments, opts: Options=None):
        opts = get_options(fragments.opts, opts)

        occ = occupancy(fragments.nfr, fragments.nuc, opts)

        if opts.smooth:
            occ = smooth(occ, opts.gaussian, opts.gaussian_norm)
        else:
            occ = occ[opts.shift:-opts.shift]

        s2 = opts.shift * 2

        r = replace(
            fragments.range.extend(-s2 * opts.bin_size),
            score = np.max(occ)
        )

        return cls(
            r,
            fragments.nfr[s2:-s2],
            fragments.nuc[s2:-s2],
            occ,
            opts,
        )
