from dataclasses import dataclass
import json

import numpy as np
from copy import copy

from .np_utils import gaussian

@dataclass
class Options:
    """
    This dataclass facilitates passing around large numbers of somewhat-related params
    to functions in fastocc. In general, you should set this once, when loading fragments,
    and let these values cascade through the program. This avoids errors where, for 
    example, longer reads are read in, but then not accounted for properly in windowing
    functions.
    """
    # bins + windows
    window_size: int = 5 # number of bins for window estimation
    bin_size: int = 10 # width of bins
    extend: int = 500 # base pairs to extend input ranges
    merge: int = 200 # merge input ranges closer than this (after extension)
    # ATAC fragment options
    max_insert: int = 300
    shift_atac: bool = True
    filter_duplicates: bool = True
    # occupancy params
    nfr_cutoff: int = 119 # threshold above which reads are considered to come from nucleosomes
    nfr_ratio: float = 2.0 # ratio of overall NFR reads to reads < nfr_cutoff
    pseudocount: int = 20 
    pseudocount_nuc: int = 1
    # smoothing params
    smooth: bool = True
    gaussian_sd: float = 1.5
    # nfr params
    nfr_threshold: float = 0.6
    # nfr_min: float = 0.8
    # nfr_extend: int = 20

    def __post_init__(self):
        self.shift = self.window_size // 2
        self.gaussian = gaussian(self.window_size, self.gaussian_sd)
        self.gaussian_norm = np.convolve(self.gaussian, np.ones(self.window_size), mode="valid")

    @classmethod
    def load(cls, f: str):
        jsn = json.load(open(f))

        return cls(**jsn)


class NucPosOptions(Options):
    """
    Alternative options to visualise nucleosome positions rather than NFRs, using
    larger windows, a lower NFR ratio, and much weaker regularization towards 0.5
    (not zero).

    The output is best visualised as -2*(occ-0.5), with 1 representing a well-positioned
    nucleosome, -1 an NFR, and 0 no data (or a disordered/atypical region).
    """
    window_size: int = 11
    nfr_ratio: float = 1.2
    pseudocount: int = 2
    pseudocount_nuc: int = 1
    smooth: bool = True
    gaussian_sd: float = 2


def get_options(*args):
    """
    Convenience function to avoid boilerplate for functions which take an Options
    argument, and another argument which supplies its own opts (e.g. BinnedFragments).

    The rightmost non-None argument is returned, or else default options are used.
    """
    opts = [arg for arg in args if arg is not None]

    if len(opts) == 0:
        return Options()

    return copy(opts[-1])
