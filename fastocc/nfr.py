from typing import Iterator

import numpy as np

from .np_utils import nonzero_intervals
from .occupancy import Occupancy
from .range import Range

def threshold(occ: Occupancy, th=0.6, th_min=0.8) -> Iterator[Range]:
    """
    Alternative th function to expand range if minimum reached.
    """
    coords = nonzero_intervals(occ.occ > th)

    r = occ.range

    i = 0

    for start, stop in coords:
        name = r.coords() if r.name == "." else r.name

        th_max = np.max(occ.occ[start:stop])

        if not th_max > th_min:
            continue

        nfr_reads = np.sum(occ.nfr[start:stop])

        i += 1

        yield Range(
            r.chr,
            r.start + occ.opts.bin_size*(start),
            r.start + occ.opts.bin_size*stop - 1,
            f"{name}:NFR{i}",
            nfr_reads,
            "."
        )
