import csv
import sys

import numpy as np

import pyBigWig

from .range import Range


class BigWigWriter:
    def __init__(self, f: str, chrom_info: dict[str, int], bin_size=10):
         self.file = pyBigWig.open(str(f), 'w')
         self.file.addHeader(list(chrom_info.items()))
         self.bin_size = bin_size

    def write(self, r: Range, data: np.ndarray):
        self.file.addEntries(
            r.chr,
            r.start-1,
            values=data,
            span=self.bin_size,
            step=self.bin_size,
        )
