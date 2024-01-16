import csv
import sys

import numpy as np

import pyBigWig

from .range import Range


class BigWigWriter:
    """
    Class to write binned range data to a BigWig file. Writes must be in order
    (chromosome and position).
    """
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


class BedGraphWriter:
    """
    Class to write binned range data to a BedGraph file.
    """
    def __init__(self, f: str, bin_size):
        self.file = open(f, 'w')
        self.writer = csv.writer(self.file, delimiter="\t")
        self.bin_size = bin_size

    def write(self, r: Range, *data):
        for i, b in enumerate(r.make_bins(self.bin_size)):
            row = [b.chr, b.start - 1, b.end] + [col[i] for col in data]

            self.writer.writerow(row)
