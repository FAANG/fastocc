import sys
from typing import Iterator
from dataclasses import dataclass

import pysam
import numpy as np

from .range import Range, BedFile
from .options import Options, get_options
from .fragments import FragmentLoader, Fragments


def load_fragments(reads: Iterator[pysam.AlignedRead], opts: Options) -> tuple[np.ndarray,np.ndarray]:
    """
    Reconstruct full-length ATAC-seq fragments from Bam reads. Filters for the first
    read of proper pairs, corrects the Tn5 insertion distance and returns the start and
    length (size) of a fragment in a named tuple. By default, fragments with identical
    start/end coordinates are considered to be duplicates and filtered.
    """

    opts = get_options(opts)

    start = []
    size = []

    for read in reads:
        if read.is_proper_pair and not read.is_reverse:
            if opts.shift_atac:
                #get left position
                left_insertion = read.pos + 4
                # insert size corrected by 8 base pairs to be insertion to insertion
                frag_size = read.template_length - 8
            else:
                left_insertion = read.pos
                frag_size = read.template_length

            if frag_size > opts.max_insert:
                continue

            start.append(left_insertion)
            size.append(frag_size)

    if opts.filter_duplicates:
        unique = np.unique(np.column_stack([start, size]), axis=0)
        start = unique[:,0]
        size = unique[:,1]

    return np.array(start, dtype=int), np.array(size, dtype=int)


def fetch_chrom_sizes(h: pysam.AlignmentHeader) -> dict[str, int]:
    """
    Extract chromosome sizes from BAM header (via pysam).

    h["SQ"] => [{"SN": "chr1", "LN", 248956422}, ...]
    output => {"chr1": 248956422, ...}
    """
    return {seq["SN"]: seq["LN"] for seq in h['SQ']}


class Bam:
    """
    Wrapper around pysam.AlignmentFile to return full-length ATAC-seq fragments
    rather than individual reads.

    Also loads genome info and checks bounds on queries.

    NB By default, Bam files are indexed if there is none present.
    """
    def __init__(self, bam: pysam.AlignmentFile, index=True):
        if index and not bam.has_index():
            print("No BAM index found, indexing BAM file...", file=sys.stderr)
            fn = bam.filename.decode()
            bam.close()

            pysam.index(fn)
            bam = pysam.AlignmentFile(fn) # reload index
        
        self.bam = bam
        self.chrom_info = fetch_chrom_sizes(bam.header)

    def fetch(self, r: Range) -> Iterator[pysam.AlignedRead]:
        """
        Wrap pysam so that extended ranges don't fail due to out-of-bounds
        """
        return self.bam.fetch(
            r.chr,
            max(r.start, 1),
            min(r.end, self.chrom_info[r.chr])
        )
    
    def fragments(self, r: Range, opts: Options) -> Fragments:
        """
        Capture all fragments less than opts.max_insert overlapping a given interval.
        """
        opts = get_options(opts)

        reads = self.fetch(r.extend(opts.max_insert))
        start, size = load_fragments(reads, opts)

        return Fragments(
            r,
            start,
            size,
            opts
        )
    
    @classmethod
    def open(cls, bam_file: str, index=True):
        return cls(pysam.AlignmentFile(bam_file), index)


@dataclass
class BamFragmentLoader(FragmentLoader):
    """
    Dataclass containing a Bam file and a associated set of ranges, providing an 
    iterator of ATAC-seq fragments at each range.

    `open` method handles creation of a BamFragmentLoader from filenames, automatically
    extracting genome info and passing this to BedFile.
    """
    bam: Bam
    ranges: Iterator[Range]
    opts: Options

    def __iter__(self) -> Iterator[Fragments]:
        for r in self.ranges:
            r_ext = r.extend(2*self.opts.shift*self.opts.bin_size) # extend for windowing // smoothing

            yield self.bam.fragments(r_ext, self.opts)
    
    @classmethod
    def open(cls, bam: str, bed: str, opts: Options, index=True):
        bam_file = Bam.open(bam, index=index)

        ranges = BedFile.open(bed, bam_file.chrom_info, opts)

        return cls(bam_file, ranges, opts)
