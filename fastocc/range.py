from dataclasses import dataclass, replace, astuple
from typing_extensions import Self
from tempfile import NamedTemporaryFile
from typing import Iterator

from pybedtools import BedTool

from .options import Options, get_options


@dataclass(frozen=True)
class Range:
    """
    Class representing a genomic range or interval.

    NB this uses "GRanges" style inclusive ranges, not BED-style insanity.
    """

    chr: str
    start: int
    end: int
    name: str = "."
    score: float = 0
    strand: str = "."

    def extend(self, bp: int) -> Self:
        """
        Return a range extended by *bp* on either end. Does not check bounds.
        """

        s = self.start - bp
        e = self.end + bp

        return replace(self, start=s, end=e)

    def align(self, bin_size: int) -> Self:
        """
        Return a window aligned to fixed-width bins.

        e.g. with 10bp bins: chr2:54-89 => chr2:51-90
        """

        s = (self.start - 1) // bin_size * bin_size + 1
        e = (self.end - 1) // bin_size * bin_size + bin_size

        return replace(self, start=s, end=e)

    def make_bins(self, bin_size: int) -> list[Self]:
        """
        Return a list of Ranges representing individual bins.

        NB These are aligned to a genomic 'grid' based on bin_size
        e.g. make_bins(chr2:54-84) => [Range(chr2, 51, 60), ..., Range(chr2, 81, 90)]
        """

        r = self.align(bin_size)
        starts = range(r.start, r.end, bin_size)

        return [
            Range(self.chr, start=s, end=s + bin_size - 1)
            for s in starts
        ]

    def to_bed(self) -> tuple:
        return astuple(replace(self, start=self.start-1))

    def __len__(self) -> int:
        """
        Return the number of base pairs included in the range. Behaviour is similar to
        GRanges, so the interval "chr1:31-40" has length 10, "chr1:25-25" has length 1.
        """
        return self.end - self.start + 1
    
    def coords(self) -> str:
        """
        String epresentation of range coordinates e.g. "chr1:53-70"
        """
        return f"{self.chr}:{self.start}-{self.end}"


def write_chrom_sizes(chrom_sizes: dict[str, int], f):
    """
    Write chromosome sizes to a file in UCSC format. Format is the same as
    the output of `fetch_chrom_sizes`, ie [("chr1", 248956422), ...]
    """
    for chr, size in chrom_sizes.items():
        f.write(f"{chr}\t{size}\n")
    
    f.flush()


@dataclass
class BedFile:
    """
    Dataclass wrapping `pybedtools.BedTool`, iterating over the class will yield
    Range instances (with the correct start coordinate).
    
    A constructor from a bed file is provided to ensure ranges are extended, sorted and 
    non-overlapping, and aligned to genomic bins, based on parameters provided in
    Options. This requires chromosome information as well.

    BedFile should be possible to replace with any classes which provides an iterator
    over ranges.
    """
    bed: BedTool
    opts: Options

    def __iter__(self) -> Iterator[Range]:
        for i in self.bed:
            score = i.score if i.score != "." else 0
            r = Range(i.chrom, i.start+1, i.end, i.name, score, i.strand)
            yield r.align(self.opts.bin_size)

    @classmethod
    def open(cls, bed: str, chrom_info: dict[str, int], opts: Options=None):
        opts = get_options(opts)

        with NamedTemporaryFile('w') as genome:
            write_chrom_sizes(chrom_info, genome)
        
            fixed = BedTool(bed).slop(l=opts.extend, r=opts.extend, g=genome.name).sort(g=genome.name).merge(d=opts.merge)

        return cls(
            fixed,
            opts
        )

    @classmethod
    def make_windows(cls, chrom_info: list[tuple[str, int]], size=10000, opts=None):
        """
        Wrap window_maker function to use chromsizes.

        NB chromsizes are created in order, so no need to pipe into sort and so we
        avoid PR#380
        """
        opts = get_options(opts)

        g = {chr: (0, end) for chr, end in chrom_info}
        windows = BedTool().set_chromsizes(g).window_maker(w=size)

        return cls(
            windows,
            opts
        )
    
