from typing import Any
import h5py
import numpy as np

from itertools import groupby
from dataclasses import asdict
from pathlib import Path

from ncls import NCLS

from .options import Options
from .range import Range
from .io import BigWigWriter, BedGraphWriter


def nbins(r: Range, size):
    return len(r) // size


def write_datasets(grp: h5py.Group, datasets: dict[str, np.ndarray], **kwargs):
    for k, d in datasets.items():
        grp.create_dataset(k, data=d, **kwargs)


class HDF5Writer:
    compression = "gzip"
    chunks = 10000
    scaleoffset = 2 # decimal places for floats

    def __init__(self, file: h5py.File, schema: dict, opts: Options=None, genome_info=None):
        self.file = file
        self.file.create_group("chromosomes")

        meta = self.file.create_group("metadata")
        
        meta.create_group("schema")
        meta.create_group("opts")

        if opts:
            meta["opts"].attrs.update(asdict(opts))
            schema["bin_size"] = opts.bin_size

        if genome_info:
            meta.create_group("genome_info")
            names, sizes = zip(*genome_info.items())
            write_datasets(meta["genome_info"], {"names": names, "sizes": sizes})

        self.cols = schema["cols"]
        self.bin_size = schema["bin_size"]

        meta["schema"].attrs.update(schema)
        self.schema = schema

    def write(self, range_data):
        data_by_chr = groupby(range_data, key=lambda fs: fs.range.chr)

        for chr, windows in data_by_chr:
            self.write_chr(chr, windows)

    def write_chr(self, chr, windows):
        """
        A chromosome is an atomic write because of how indexes are calculated and
        the fact that windows are required in order.
        """
        grp = self.file["chromosomes"].create_group(chr)
        grp.create_group("ranges")
        grp.create_group("data")

        windows = list(windows) # don't exhaust iterator, mem in tiny anyway

        ranges, data = self.encode_as_arrays(windows)

        # fix for empty // small chrs:
        chunks = min(len(data[self.cols[0]]), self.chunks)

        print(chunks.__repr__())

        write_datasets(grp["ranges"], ranges)

        # write float arrays with scaleoffset
        for col, type in zip(self.cols, self.schema["types"]):
            if type[:5] == "float" and self.scaleoffset:
                grp["data"].create_dataset(col, data=data[col], chunks=chunks, scaleoffset=self.scaleoffset)
                data.pop(col)

        write_datasets(grp["data"], data, chunks=chunks, compression="gzip")

        # print(f"{chr} written.")

    def encode_as_arrays(self, windows) -> (dict, dict):
        """
        Encode a list of binned range data (e.g. BinnedFragments, Occupancy) into
        numpy arrays representing the ranges and the data.
        """
        rs = (w.range for w in windows)
        start, end, bins = zip(
            *[(r.start, r.end, nbins(r, self.bin_size)) for r in rs]
        )

        total = sum(bins)

        idx = np.concatenate(([0], np.cumsum(bins)[:-1]), axis=0)

        ranges = {
            "start": start,
            "end": end,
            "idx": idx
        }

        data = {}
        
        for col, typ in zip(self.schema["cols"], self.schema["types"]):
            data[col] = np.zeros(total, dtype=np.dtype(typ))

        for i, w in enumerate(windows):
            end = idx[i] + bins[i]
            for col in self.schema["cols"]:
                # TODO add proper error for wrong size: current check is implicit
                data[col][idx[i]:end] = getattr(w, col)

        return ranges, data

    @classmethod
    def open(cls, f: str, *args, **kwargs):
        if Path(f).exists():
            raise IOError(f"File {f} exists already, open manually to overwrite.")

        return cls(h5py.File(f, 'w'), *args, **kwargs)


class HDF5Loader:
    data_group = "chromosomes"

    def __init__(self, f: h5py.File, cls=None):
        self.file = f
        self.cls = cls

        meta = self.file["metadata"]

        self.schema = dict(meta["schema"].attrs)
        self.bin_size = self.schema["bin_size"]
        self.cols = self.schema["cols"]

        self.index = {}

        # only load relevant fields
        if self.cls and hasattr(self.cls, "__dataclass_fields__"):
            self.cols = [col for col in self.schema["cols"] if col in self.cls.__dataclass_fields__]
        else:
            self.cols = self.schema["cols"]

        if "genome_info" in meta:
            self.chrom_info = {
                str(k, "utf-8"): v for k, v in
                zip(
                    meta["genome_info"]["names"],
                    meta["genome_info"]["sizes"]
                )
            }
        else:
            self.chrom_info = None

        if "opts" in meta:
            self.opts = Options(**meta["opts"].attrs)
        else:
            self.opts = None

    @classmethod
    def open(cls, f: str, *args, **kwargs):
        return cls(h5py.File(f, 'r'), *args, **kwargs)

    def __iter__(self):
        # this maintains correct order if chrom_info is available
        chrs = self.chrom_info or self.file[self.data_group]

        for chr in chrs:
            if chr in self.file[self.data_group]:
                yield from self.iter_chrom(chr)

    def build_index(self, chrs=None):
        if not chrs:
            chrs = self.file[self.data_group]

        for chr in chrs:
            grp = self.file[self.data_group][chr]

            self.index[chr] = NCLS(
                grp["ranges"]["start"][:],
                grp["ranges"]["end"][:],
                np.arange(len(grp["ranges"]["start"]))
            )

    def find_overlaps(self, r: Range):
        if not r.chr in self.index:
            self.build_index([r.chr])

        return self.index[r.chr].find_overlap(r.start, r.end)

    def iter_chrom(self, chr):
        grp = self.file[self.data_group][chr]

        ranges = [
            Range(chr, s, e) for s, e in zip(grp["ranges"]["start"], grp["ranges"]["end"])
        ]

        # load info to avoid chunk nonsense
        data = {col: grp["data"][col][:] for col in self.cols}

        idx = np.concatenate((grp["ranges"]["idx"], [len(data[self.cols[0]])]), axis=0)

        for i, r in enumerate(ranges):
            frm = idx[i]
            to = idx[i+1]

            result = {
                "range": r,
                "data": {col: data[col][frm:to] for col in self.cols}
            }

            yield self.encode(result)

    def iter_range(self, q: Range):
        grp = self.file[self.data_group][q.chr]

        for s, e, i in self.find_overlaps(q):
            r = Range(q.chr, s, e)
            frm = grp["ranges"]["idx"][i]
            to = frm + nbins(r, self.bin_size)

            yield self.encode({
                "range": r,
                "data": {col: grp["data"][col][frm:to] for col in self.cols}
            })

    def encode(self, result: dict):
        if self.cls:
            return self.cls(range=result["range"], opts=self.opts, **result["data"])
        else:
            return result

    def export_bigwig(self, bw: str, track: str):
        bigwig = BigWigWriter(bw, self.chrom_info, self.bin_size)

        for r in self:
            bigwig.write(r["range"], r["data"][track].astype(np.float32))

        bigwig.file.close()

    def export_bedgraph(self, bg: str, regions: list[Range]=None):
        bedgraph = BedGraphWriter(bg, self.bin_size)

        if regions:
            for r in regions:
                for q in self.iter_range(r):
                    bedgraph.write(q["range"], *q["data"].values())
        else:
            for r in self:
                bedgraph.write(r["range"], *r["data"].values())
