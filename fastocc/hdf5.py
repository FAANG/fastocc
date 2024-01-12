import h5py
import numpy as np

from itertools import groupby
from dataclasses import asdict
from pathlib import Path

from .options import Options
from .range import Range
from .io import BigWigWriter


def count_bins(r: Range, size):
    return (r.end - r.start + 1) // size


def write_datasets(grp: h5py.Group, datasets: dict[str, np.ndarray], **kwargs):
    for k, d in datasets.items():
        grp.create_dataset(k, data=d, **kwargs)


class HDF5Writer:
    compression = "gzip"
    chunks = 10000
    scaleoffset = 2 # decimal places for floats

    def __init__(self, file: h5py.File, schema: dict, opts: Options=None, genome_info=None):
        self.file = file
        self.schema = schema

        self.file.create_group("chromosomes")

        meta = self.file.create_group("metadata")
        
        meta.create_group("schema")
        meta["schema"].attrs.update(schema)


        meta.create_group("opts")

        if opts:
            meta["opts"].attrs.update(asdict(opts))
            self.bin_size = opts.bin_size
        else:
            bin_size = schema.get("bin_size", None)

            if not bin_size:
                raise ValueError("bin_size must be given in opts or schema")

            meta["opts"].attrs.update({"bin_size": bin_size})
            self.bin_size = bin_size

        if genome_info:
            meta.create_group("genome_info")
            names, sizes = zip(*genome_info.items())
            write_datasets(meta["genome_info"], {"names": names, "sizes": sizes})

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

        windows = list(windows) # don't exhaust iterator, mem in tiny anyway

        ranges, data = self.encode(windows)

        # fix for empty // small chrs:
        total = sum(ranges["bins"])

        chunks = min(total, self.chunks)

        write_datasets(grp, ranges)

        # write float arrays with scaleoffset
        for col, type in zip(self.schema["cols"], self.schema["types"]):
            if type[:5] == "float" and self.scaleoffset:
                grp.create_dataset(col, data=data[col], chunks=chunks, scaleoffset=self.scaleoffset)
                data.pop(col)

        write_datasets(grp, data, chunks=chunks, compression="gzip")

        # print(f"{chr} written.")

    def encode(self, windows) -> (dict, dict):
        """
        Encode a list of binned range data (e.g. BinnedFragments, Occupancy) into
        numpy arrays representing the ranges and the data.
        """
        rs = (w.range for w in windows)
        start, end, bins = zip(
            *[(r.start, r.end, count_bins(r, self.bin_size)) for r in rs]
        )

        total = sum(bins)

        ranges = {
            "start": start,
            "end": end,
            "bins": bins
        }

        data = {}
        
        for col, typ in zip(self.schema["cols"], self.schema["types"]):
            data[col] = np.zeros(total, dtype=np.dtype(typ))

        start = 0

        for i, w in enumerate(windows):
            end = start + bins[i]
            for col in self.schema["cols"]:
                data[col][start:end] = getattr(w, col)
            start = end

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

        if "genome_info" in meta:
            self.chrom_info = {
                str(k, "utf-8"): v for k, v in
                zip(
                    meta["genome_info"]["names"],
                    meta["genome_info"]["sizes"]
                )
            }

        if "opts" in meta:
            self.opts = Options(**meta["opts"].attrs)
        else:
            self.opts = None

    @classmethod
    def open(cls, f: str, *args, **kwargs):
        return cls(h5py.File(f, 'r'), *args, **kwargs)

    def __iter__(self):
        # this maintains correct order if chrom_info is available
        if self.chrom_info:
            for chr in self.chrom_info:
                if chr in self.file[self.data_group]:
                    yield from self.decode(chr)
        # fall back to iterating over group
        else:
            for chr in self.file[self.data_group]:
                yield from self.decode(chr)

    def decode(self, chr):
        grp = self.file[self.data_group][chr]

        ranges = [
            Range(chr, s, e) for s, e in zip(grp["start"], grp["end"])
        ]

        # only load relevant fields
        if self.cls and hasattr(self.cls, "__dataclass_fields__"):
            cols = [col for col in self.schema["cols"] if col in self.cls.__dataclass_fields__]
        else:
            cols = self.schema["cols"]

        # load info to avoid chunk nonsense
        datasets = {col: grp[col][:] for col in cols}

        start = 0

        for r, n in zip(ranges, grp["bins"][:]):
            end = start + n

            data = {"range": r}

            for col in cols:
                data[col] = datasets[col][start:end]

            if self.cls:
                yield self.cls(**(data | {"opts": self.opts}))
            else:
                yield data

            start = end

    def export_bigwig(self, bw, track):
        bigwig = BigWigWriter(bw, self.chrom_info, self.opts.bin_size)

        for r in self:
            bigwig.write(r["range"], r[track].astype(np.float32))

        bigwig.file.close()
