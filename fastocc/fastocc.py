import click
import csv
import os

from .occupancy import Occupancy
from .options import Options, NucPosOptions
from .bam import BamFragmentLoader
from .hdf5 import HDF5Writer, HDF5Loader
from .nfr import threshold
from .range import BedFile, BedTool


@click.group()
def fastocc():
    pass


@click.command()
@click.option("--bam", help="BAM file containing ATAC-seq reads", required=True, type=str)
@click.option("--bed", help="BED file containing ranges.", type=str, required=True)
@click.option("-o", "--output", help="Output file location.", type=str, required=True)
@click.option("--options", help="Option file location.", type=str, default=None)
@click.option("--nps", help="Calculate Nucleosome positions instead of NFRs", is_flag=True)
def call(bam: str, bed: str, output: str, options: str, nps: bool):
    if nps:
        defaults = NucPosOptions
    else:
        defaults = Options

    if options:
        opts = defaults.load(options)
    else:
        opts = defaults()

    if os.path.exists(output):
        IOError(f"Output {output} already exists.")

    fragments = BamFragmentLoader.open(bam, bed, opts)

    schema = {
        "cols": ["nfr", "nuc", "occ"],
        "types": ["int16", "int16", "float32"],
        "bin_size": opts.bin_size
    }

    writer = HDF5Writer.open(output, schema, opts, fragments.bam.chrom_info)

    writer.write(Occupancy.from_fragments(fs) for fs in fragments)


@click.command
@click.option("--fragments", help="HDF5 file containing `fastocc call` output.", required=True, type=str)
@click.option("--output", help="Output file location.", type=str, required=True)
@click.option("--options", help="Option file location.", type=str)
def nfr(fragments: str, output: str, options: str):
    if options:
        opts = Options.load(options)
    else:
        opts = Options()

    h5 = HDF5Loader.open(fragments, Occupancy)

    fh = open(output, 'w')
    writer = csv.writer(fh, delimiter="\t")

    for occ in h5:
        nfrs = threshold(occ, opts.nfr_threshold)

        writer.writerows(nfr.to_bed() for nfr in nfrs)


@click.command
@click.option("--fragments", help="HDF5 file containing `fastocc call` output.", required=True, type=str)
@click.option("--output", help="Output file location.", type=str)
@click.option("--track", help="Track to export.", type=str, default="occ")
def export(fragments: str, track: str, output: str):
    h5 = HDF5Loader.open(fragments)

    h5.export_bigwig(output, track)

@click.command
@click.option("--fragments", help="HDF5 file containing `fastocc call` output.", required=True, type=str)
@click.option("--bed", help="Optional list of regions to export.", type=str, default=None)
@click.option("--output", help="Output file location.", type=str, required=True)
def dump(fragments: str, bed: str, output: str):
    h5 = HDF5Loader.open(fragments)

    opts = Options(bin_size=h5.bin_size)

    if bed:
        regions = BedFile(BedTool(bed), opts)
    else:
        regions = None

    h5.export_bedgraph(output, regions)


fastocc.add_command(call)
fastocc.add_command(nfr)
fastocc.add_command(export)
fastocc.add_command(dump)
