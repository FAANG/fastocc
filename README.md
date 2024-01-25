# fastocc

Rapid and robust estimation of nucleosome-free regions and nucleosome positioning.

<!-- ## Introduction -->

<!-- TODO very nice picture -->

## Quickstart

Most of the core functionality in `fastocc` is accessible through the CLI with a small
number of commands.

The bulk of processing is done by the `call` command, which takes a `.bam` file containing
paired ATAC-seq reads, processes it into fragments, and then counts those fragments into
nucleosomal (longer) and nucleosome-free (shorter) bins along the genome (default is 10bp).

`fastocc` can be run in tiles over the whole genome, but the algorithm will generally
be uninformative outside of regions of high ATAC-seq signal (ie peaks), specified with
the `--bed` argument.

```
fastocc call --bam my_atac.bam --bed peaks.bed --output fragments.h5wig
```

`.h5wig` is a novel format (built on HDF5), which stores a set of regions (e.g. peaks),
and acro each region tracks multiple signals in consecutive bins across the whole region.
This concept of binned data in ranges is used internally in many places in `fastocc`,
and is a significant reason behind its efficiency.

We can then call nucleosome-free regions using the `nfr`.

```
fastocc nfr --fragments fragments.h5wig --output nfrs.bed
```

If we want to visualise the predicted occupancy in a genome browser (or the counts of
NF for nucleosomal reads), we can export these from the fragments file as a `.bigWig`:

```
fastocc export --fragments fragments.h5wig --track occ --output occ.bw
```

Lastly, we can export everything to an extended `.bedGraph` format. This isn't recommended
for the whole dataset due to the resulting file size, so a set of regions can be provided:

```
fastocc dump --fragments fragments.h5wig --bed regions_of_interest.bed --output fragments.bg
```

If you do need programmitic random access to the fragment-level data, functions are
provided in Python (exported in the package) and R (in `extra/`). Implementing reading
should be straightforward in any other language with a mature HDF5 library.

## What happens next?

More information about the package can be found in the documentation and source code,
and the Jupyter notebooks included in the repository.

Some things you might be interested in:
- Merge NFRs across cell lines to generate precises consensus cis-regulatory elements.
- Scan NFRs for enriched motifs.
- Predict TF-binding based on ChIP-seq and NFRs
- Use the internal API to analyse scATAC data directly from `.arrow` files.

If you do any of this, please let me know! Github DMs are the best way to reach me,
see [here](https://github.com/mgperry). If you have any issues or feature requests,
I will try to watch issues here but may require a nudge. I am no longer working on this
full time, so you may be asked to contribute, please don't be scared to do this!

## Acknowledgements

The development of `fastocc` took place within the 
[Ensembl Regulation](http://www.ensembl.org/info/genome/funcgen/index.html) team at
EMBL-EBI. The project received funding from the European Union’s Horizon 2020
programme under the [AQUA-FAANG](https://www.aqua-faang.eu/)
and [GENE-SWitCH](https://www.gene-switch.eu/) projects.

The tool itself was based on ideas originally implemented in the NucleoATAC package, if
you're interested in the biological reasoning behind `fastocc` you sohuld read their
paper in [Genome Research](https://genome.cshlp.org/content/25/11/1757) and check out
NucleoATAC on [GitHub](https://github.com/GreenleafLab/NucleoATAC).

`fastocc` is available under the Apache v2.0 license, included in the repository.
