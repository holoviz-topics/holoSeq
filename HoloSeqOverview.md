# HoloSeq

This project uses a precomputed mapping data format for sequence annotation that allows large scale data to be viewed
as 1D charts or 2D heatmaps using a generic visualisation infrastructure built using the [Holoviews ecosystem](https://holoviews.org/).

Preparing the coordinates for a 60GB HiC paf file with 720 million pairs takes a couple of hours so it makes sense to save them in
an intermediate precomputed format. Displaying the resulting 1.2GB compressed coordinates for 200M pairs only involving H1 takes about 10 minutes in comparison.

The presentation layer supports genome scale feature data associated with a genomic reference or other sequence. The user can pan and zoom 
smoothly from whole genomes down to individual points, with tens of millions of rows of data, in a web browser running on a suitable laptop or in Galaxy. 
The data format provides all the information needed for recreating a plot using the inbuilt reference sequence coordinates as axes.

The design isolates the complexities of displaying many different kinds of annotation at genomic scale, from the messy challenges of converting
complex existing data in standard formats. The intention is to allow any number of precomputed track coordinate files to be supplied to the generic display 
application, where they are automatically organised and displayed, using hints on layout supplied on the command line. 

The main use case envisioned is a central repository of pre-computed plots to make annotation of the VGP genomic data easily accessible.
Precomputed plot tracks can be reused indefinitely, mixed and matched by the user to suit their needs.

Each species has different coordinate systems so cannot share a reference sequence axis, but can each be displayed in
tracks with separate reference sequences side by side or stacked.

### Potential sources of annotation for display

HiC data in PAF format was used for the proof of concept 2D heatmaps. Mashmap approximate mapping PAF works well.

Bigwig, bed and GFF are the major formats for 1D annotation tracks.

### Coordinate system

Genomic reference data forms the backbone for any annotation browser.
New genomes are assembled into multiple "contigs", that are refined into chromosomes in reference genomes.
Genomes are typically handled in fasta format. A newly assembled haplotype may have thousands of contigs
that have not yet been merged into chromosomes. Contig names must be unique to each genome or haplotype.

For an interactive genome browser, tracks typically run horizontally, from the start of the reference sequence on the left to the last nucleotide of the last contig.

Holoviews dynamic maps require ordinal axis coordinates. Contigs must be ordered on the axes, usually by name or by length,
so they can be mapped as axes, with tick marks and labels.

Features must have a position to locate them in the reference sequence used to create the axis. 
Position is unambiguously described by the name of the contig, and the number of bases from the start of the contig to the start of the feature. 
Some features have a length while others for 2D grids are points. Many features may have optional annotation.

To convert feature positions into plot axis coordinates, contig lengths are cumulated in the order given, and a zero is inserted,
to give the ordinal position on the axis, for the first base of each ordered contig. 

When a feature is mapped, the appropriate contig's cumulated start is added to the feature offset, to give the 
ordinal start coordinate on that axis for the start of the feature.

Additional annotation values may be optionally displayed as hover tooltips.

## Input format for 1D and 2D features on pre-mapped axis coordinates

The converters produce gzip compressed text files. 

The text file must start with a header section, where every row begins with `@`.

The first row of the must be either `@v1HoloSeq1D [bar|scatter|line]` or `@v1HoloSeq2D`, or the data will not be processed.

For 1D data, the chart type may be one of `bar`, `scatter` or `line`. Default is `bar`. Regions with 4 or more SD above or below the global mean are 
emphasised

2D data will be presented as an autoscaling density heatmap. The header and data might only have 1 axis name, for example where HiC pairs from one haplotype are plotted
with that sequence on both axes, or 2 axis names, if HiC pairs involving both haplotypes, one on each axis, are being plotted.

The subsequent header rows must have the plot title, plot type, axis names, contig names and their cumulated lengths, delimited by whitespace, and starting with `@` such as

```
@v1HoloSeq1D bar
@@title a very small bar plot
@@xclenfile hg002_H1suffixed.len
@H1 chr1 0
@H1 chr2 2500
@H1 chr3 7000
@H1 chr4 9500
```

In this case, H1 will show four chromosomes starting at each of the positions shown.

Metadata such as the name of the chromosome lengths file and the plot title, is prefixed with `@@`


Data rows for 1D data must have the x ordinal axis coordinate, a feature length, and an annotation value to show for that length, such as:
`2455453443 128 99.8`

Data rows for 2D data must have the x and y ordinal grid coordinates that correspond to the contigs in the header, and 
may have a heatmap value or other annotation such as:

`3999543 58898548  2`

If there are multiple rows at the same coordinates without annotation, the count is used to provide heatmap values.

The header allows the cumulated start axis coordinate of each ordered contig to be recovered, and used to create the axis tick marks 
corresponding to the start of each contig.

## Visualisation

Input files with identical contig names are grouped into a vertical stack with linked axes.
For each of these stacks, the header contig/length values are used to prepare the axis tick marks.

For each 1D input, the chart type is prepared as a panel row.
2D inputs are turned into datamaps with optional scale bars and tap for location coordinates.

Groups of tracks are stacked into a panel to be served.

All plots have optional scale bars for reference and a mouse click anywhere on a plot will show the coordinates above the plot.
Tap coordinates are calculated from a stream giving the tap x and y coordinates, using a binary search on the contig coordinate starts.

## Deployment

A visualisation can be run locally and viewed in a desktop browser by stacking any number of compressed coordinate hseq 
files 

`panel serve holoseq_display.py --show --args --inFile foo.gz bar.gz baz.gz zot.gz --title my holoseq plot`

The Dockerfile is a work in progress as part of an interactive Galaxy tool.


