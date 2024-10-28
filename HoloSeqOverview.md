## HoloSeq

This document describes a data format, containing pre-gridded sequence annotation, that allows large scale data to be viewed
as 1d charts or 2d heatmaps using a generic visualisation infrastructure built using the Python [Holoviews ecosystem](https://holoviews.org/).

The presentation layer is designed for large scale data associated with a genomic reference or other sequence. It can pan and zoom smoothly from whole genomes down to 
individual points with tens of millions of rows of data, in a web browser running on a suitable laptop or in Galaxy.

### Coordinate system

For an interactive genome browser, tracks typically run horizontally, from the start of the reference sequence on the left to the last nucleotide of the last contig.

Holoviews dynamic maps require ordinal axis coordinates while a typical assembly haplotype or reference sequence consists of
an arbitrary number of chromosomes, or more generally `contigs`, each of a fixed length. The contigs must be ordered, usually by name or by length,
so they can be placed in a row to create the horizontal axis and for 2D grids, the vertical axis, each with tick marks indicating the start of each contig.

The features annotating any track must have a position in the reference sequence used to create the axis. 
It is usually described by the name of the contig, and the number of bases from the start of the contig to the start of the feature. 
Some features have a length while others for 2D grids are points. Features may have additional annotation for user toggled hover tool tips. 

To convert feature positions into plot axis coordinates, contig lengths are cumulated in the order given, to give the ordinal axis value for the first base
of each ordered contig. When a feature is mapped, the contig cumulated start is added to the feature offset, to give the 
ordinal coordinate in that axis.

In addition to either one or two coordinates, features may have additional annotation values, optionally displayed as tooltips.

## Input format for 1D and 2D features on pre-mapped axis coordinates

The input file may be a gzip in which case it will be uncompressed.

The text data must start with a header section describing the data, where every row starts with `@`.

The first row of the (uncompressed) header must be either `@v1HoloSeq1D [chart type | bar]` or `@v1HoloSeq2d` or the data will not be processed.

For 1D data, the chart type may be one of bar, scatter or line. Default is bar. Regions with 4 or more SD above or below the global mean are 
emphasised

2D data will be presented as an autoscaling density heatmap, and the data may have 1 axis name if HiC pairs from 1 haplotype are being plotted
where that sequence is on both axes, or 2 axis names if HiC pairs that involve both haplotypes are being plotted.

The subsequent header rows must have the axis names, contig names and their lengths, delimited by whitespace, and starting with `@` such as

```
@H1 chr1 2000
@H2 chr1 2500
@H1 chr2 1000
@H2 chr2 1500
```

In this case, H1 will show chr1 and chr2 tick marks on the 3000 base horizontal axis, while H2 coordinates will be placed on the 4000 base vertical axis
with chr1 and chr2 tick marks.

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

The visualisation can be run locally and viewed in a desktop browser using 

`panel serve [notebookname].ipynb --args foo=23 bar=./readme.paf`

