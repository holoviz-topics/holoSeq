## HoloSeq

This document describes a data format, containing pre-gridded sequence annotation, that allows large scale data to be viewed
as 1d charts or 2d heatmaps using the Python Holoviews ecosystem.

The presentation layer is designed for data associated with a genomic reference or other sequence.

### Coordinate system

Holoviews dynamic maps require ordinal axis coordinates. A reference sequence typically consists of
an arbitrary number of chromosomes, or more generally `contigs`, each of a fixed length.

The annotation for any track to be displayed, has a position in the sequence, given by the name of the contig and the number of 
bases from the start of the contig to the start of the feature. 

To convert these into ordinal coordinate, contig lengths are cumulated, to give the ordinal axis value for the first base
of each of the ordered contigs. When a feature is mapped, the contig cumulated start is added to the feature offset, to give the 
ordinal coordinate in that axis.

In addition to either one or two coordinates, features may have additional annotation values to be displayed as tooltips.


## Input format

The input file may be a gzip in which case it will be uncompressed.

The text data must start with a header section describing the data, where every row starts with `@`.

The first row of the (uncompressed) header must be either `@v1HoloSeq1D [chart type]` or `@v1HoloSeq2d` or the data will not be processed.

For 1D data, the chart type may be one of bar, scatter or line. Default is bar. Regions with 4 or more SD above or below the global mean are 
emphasised
.
2D data will always be presented as a heatmap, and the data may have 1 axis name if HiC pairs from 1 haplotype are being plotted
where that sequence is on both axes, or 2 axis names if HiC pairs that involve both haplotypes are being plotted.

The subsequent header rows must have the axis names, contig names and their lengths, delimited by whitespace, and starting with `@` such as
`
@H1 chr1 199943902
@H2 chr1 199942264
...
`

Data rows for 1d data must have the x ordinal axis coordinate, a feature length, and an annotation value to show for that length, such as
`2455453443 128 99.8`

Data rows for 2d data must have the x and y ordinal grid coordinates and optional annotation such as a heatmap value
`3999543 58898548  2`

If there are multiple rows at the same coordinates without annotation, the count is used to provide heatmap values.

## Visualisation

Input files with identical contig names are grouped into a vertical stack with linked axes
For each of these stacks, the header contig/length values are used to prepare the axis tick marks.

For each 1D input, the chart type is prepared as a panel row.
2D inputs are turned into datamaps with optional scale bars and tap for location coordinates.

All the tracks are stacked into a panel to be served.

## Deployment

The visualisation can be run locally and viewed in a desktop browser using 

`panel serve [notebookname].ipynb --args foo=23 bar=./readme.paf`


