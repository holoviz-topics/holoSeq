# HoloSeq

This repository contains an application that allows large scale genomic data with tens of millions of features to be viewed interactively,
as 1D charts or 2D heatmaps, using generic visualisation infrastructure built using the [Holoviews ecosystem](https://holoviews.org/). 

The presentation layer supports genome scale feature data associated with a genomic reference or other sequence. The user can pan and zoom 
smoothly from whole genomes down to individual points, in a web browser running on a suitable laptop or in Galaxy. At present, 2D plots can be
made from HiC contact pairs as PAF, and from mashmap sequence similarity as PAF. Bigwig and gff format tracks can be converted into coordinates ready for 1D plotting.

The intermediate file containing pre-computed coordinates and metadata for reconstructing any converted track, can be saved as a compressed `hseq.gz` 
file and re-displayed, saving the time taken to compute potentially millions of coordinates. This separation was designed for use in a Galaxy interactive
tool, where it is easy to save all the intermediate coordinate files as a collection for re-display without recalculating all the millions of coordinate pairs.

The main use case envisioned is a central repository of pre-computed plots to make annotation of the VGP genomic data easily accessible.
Precomputed plot tracks can be reused indefinitely, mixed and matched by the user to suit their needs.

Obviously, each species reference genome has a different coordinate system since they have different genomes.
It's important to emphasise that the two haplotypes being assembled from each species also have a different coordinate system because 
contig lengths nearly always differ. As a result plots cannot share a reference sequence axis unless the data were all mapped to exactly the same set of contigs.
Different haplotypes can be displayed with tracks made from that haplotype reference sequence side by side or stacked, but otherwise, linking will always get 
out of sync. Linking is automatic between plots with identical axes in Holoviews, so the metadata must always identify what was used to map any track so it
can be displayed with other tracks from the same sequence. 

Futzing to try to make things line up will almost certainly fail at some point - let's not bother. It might be possible to make two different reference sequence tracks display the same 
chr:start-end region if they share contig names, but even then, there will be tricky differences - important to respect the original mapping for reliable synchronised displays.

### Coordinate system for genomic assembly data

Genomic reference data forms the backbone for any annotation browser, and each genomic region or position must be translated from something like CHR3:10000 into an internally
consistent axis coordinate system for reliable display on a holoviews datamap. The first time a new species genome is assembled from sequence data, it's like taking all the pieces from a woodchipper and gluing
them back together into the pre-existing tree, but this tree has never been seen before, making the task even harder. New genomes are painstakingly assembled into multiple scaffolds 
and super-scaffolds as curation proceeds, referred to generally as "contigs". These are eventually placed into a parsimonious set of chromosomes, in the final "reference" genome. A reference genome such as the 
current Human [release](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/) is the final product of many iterations of assemblies and polishing.

Whatever the stage of assembly and polishing, genomes are typically handled in fasta format. A newly assembled haplotype may have thousands of contigs that have not yet been merged into chromosomes.

For an interactive genome browser, tracks typically run horizontally, from the start of the reference sequence on the left to the last nucleotide of the last contig.

Holoviews dynamic maps require ordinal axis coordinates. Contigs must be ordered on the axes, usually by name or by length,
so they can be mapped as axes, with tick marks and labels.

Features must have a position to locate them in the reference sequence used to create the axis. 
Position is usually described by the name of the contig, and the number of bases from the start of the contig to the start of the feature. 
Some features have a length while others for 2D grids are points. Many features may have optional annotation.

To convert feature positions into plot axis coordinates, contig lengths are cumulated in the order given, and a zero is inserted,
to give the ordinal position on the axis, for the first base of each ordered contig. The ordering makes a big difference to the visualisation.
Typically, contigs are length sorted, but sometimes it makes more sense to use the names allocated at assembly - these depend on the assembly team
and for the VGP, are typically SUPER for large, chromosome like features and SCAFFOLDS for less complete pieces. Then there are UNLOC versions
of both of these where the actual ordering is not yet known well enough to assemble into a chromosome.

When a feature is mapped, the appropriate contig's cumulated start is added to the feature offset, to give the 
ordinal start coordinate on that axis for the start of the feature. 

A linear model makes sense for a contig, but the idea of making a single linear axis from many contigs is a fiction - they exist as multiple chromosomes, each with a
maternal and paternal haplotype strand. It is important to understand that this fiction is used here to make it easy to "see" all the data. It is a model,
and like all models, wrong, but in this case, useful if the limitations are understood. Sequential contigs on any axis are not sequential in the nucleus where the entire
genome is rolled up into a complicated kind of hairball, confined to a volume of about 10 cubic microns. 

## Metadata and data format for 1D and 2D features on pre-mapped axis coordinates

The converters produce gzip compressed text files, ready for the display application.

The text file must start with a header section, where every row begins with `@`.

The first row of the must be either `@v1HoloSeq1D [bar|scatter|line]` or `@v1HoloSeq2D`, or the data will not be processed.

Rows starting with '@@' are metadata, such as the plot title, and the URI and name of the reference sequence or haplotype.

For 1D data, the chart type may be one of `bar`, `scatter` or `line`. Default is `bar`. Regions with 4 or more SD above or below the global mean are 
emphasised

2D data will be presented as an autoscaling density heatmap. The header and data might only have 1 axis name, for example where HiC pairs from one haplotype are plotted
with that sequence on both axes, or 2 axis names, if HiC pairs involving both haplotypes, one on each axis, are being plotted.

The subsequent header rows must have the plot title, plot type, axis names, contig names and their cumulated lengths, delimited by whitespace, and starting with `@` such as

```
@v1HoloSeq1D bar
@@title a very small bar plot
@@refURI https://foo.bar.zot.fasta
@H1 chr1 0
@H1 chr2 2500
@H1 chr3 7000
@H1 chr4 9500
```

In this case, H1 will show four chromosomes starting at each of the positions shown.

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

`panel serve holoseq_display.py --args --inFile foo.gz bar.gz baz.gz zot.gz --title my holoseq plot`

An interactive tool for Galaxy will be created.


