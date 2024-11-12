# holoSeq: Interactive viewing of genomic annotation

### Interactively browse millions of features in a laptop browser window 
*Browse 1D charts and 2D heatmap plots, that scale themselves when zoomed from the entire genome down to individual points and back.* 
*Designed for Galaxy interactive tools. Works on a laptop. Built on the [Holoviews](https://holoviews.org/) and IPython notebook ecosystem.*

<img src="https://github.com/fubar2/hv-notebooks/blob/main/h2.gif" alt="zoom demo" width="125"/>       <img src="https://github.com/fubar2/hv-notebooks/blob/main/h1.gif" alt="zoom demo" width="125"/>

## Project status

*Help wanted. PR and suggestions welcomed.*

This is new work in progress. 

Development started in late October 2024. 
A draft framework [specification is here.](https://github.com/fubar2/holoSeq/blob/main/HoloSeqOverview.md). 

## Core idea: Features on intervals arranged along linear axes for browsing

This proof of concept runs in a notebook, or if the dependencies are available, can be served from
this repository's root, as:
> panel serve holoSeq_random.py --show

Edit the default 10000 xmax value to get a sense of scale capacity - 10M is not a problem.
Most of the code is needed to create some sample contigs of fixed length into some arbitrary ordering along the axes.
Their lengths are cumulated and the resulting array represents axis offsets to the starting nucleotide of each contig,
used for calculating the x and y coordinates for randomly generated points so the are correctly located
in an internally consistent 2D space. 

This is proposed as a generalisable model for a linear display of genomic features located on a set of independent contigs in
holoSeq displays.

```
# see https://github.com/fubar2/holoSeq
# Needs dependencies - uncomment and run the next line to install them in a notebook
# ! pip install datashader dask[dataframe] holoviews[recommended] pandas matplotlib bokeh
#

from bisect import bisect_left
from collections import OrderedDict
import itertools
import numpy as np

import holoviews as hv
import panel as pn

from holoviews.operation.datashader import (
    rasterize,
    dynspread,
)

hv.extension("bokeh")
pn.extension()

def showTap(x, y):
    # populate a string widget with tap coordinates translated into contig and offset
    if np.isnan(x) or np.isnan(y):
        s = "Click on plot for coordinates as contig:offset"
    else:
        i = bisect_left(hstarts, x)
        chrx = hnames[i - 1]
        offsx = x - hstarts[i - 1]
        i = bisect_left(hstarts, y)
        chry = hnames[i - 1]
        offsy = y - hstarts[i - 1]
        s = "X axis %s:%d Y axis %s:%d" % (chrx, offsx, chry, offsy)
    str_pane = pn.pane.Str(
        s,
        styles={"font-size": "10pt", "color": "darkblue", "text-align": "center"},
        width=width,
    )
    return str_pane

xwidth = 3000000
xmax = 10000
width = 1000
rng = np.random.default_rng(1) # all plots will be identical !
hlen = [xwidth/2, xwidth/4, xwidth/8, xwidth/8]
hstarts = list(itertools.accumulate(hlen))
hstarts.insert(0,0)
hnames = ['chr%d' % i for i in range(1,5)]
xcoords = np.array([rng.uniform(0,xwidth) for i in range(xmax)])
yoffset = rng.normal(0,100000,xmax)
ycoords = np.array(xcoords + yoffset)
title = "%d random points" % xmax
ticks = [(hstarts[i], hnames[i]) for i in range(4)]
points = hv.Points((xcoords,ycoords))
stream = hv.streams.Tap(source=points, x=np.nan, y=np.nan)
showloc = pn.bind(showTap, x=stream.param.x, y=stream.param.y)
pnc = pn.Column(
    showloc,
    pn.pane.HoloViews(
        dynspread(
            rasterize(points)
            .relabel("%s" % title)
            .opts(
                cmap="inferno",
                cnorm="log",
                colorbar=True,
                width=width,
                height=width,
                xticks=ticks,
                yticks=ticks,
                xrotation=45,
                fontsize={"xticks": 7, "yticks": 7},
                tools=["tap"],
                scalebar=True,
                scalebar_range="x",
                scalebar_location="top_left",
                scalebar_unit=("bp"),
            )
        )
    ),
)


pnc.servable(title=title, )
```

## Demonstration with millions of pairs from HiC data

Requires both Python and the python-venv module installed. 
These steps ensure that the virtual environment directory can be made and deleted without any other changes to your system.

1. Clone this repository in a convenient and disposable directory:

`git clone https://github.com/fubar2/holoSeq`

2. Change directory to that new holoSeq subdirectory, prepare the python virtual environment needed, and start the demonstration

```
cd holoSeq
python -m venv venv
. venv/bin/activate
pip install -r requirements.txt 
panel serve holoseq_display.py --show --args --inFile mUroPar1.paf_cis1.hseq.gz --size 1000

```

Expect to see output shown below, and a web browser window should pop open. Takes 10-20 seconds to read the 3.4M pairs and to show the interactive visualisation:

```
(venv) ross@pn50:~/rossgit/holoSeq$ panel serve holoseq_display.py --args --inFile mUroPar1_cis1.hseq.gz --size 1000
2024-10-29 17:06:33,645 Starting Bokeh server version 3.6.0 (running on Tornado 6.4.1)
2024-10-29 17:06:33,646 User authentication hooks NOT provided (default user enabled)
2024-10-29 17:06:33,649 Bokeh app running at: http://localhost:5006/holoseq_display
2024-10-29 17:06:33,649 Starting Bokeh server with process id: 388673
```

- `--show` option opens a browser window and immediately starts processing the script.
- `--size` option determines both x and y plot window size.
- `--inFile` must be the path of a prepared hseq compressed coordinate data file, such as the supplied demonstration HiC heatmap.
- Delay will vary with processing power and the number of coordinates

When the plot appears, the view is controlled by the usual Bokeh tools in the toolbar on the left of the plot.

- Mouse click anywhere on the plot to see the coordinates.
- Zoom with the mouse scroll wheel
- Pan by grabbing with the left mouse button.
- Only pairs involving H1 contigs (H1 cis) are used in the demonstration.

Briefly, the framework creates the [minimum data required](https://github.com/fubar2/holoSeq/blob/main/HoloSeqOverview.md) to create a plot.
A genome lengths file is required, and the named contigs can be reordered by name or length. The axes are defined by the ordering. The
lengths are cumulated to give an offset to the first nucleotide of each contig, so the track can be read and feature locations 
converted into the plot coordinate system, and stored as a compressed intermediate file. The display application reads these pre-computed plot 
coordinate files, with enough metadata about the reference sequence to add tic marks to the axes and to back-calculate the stream of user tap coordinates. 
A converter for PAF to compressed hseq format for input is available and was used to generate the demonstration. 
Bigwig is working and other common genomic annotation formats, such as gff and vcf will follow. 

Multiple input files will produce a stack of plots that work independently:

`panel serve holoseq_display.py --show --args --inFile mUroPar1.paf_cis1.hseq.gz small.paf_cis1.hseq.gz  --size 1000`


## holoSeq data format disk sizes

Data size on disk using the hseq format provides more than an order of magnitude shrinkage, from about 300MB
of input PAF to 23MB of hseq coordinate format data.

The original PAF is 1.2GB and contains about 14M rows. About 3.6M pairs of points had both contigs on the mUroPar1 paternal haplotype - so about
1/4 of all rows. The sample used in the demonstration is a 23M gzip containing all the information needed to plot 
these 3.6M pairs.

## Prepare a PAF file containing the points to be plotted

The demonstration holoSeq plot hseq format coordinate data was prepared from a PAF file. That was output from a Galaxy VGP workflow, mapping the arima HiC reads
against the assembled haplotypes, running Bellerophon to remove chimeric Arima reads and merging into a bam containing all HiC pairs. 
These were converted to a sam file with header using `samtools view`, then processed into a PAF using this AWK script:

```
#!/bin/awk -f
{ 
if ($1=="@SQ") { chromlen[substr($2,4)] = substr($3,4) }
else
if (substr($1,0,1) == "@") {foo = 1}
else
{
if ($2 AND 16) { sign1="-"};
 
if (! ($2 AND 16)) { sign1="+"};

if ($2 AND 32) { sign2="-"};

if (! ($2 AND 32)) { sign2="+"};

if ($7 == "=") 
{ printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $3, chromlen[$3], $4, $4+length($10), sign1, $3, chromlen[$3], $8, $8 + length($10), length($10), length($10), 255)}
else
{ printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $3, chromlen[$3], $4, $4+length($10), sign1, $7, chromlen[$7], $8, $8 + length($10), length($10), length($10), 255)}
}
}

```

This repository includes a python script conversion utility for PAF inputs, `holoSeq_prepare_paf.py`, that works with the awk PAF output and converts it into
a compressed coordinate file. The compressed demonstration plotting data were prepared using:

`python holoSeq_prepare_paf.py --inFile mUroPar1.paf --title "VGP Arctic Ground Squirrel arima HiC contact matrix, paternal haplotype" `

for Mashmap paf inputs, use: 

`python holoSeq_prepare_paf.py --inFile  hg002_2k99.paf --title "hg002 Mashmap" --hap_indicator None --contig_sort length`

This step produces outputs containing subsets of contact point pairs.

These can be viewed like the supplied local demonstration example using `panel` as described above

### Gallery of screenshots

The images are highly interactive in a Jupyter notebook but these static screenshots give some idea of what they can do.

1. ChrY 

Zooming in to the H1 Y chromosome shows

<img src="https://github.com/user-attachments/assets/3f94291c-8905-40d4-aa5a-ba379812d67b" alt="h1-h1, h2-h2 pairs" width="800"/>

There is an obvious yellowish diagonal line, where the majority of contact pairs have their two ends close together along the Y chromosome. 
There is a large, dense and chaotic rectangular pattern at one end, and some dense linear features at an angle to the main diagonal.
There is a variable, thin scattering of off-diagonal dots, representing the less common pairs with their two ends much further apart on the Y chromosome.

Perhaps there is some biology there because it seems to correspond to the very large region of non B DNA repeats in the Human T2T HG002 Y chromosome, described 
in https://www.nature.com/articles/s41586-023-06457-y and seen in this figure from the paper.

![image](https://github.com/user-attachments/assets/606a1b3e-be5e-4915-a335-7c300601392e)

2. Some validation of point coordinates
   
Here are the first few pairs on chr1 for H1 and H2 - they are about 20k from the start of the assembly - presumably no HiC pairs in the actual telomeres.

<img src="https://github.com/user-attachments/assets/4d7a62b9-a98b-4854-a737-f8d68b0f8c4b" alt="h1-h1 and h2-h2 pairs" width="800"/>

and the H1-H2 pairs

<img src="https://github.com/user-attachments/assets/84f0066c-ce9e-4f14-92f9-4a8f019f2f9d" alt="h1-h2 pairs" width="500"/>


The first records in the PAF for a pair nearest offset 0 are at 19087/19310 and 19310/19087 so expect to see the very first points in SUPER1 of H1 somewhere near those coordinates.
```
SUPER_1H1       284260672       19087   19238   -       SUPER_1H1       284260672       19310   19461   151     151     255
SUPER_1H1       284260672       19310   19461   -       SUPER_1H1       284260672       19087   19238   151     151     255


```
Zooming right in to the start of SUPER1 on H1 shows that the single point shown above resolves into a pair of points at the expected coordinates, confirming that the image may be showing what's in the data... 

3. Interesting features

Below are the points in H1 and H2 in the same region of the trans plot where I noticed an odd looking
break in the main diagonal as an example of what's possible in the notebook.

<img src="https://github.com/user-attachments/assets/20252a6d-2058-484d-a70d-413cc22cd706" alt="h1-h2 pairs" width="500"/>

### Historical material 

#### Original proof of concept application: hapsHiCpafHoloview.ipynb is a plotter for paired HiC contact pointss.ss, input as a PAF.

This uses HiC data from the [Arctic Ground Squirrel mUroPar1](). It shows the three images and about 14 million points.

<img src="https://github.com/user-attachments/assets/e288f295-84b0-4121-9099-db5007445a27" alt="h1-h1, h2-h2 pairs" width="800"/>


The matrix is symmetrical so only one half is needed but it is visually more impressive as shown

The notebook makes 3 interactive plots showing H1/H1, H2/H2 "cis" pairs, and H1/H2 "trans" pairs. 
These work well with 14 million points. Occasional pauses before rescaling and refresh, but can pan and wheel zoom down from all ponts to single pairs easily. 

Click anywhere to see the coordinates or drag, and use the mouse wheel to zoom. When run locally with 4.2.5, it is astounding. 14 million points. The old version does not do datashader properly so goes all blocky as you zoom.

There's ~250 lines of code to make the axes - must be ordinal for datashader to work so have to map the contigs end to end before can assign x/y grid coordinates to the pairs read from the paf input file.

That code already contains a PAF converter, and the generic IPython notebook visualiser will be built using the existing visualisation components, adding 1D converters and tracks, with 
grouping by common reference sequence in the header and vertical stack layout.

Contigs and sizes have to be inferred from the data so it's 2 passes. The actual holoviews/panel code is ~20 lines once the data are in a grid. 
Holoviews is dynamite.


   
