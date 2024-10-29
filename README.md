# holoSeq

#### An interactive visualisation framework for genomic data, scaling to tens of millions of points.
*Built using the [Holoviews](https://holoviews.org/) and IPython notebook ecosystem.*

<img src="https://github.com/fubar2/hv-notebooks/blob/main/h2.gif" alt="zoom demo" width="125"/>       <img src="https://github.com/fubar2/hv-notebooks/blob/main/h1.gif" alt="zoom demo" width="125"/>

A draft framework [specification is here.](https://github.com/fubar2/holoSeq/blob/main/HoloSeqOverview.md). 

#### Quick start: Install and test the local demonstration

From the directory where you have cloned this repository, prepare the python virtual environment needed, then run the panel application.

```
python -m venv venv
. venv/bin/activate
pip install -r requirements.txt 
panel serve holoseq_display.py --args --inFile mUroPar1_cis1.hseq.gz --size 1000
```

Expect to see:

```
(venv) ross@pn50:~/rossgit/holoSeq$ panel serve holoseq_display.py --args --inFile mUroPar1_cis1.hseq.gz --size 1000
2024-10-29 17:06:33,645 Starting Bokeh server version 3.6.0 (running on Tornado 6.4.1)
2024-10-29 17:06:33,646 User authentication hooks NOT provided (default user enabled)
2024-10-29 17:06:33,649 Bokeh app running at: http://localhost:5006/holoseq_display
2024-10-29 17:06:33,649 Starting Bokeh server with process id: 388673
```
No processing will take place until a browser window is opened at the address shown, and it will take 10-20 seconds to read the 3.4M pairs and to show the interactive visualisation.
When it appears:
- Mouse click anywhere on the plot to see the coordinates.
- Zoom with the mouse scroll wheel
- Pan by grabbing with the left mouse button.
- Usual `Bokeh` display controls are available on the sidebar.
- Only pairs involving H1 contigs (H1 cis) are used in the demonstration.

Briefly, the framework is built around an optionally gzipped text data format, that supplies pre-computed plot coordinates, with a header 
describing the reference sequence or sequences forming the plot axes. Converters for common genomic annotation formats, such as PAF, bigwig and bed
will be provided. 

This design isolates the complexities of displaying many different kinds of annotation at genomic scale, from the details of converting a zoo of data
from external sources in important standard formats. Any number of input files can be supplied to the generic IPython notebook, where
they are automatically organised and displayed. Scaling and zooming rely on datashader running on a Holoviews or Bokeh server. 
Static images can be captured. Interactive HTML can be exported, but without a datashader provider, zoomed detail is limited. 

As proof of concept, a genomic HiC contact pair viewer that can show the entire map of 14 million pairs and any level of
zoom down to individual contact pair points is provided. Screenshots below use Arima HiC reads from the VGP *mUroPar1* Arctic Ground Squirrel, 
processed through a Galaxy workflow.

#### Proof of concept holoSeq data format disk sizes

Data size on disk using the format in the specification gives more than an order of magnitude, from about 300MB
of input PAF to 23MB of hseq format data.

Open a browser window at the address shown - default is  http://localhost:5006/holoseq_display
It will take 20 seconds or so to prepare and show the interactive visualisation.
The original PAF is 1.2GB and contains about 14M rows.
About 3.6M pairs of points had both contigs on the mUroPar1 paternal haplotype - so about
1/4 of all rows. The sample used in the demonstration is a 23M gzip containing all the information needed to plot 
these 3.6M pairs.

#### Prepare a PAF file containing the points to be plotted

The demonstration compressed holoSeq plot information was prepared with a PAF file prepared in a Galaxy VGP workflow, mapping the arima HiC reads
against the assembled haplotypes, running Bellerophon to remove chimeric Arima reads into a bam of pairs. These are converted to a sam file with header
using samtools view, and processed into a PAF with an AWK script

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
A conversion utility is provided for PAF inputs. The demonstration data were prepared using:

`python holoSeq_prepare_paf.py --inFile mUroPar1.paf --title "VGP Arctic Ground Squirrel arima HiC contact matrix, paternal haplotype" `

This step produces outputs containing subsets of contact point pairs such as

`mUroPar1.paf_cis1.hseq.gz`

These can be viewed like the supplied local demonstration example using panel.

### 1. hapsHiCpafHoloview.ipynb is a plotter for paired HiC contact pointss.ss, input as a PAF.

This is from the Arctic Ground Squirrel mUroPar1 HiC VGP data. It shows the three images and about 14 million points.

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

#### Gallery of screenshots

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

   
