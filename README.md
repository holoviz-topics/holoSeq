# holoSeq

#### An interactive visualisation framework for genomic data, scaling to tens of millions of points.
*Built using the Holoviews and IPython notebook ecosystem.*

<img src="https://github.com/fubar2/hv-notebooks/blob/main/h2.gif" alt="zoom demo" width="125"/>       <img src="https://github.com/fubar2/hv-notebooks/blob/main/h1.gif" alt="zoom demo" width="125"/>

A detailed framework [specification is here.](https://github.com/fubar2/holoSeq/blob/main/HoloSeqOverview.md). Briefly, the framework relies an optionally gzipped text data format containing
pre-computed plot coordinates, and a header describing the reference sequence or sequences that those coordinates refer to. Converters for common genomic annotation formats, such as PAF, bigwig and bed
will be supplied. Any number of input files in that data format can be supplied to a generic IPython notebook, for organisation and interactive display.

Scaling and zooming rely on datashader running on a Holoviews or Bokeh server. Static images can be captured. Interactive HTML can be exported, but without a datashader provider, zoomed detail is 
very limited.

As proof of concept, a genomic HiC contact pair viewer that can show the entire map of 14 million pairs and any level of
zoom down to individual contact pair points is provided. Screenshots below use Arima HiC reads from the VGP mUroPar1 Arctic Ground Squirrel, 
processed through a Galaxy workflow.

For all these notebooks, first, run 

> ! pip install datashader dask[dataframe] holoviews[recommended] pandas matplotlib bokeh

to load all the dependencies, then run the main script. 
Both take about 90 seconds each. 

#### 1. hapsHiCpafHoloview.ipynb is a plotter for paired HiC contact pointss.ss, input as a PAF.

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

   
