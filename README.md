# hv-notebooks

### Jupyter visualisation notebooks supporting millions of points smoothly

This collection of notebooks use datashader so must run in Jupyter.
Static png or html can be exported but they are not interactive.

For all these notebooks, first, run 

> ! pip install datashader dask[dataframe] holoviews[recommended] pandas matplotlib bokeh

to load all the dependencies, then run the main script. 
Both take about 90 seconds each. 

#### 1. hapsHiCpafHoloview.ipynb is a plotter for paired HiC contact points, input as a PAF.

This is from the Arctic Ground Squirrel mUroPar1 HiC VGP data
<img src="https://github.com/user-attachments/assets/3f94291c-8905-40d4-aa5a-ba379812d67b" alt="h1-h1, h2-h2 pairs" width="800"/>

The strangely disordered block of HiC contact points on chrY looks very unusual.
Interestingly it corresponds to a region of non B DNA repeats in the T2T HG002 Y chromosome described in https://www.nature.com/articles/s41586-023-06457-y

![image](https://github.com/user-attachments/assets/606a1b3e-be5e-4915-a335-7c300601392e)


The notebook makes 3 interactive plots showing H1/H1, H2/H2 "cis" pairs, and H1/H2 "trans" pairs. 
These work well with 14 million points. Occasional slow patches but can pan and wheel zoom down from all ponts to single pairs easily. 

Click anywhere to see the coordinates or drag, and use the mouse wheel to zoom. When run locally with 4.2.5, it is astounding. 14 million points. The old version does not do datashader properly so goes all blocky as you zoom.

There's ~200 lines of code to make the axes - must be ordinal for datashader to work so have to map the contigs end to end before can assign x/y grid coordinates to the pairs read from the paf input file.

Contigs and sizes have to be inferred from the data so it's 2 passes. The actual holoviews/panel code is ~20 lines once the data are in a grid. 
Holoviews is dynamite.

Here are the first few pairs on chr1 for H1 and H2 - they are about 20k from the start of the assembly - presumably no HiC pairs in the actual telomeres.

<img src="https://github.com/user-attachments/assets/4d7a62b9-a98b-4854-a737-f8d68b0f8c4b" alt="h1-h1 and h2-h2 pairs" width="800"/>

and the H1-H2 pairs

<img src="https://github.com/user-attachments/assets/84f0066c-ce9e-4f14-92f9-4a8f019f2f9d" alt="h1-h2 pairs" width="500"/>


