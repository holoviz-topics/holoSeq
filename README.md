# hv-notebooks
Jupyter visualisation notebooks supporting millions of points smoothly - must run in Jupyter for datashader

Run 

> ! pip install datashader dask[dataframe] holoviews[recommended] pandas matplotlib bokeh


to load all the dependencies, then run the main script. 
Both take about 90 seconds each. 

There are 3 plots showing H1/H1, H2/H2 "cis" pairs, and H1/H2 "trans" pairs. Click anywhere to see the coordinates or drag, and use the mouse wheel to zoom. When run locally with 4.2.5, it is astounding. 14 million points. The old version does not do datashader properly so goes all blocky as you zoom.

There's ~200 lines of code to make the axes - must be ordinal for datashader to work so have to map the contigs end to end before can assign x/y grid coordinates to the pairs read from the paf input file.

Contigs and sizes have to be inferred from the data so it's 2 passes. The actual holoviews/panel code is ~20 lines once the data are in a grid. This stuff is dynamite.

