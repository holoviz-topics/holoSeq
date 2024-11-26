# Base Image
FROM ubuntu

# Maintainer
MAINTAINER Ross Lazarus <ross.lazarus@gmail.com>
USER root
RUN apt-get update \
&& apt-get install -y --no-install-recommends python3 python3-venv python3-pip python3-wheel wget unzip nano git \
&& mkdir -p /data \
&& useradd bio \
&& chown -R bio /data /home/bio
USER bio
RUN cd /data \
&& python3 -m venv venv \
&& . venv/bin/activate \
&& python3 -m pip install --upgrade pip numpy holoviews[recommended] datashader dask[dataframe] pandas bokeh 
WORKDIR /data/
