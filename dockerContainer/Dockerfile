FROM debian:8 

MAINTAINER Tyler W. H. Backman "tbackman@lbl.gov"

USER root

ENV LANG C.UTF-8
WORKDIR /root

# Install debian packages
RUN apt-get update \
	&& apt-get install -y \
		screen \
		wget \
		vim \
		man \
		python-tk \
		pkg-config \
 		ipython \
		python-pip \ 
		python-dev \
		cython \
		python-tk \
		python-libxml2 \
		python-lxml \
		tcsh \
		git \
	&& apt-get autoremove -y \
	&& rm -rf /var/lib/apt/lists/*

# install pip packages
RUN pip install -U pip setuptools
ADD ./requirements.txt /root/
RUN pip install -r requirements.txt

# install script to download GAMS
ADD ./install_gams.sh /root/
ADD ./gamsbatch /root/
ENV PATH="${PATH}:/opt/gams/gams24.5_linux_x64_64_sfx"

# install jQMM
RUN git clone https://github.com/JBEI/jqmm.git \
    && mkdir -p /users/hgmartin/libraries/ \ 
    && mkdir -p /scratch/hgmartin_scratch/tests \ 
    && ln -s /root/jqmm /users/hgmartin/libraries/quantmodel

# Add Tini
ENV TINI_VERSION v0.14.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini
ENTRYPOINT ["/tini", "--"]

# run jupyter server
CMD ["jupyter", "notebook", "--no-browser", "--ip=*", "--notebook-dir=/root/jqmm"]
