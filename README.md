# JBEI Quantitative Metabolic Modeling (jQMM) Library

The jQMM library is an open-source, python-based library for bioengineering. The jQMM library provides tools for modeling internal metabolic fluxes and leveraging other 'omics data for metabolic engineering for the production of biofuels and other purposes. It is intended to facilitate the design and metabolic engineering of organisms for biofuels and other chemicals, as well as investigations of cellular metabolism and ’omics data leveraging.  As an open source software project, JBEI hopes it will attract additions from the community and grow with the rapidly changing field of metabolic flux analysis.

Specifically, the jQMM library provides tools for the following capabilities:

* The jQMM is a complete toolbox for performing types of flux analysis that are typically disjoint: Flux Balance Analysis and 13C Metabolic Flux Analysis.
* The jQMM can use 13C labeling experimental data to constrain comprehensive genome-scale models through a technique called two-scale 13C Metabolic Flux Analysis (2S-13C MFA).
* The jQMM demonstrates a method that uses proteomics data to produce actionable items that increase biofuel production.

All code is illustrated in a variety of Jupyter notebook demos (see /jupyterNotebooks folder) that enhance reproducibility and provide the capability to be modified by the user. A paper explaining the different capabilities is now in review, and a link to it will be posted as soon as it is published.

Code can be found in the /code folder whereas data used for the demonstratiosn can be found in the /data folder.

This library requires GAMS (https://www.gams.com/) and CONOPT (http://www.conopt.com/) licenses to work. Jupyter notebooks have been created with GAMS Release: 24.5.4 r54492 LEX-LEG x86 64bit/Linux and CONOPT version 3.17A .

## Getting Started

Clone this github repository and test the library by running the jupyter notebooks (http://jupyter.org/). Run first the "Quick tests" in notebook A0, and then choose tests from the other jupyter notebooks as desired.  

We also provide a Docker version of this library which automates the installation, including all required software prerequsites. To use the Docker version, first download and install Docker from https://www.docker.com, and then follow the jQMM docker instructions here: https://hub.docker.com/r/mhgarci1/jqmm/

## License

JBEI Quantitative Metabolic Modeling library (JQMM),
Copyright (c) 2016, The Regents of the University of California,
through Lawrence Berkeley National Laboratory (subject to receipt of
any required approvals from the U.S. Dept. of Energy).  All rights
reserved.

Redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following
conditions are met:

(1) Redistributions of source code must retain the above
copyright notice, this list of conditions and the following
disclaimer.

(2) Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials provided
with the distribution.

(3) Neither the name of the University of California, Lawrence
Berkeley National Laboratory, U.S. Dept. of Energy nor the names
of its contributors may be used to endorse or promote products
derived from this software without specific prior written
permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
