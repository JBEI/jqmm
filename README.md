# JBEI Quantitative Metabolic Modeling (jQMM) Library

The jQMM library is an **open-source**, **python**-based library for **bioengineering**. The jQMM library provides tools for **modeling internal metabolic fluxes** and **leveraging** other **'omics data** for metabolic engineering for the production of biofuels and other purposes. 

You can find the paper describing the jQMM library here:

https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1615-y

The jQMM library is intended to facilitate the design and metabolic engineering of organisms for biofuels and other chemicals, as well as investigations of cellular metabolism and ’omics data leveraging.  As an open source software project, JBEI hopes it will attract additions from the community and grow with the rapidly changing field of metabolic flux analysis.

Specifically, the jQMM library provides tools for the following capabilities:

* The jQMM is a complete toolbox for performing types of flux analysis that are typically disjoint: **Flux Balance Analysis** and **13C Metabolic Flux Analysis**.
* The jQMM can use 13C labeling experimental data to constrain **comprehensive genome-scale models** through a technique called two-scale 13C Metabolic Flux Analysis (2S-13C MFA).
* The jQMM demonstrates a method that uses **proteomics data** to produce actionable items that **increase biofuel production**.

All code is illustrated in a variety of **Jupyter notebook demos** (see /jupyterNotebooks folder) that enhance reproducibility and provide the capability to be modified by the user.

Code can be found in the /code folder whereas data used for the demonstrations can be found in the /data folder.

This library requires **GAMS** (https://www.gams.com/), **CPLEX** (https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/) and **CONOPT** (http://www.conopt.com/) licenses to work. Jupyter notebooks have been created with GAMS Release: 24.5.4 r54492 LEX-LEG x86 64bit/Linux and CONOPT version 3.17A .

## Getting Started

Clone this github repository and test the library by running the jupyter notebooks (http://jupyter.org/). To check that everything is installed properly, run notebook B1, or the "Quick tests" in notebook A0. Then choose tests from the other jupyter notebooks as desired. 

For ease of installation, a **docker version** of this library can be downloaded at:
https://github.com/JBEI/jqmm/tree/master/dockerContainer. 

Library dependencies can be found in: 
https://github.com/JBEI/jqmm/blob/master/dockerContainer/requirements.txt

## Computational requirements
For optimal performance in solving large genome-scale models with jQMM we recommend a computer system with at least 32 cpu cores (64 bit), 32 GB or more of ram, and 32 GB of solid state disk storage. The absolute minimum requirements are 2 cpu cores (64 bit), 4 GB of ram, and 16 GB of disk storage.


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
