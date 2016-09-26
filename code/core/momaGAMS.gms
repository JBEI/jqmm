$offempty

$offlisting
$offsymlist
$offuellist
$offinclude
$offsymxref
$offuelxref
$onempty

option profile=1;
OPTION decimals = 8
       sysout = off
       solprint = on
       reslim = 1000000
       iterlim = 10000000
       domlim = 10
       limcol = 1000
       limrow = 1000
       optca = 0.0
       optcr = 0.0
       work = 10000000
       nlp = conopt 
       mip = cplex
       solprint = on;


***********Genome-scale (GEM) SETS&PARS***************

sets
	iGEM    set all metabolites for GEM
$include "hybridmets.txt"
	jGEM    set all reactions for GEM
$include "hybridreacs.txt"
;

parameters
Sbig(iGEM,jGEM)  Stochiometry matrix for GEM
bvec(iGEM)       b vector for FBA problem (metabolite accumulation)
lbvec(jGEM)      lower bound for flux vector
ubvec(jGEM)      upper bound for flux vector
vGEMini(jGEM)    Target fluxes for MOMA
;

$include "Sbig.txt"
$include "bvec.txt"
$include "lbvec.txt"
$include "ubvec.txt"
$include "vGEMini.txt"






***********VARIABLES***************

* HYBRID
Free Variables 
	OFF                       Objective function for MOMA procedure
;

Positive Variables
	VGEM(jGEM)                GEM fluxes for large problem
;








***********EQUATIONS**************

* INIT
Equations
	STOICHIOMETRYFBA(iGEM)      Stoichiometric flux conservation constraint for GEM
	OBJFUN                      Objective function for initialization (as close as possible to Vini meeting stoich constraints)
;

STOICHIOMETRYFBA(iGEM)                                        .. bvec(iGEM)  =e= sum((jGEM),Sbig(iGEM,jGEM)*VGEM(jGEM)); 
OBJFUN                                                        .. OFF =e= sum(jGEM, sqr(VGEM(jGEM)-vGEMini(jGEM)) );





***********UPPER AND LOWER CONSTRAINTS*********

VGEM.lo(jGEM)          = lbvec(jGEM);
VGEM.up(jGEM)          = ubvec(jGEM);

***********INITIAL CONDITIONS*************** 
VGEM.l(jGEM)           = vGEMini(jGEM);






***********MODEL DECLARATION****************

model MOMA /
STOICHIOMETRYFBA
OBJFUN
/;
MOMA.OPTFILE = 1


	
** Solving MOMA
Solve MOMA using nlp minimizing OFF;

***********OUTPUT****************************


Parameters
        solvestat                       Exit status of solver
        modelstat                       Exit status of model
        numequ                          Number of equations in model 
        numnz                           Number of nonzero elements in model
        numvar                          Number of variables in model
        resusd                          Seconds used to reach final model state
;

* General info output
solvestat                   = MOMA.solvestat;
modelstat                   = MOMA.modelstat;
numequ                      = MOMA.numequ;
numnz                       = MOMA.numnz;
numvar                      = MOMA.numvar;
resusd                      = MOMA.resusd;

* Objective function and general info
file OFGenInfo/"OFGenInfo.txt"/;
put OFGenInfo
put "OF = ",OFF.l:10:5," ";     put /;   
put "solvestat = ",solvestat," "; put /;
put "modelstat = ",modelstat," "; put /;   
put "numequ    = ",numequ," ";    put /;
put "numnz     = ",numnz," ";     put /;
put "numvar    = ",numvar," ";    put /;
put "resusd    = ",resusd," ";    put /;    
putclose OFGenInfo;

* FBA Fluxes output
file VGEMoutfile/"MOMAout.txt"/;
put VGEMoutfile
loop(jGEM,
	put "vGEM('",jGEM.tl:30:10,"') = ",VGEM.l(jGEM):20:15," "; 	put /;	
     );
putclose VGEMoutfile;