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







***********HYBRID SETS&PARS***************

sets
	iFBA    set all metabolites for large FBA problem
$include "hybridmets.txt"
	jFBA    set all reactions for large FBA problem
$include "hybridreacs.txt"
        fluxes(jFBA)  reactions for which FVA is performed
$include "studiedfluxes.txt"	
;



parameters
	Sbig(iFBA,jFBA)  Stochiometry matrix for large FBA problem
	bvec(iFBA)       b vector for FBA problem (metabolite accumulation)
	lbvec(jFBA)      lower bound for flux vector
	ubvec(jFBA)      upper bound for flux vector
	cvec(jFBA)       optimization vector for FBA
	mvec(jFBA)       optimization vector for FVA
;

$include "Sbig.txt"
$include "bvec.txt"
$include "lbvec.txt"
$include "ubvec.txt"
$include "cvec.txt"

parameters
	store_maxs(jFBA)  stores maximum value
	store_mins(jFBA)  sotres minimum value


***********VARIABLES***************
* HYBRID
Free Variables 
	OFFBA                     Initialization objective function
	OFFVA
	OFFBAgrowth               Max growth objective function
;

Positive Variables
	VFBA(jFBA)                FBA fluxes for large problem
;

***********EQUATIONS**************

* HYBRID
Equations
	STOICHIOMETRYFBA(iFBA)         Stoichiometric flux conservation constraint for large problem
	OBJFUNFBA                      Standard FBA objective function for large problem
	OBJFUNFVA                      FVA objective function (changes in order to minimize and maximize each flux)
;

STOICHIOMETRYFBA(iFBA)             ..         bvec(iFBA)  =e= sum((jFBA),Sbig(iFBA,jFBA)*VFBA(jFBA)); 
OBJFUNFBA                          ..         OFFBA       =e= sum((jFBA),cvec(jFBA)*VFBA(jFBA));
OBJFUNFVA                          ..         OFFVA       =e= sum((jFBA),mvec(jFBA)*VFBA(jFBA));

***********UPPER AND LOWER CONSTRAINTS*********

VFBA.lo(jFBA)          = lbvec(jFBA);
VFBA.up(jFBA)          = ubvec(jFBA);


***********INITIAL CONDITIONS*************** 
VFBA.l(jFBA)           = uniform(0,1);

***********MODEL DECLARATION****************

model INITFBA /
STOICHIOMETRYFBA
OBJFUNFBA
/;
INITFBA.OPTFILE = 1

model FVA /
STOICHIOMETRYFBA
OBJFUNFBA
OBJFUNFVA
/;
FVA.OPTFILE = 1

***********MODEL SOLUTION*********************

** FBA Problem
Solve INITFBA using lp maximizing OFFBA;

* FBA Fluxes output
file VFBAoutfile/"VFBAout.txt"/;
put VFBAoutfile
loop(jFBA,
	put "vFBA('",jFBA.tl:30:10,"') = ",VFBA.l(jFBA):20:15," "; 	put /;	
     );
putclose VFBAoutfile;

** FVA problem
*VFBA.fx('BiomassEcoli') = OFFBA.l;
*OFFBA.fx = OFFBA.l;
* The 0.95 should be a parameter
OFFBA.lo = 0.95*OFFBA.l;  
OFFBA.up = OFFBA.l;
mvec(jFBA) = 0;
*alias(jFBA,fluxes)
* Starting jobs	
loop(fluxes,mvec(fluxes)=1;
	Solve FVA using lp maximizing OFFVA;
	store_maxs(fluxes) = OFFVA.l;
	Solve FVA using lp minimizing OFFVA;
	store_mins(fluxes) = OFFVA.l;
	mvec(fluxes) = 0;
);

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
solvestat                   = INITFBA.solvestat;
modelstat                   = INITFBA.modelstat;
numequ                      = INITFBA.numequ;
numnz                       = INITFBA.numnz;
numvar                      = INITFBA.numvar;
resusd                      = INITFBA.resusd;

    
* Objective function and general info
file OFGenInfo/"OFGenInfo.txt"/;
put OFGenInfo
put "OF = ",OFFBA.l:10:5," ";     put /;   
put "solvestat = ",solvestat," "; put /;
put "modelstat = ",modelstat," "; put /;   
put "numequ    = ",numequ," ";    put /;
put "numnz     = ",numnz," ";     put /;
put "numvar    = ",numvar," ";    put /;
put "resusd    = ",resusd," ";    put /;    
putclose OFGenInfo;




* Max and mins output
file Vmaxoutfile/"Vmaxout.txt"/;
put Vmaxoutfile
loop(fluxes,
	put "vmax('",fluxes.tl:30:10,"') = ",store_maxs(fluxes):20:15," "; 	put /;	
     );
putclose Vmaxoutfile;

file Vminoutfile/"Vminout.txt"/;
put Vminoutfile
loop(fluxes,
	put "vmin('",fluxes.tl:30:10,"') = ",store_mins(fluxes):20:15," "; 	put /;	
     );
putclose Vminoutfile;