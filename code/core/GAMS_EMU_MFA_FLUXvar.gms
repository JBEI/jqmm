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
       reslim = 4000
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
	i     set of all metabolites
$include "mets.txt"
        exchangemet(i) set of exchange metabolites
$include "ExchaMets.txt"
        excludedmet(i)  excluded metabolites
$include "ExcMets.txt"
        sourcemet(i)  excluded metabolites
$include "SourceMets.txt"	
        j     set of all reactions
$include "rxns.txt"
        jmeas(j) set of measured reactions
$include "jmeas.txt"
        jsugg(j) set of suggested reactions
$include "jsugg.txt"
;


parameters
	S(i,j)        Stoichiometry matrix
	vUPPER(jmeas) upper boundary for fluxes / /
	vLOWER(jmeas) lower boundary for fluxes / /
	vSUGG(jsugg)  suggested fluxes / /
	Vini(j)       Initialization fluxes
	maxflow13CU    maximum flux value for 13C problem (upper value)	
;




$include "StoichNetwork.txt"
$include "vUPPER.txt";
$include "vLOWER.txt";
$include "vSUGG.txt"

***********EMU system SETS&PARS***************

sets
	e   set of all emus	
$include "all_emu.txt"
	dummy(e) dummy emus for condensation reacs
$include "dummy_emus.txt"
        frag   fragment names for MS data (e.g. Ala159)
$include "fragmentsMS.txt"
        fragFit(frag) fragment names used for the fit
$include "fitFrags.txt"
	m values of m for emu MDVs
$include "mset.txt"
	n values of n for MS MDVs
$include "nset.txt"
;

alias(i,ii);
alias(e,ee);
alias(m,mm);
alias(n,nn);

parameters
	ncarbons(e)        number of carbons for each emu
	ncarbonsf(frag)    number of carbons for each fragment
	labelexp           Table with MS labeling data
	labelstd           Table with MS labeling error
	labelnorm          Table with MS labeling normalization
        aacorr             Correspondance of aminoacids and emus
	gamma              gamma matrix for derivatization effects   
	EMM                EMU mapping matrix
	metemu(i,e)        1 if e is an emu of i and zero otherwise 
;
$include "emu_carbon_numbers.txt";
$include "frag_carbon_numbers.txt";
$include "GCMSout.txt";
$include "GCMSerr.txt";
$include "GCMSnorm.txt";
$include "aacorr.txt";
$include "gamma_mat.txt";
$include "emm_simplified.txt";
$include "met_emu.txt";

***********2-SCALE SETS&PARS***************

sets
	jGEM    set all reactions for GEM (level 1 description of fluxes)
$include "linkedreacs.txt"
;










parameters
	reacmap(jGEM,j)  reaction map between GEM reactions and EMU transitions
;

$include "reacmap.txt"


























***************EXTRA SETS&PARAMETERS******************

sets
	fluxFocus(jGEM)  flux focused on for variability analysis
$include "fluxFocus.txt"
        bounds bounds for labeling focused on /'min','max'/
;

Scalar labelCompNormMin;
$include "labelCompNormMin.txt"


* Random start parameters and others
Parameters
	randparam(j)
	fluxBounds(bounds)
	Vmax(j)    GEM flux values for maximum of fluxfocus
	Vmin(j)    GEM flux values for minimum of fluxfocus
	OFEMUMax         EMU objective function for max fluxfocus
	OFEMUMin         EMU objective function for max fluxfocus
	labelcompMax(frag,n)       Simulated MS data for maximum fluxfocus
	labelcompMin(frag,n)       Simulated MS data for mininum fluxfocus
;

Scalar randseed;
$include "randseed.txt"
execseed = randseed;


loop(j,
	randparam(j)= uniform(0,1)
);







***********VARIABLES***************
* INIT
Free Variables 
	OFFINIT                     Initialization objective function
;

Positive Variables
	V(j)                       GEM fluxes
;

* EMU
Free Variables
	OFEMU                      EMU objective function
	OFvar                      Flux objective function
	bal1
	bal2
;

Positive Variables
	f(e,m)                     EMU MDVs
	labelcomp(frag,n)          Simulated MS data
	labelcompnorm(frag)        Norm for simulated MS data
;

* 2-SCALE
Free Variables
	VGEM(jGEM)                GEM fluxes 
;






***********EQUATIONS**************

* INIT
Equations
	STOICHIOMETRY(i)            Stoichiometric flux conservation constraint
	OBJFUNINIT                  Objective function for initialization (as close as possible to Vini meeting stoich constraints)
	MAPPING(jGEM)                  Mapping between EMU and large FBA problem
;


STOICHIOMETRY(i)$(i(i) - exchangemet(i) - excludedmet(i))     .. 0 =e= sum(j,S(i,j)*V(j));
OBJFUNINIT                                                    ..  OFFINIT =e= sum(j, sqr(V(j)-Vini(j)) );
MAPPING(jGEM)     ..         VGEM(jGEM)  =e= sum((j),reacmap(jGEM,j)*v(j));  


* EMU problem
Equations
	EMUNORM(e)                  Normalization of EMUs
	emubalance1(i,e,m)          Left side of emu balances
        emubalance2(i,e,m)          Right side of emu balances
        emubalanceB(i,e,m)          Both sides of emu balances
$include "condensation_rxns_names.txt";
        LABELCOMPNORMALIT(frag)     Normalization of label comp
        EMU2FRAG(frag,e,n)          Correspondence between emus and fragments
        OBJFUNEMU                   objective function to be minimized
        LABELLIM(frag,n)            Limits of labeling with respect to measured data
        OBJFUNvar                   real objective function to be minimized
;
	
EMUNORM(e)$(not dummy(e))..sum(m$(ord(m) le (ncarbons(e)+1)),f(e,m))=e=1;
emubalance1(i,e,m)$(not dummy(e) and not excludedmet(i) and metemu(i,e) eq 1 and ord(m) le (ncarbons(e)+1) )..sum((j,ee)$(EMM(j,ee,e) gt 0), EMM(j,ee,e)*S(i,j)*V(j)*f(ee,m)) =e= bal1(i,e,m);
emubalance2(i,e,m)$(not dummy(e) and not excludedmet(i) and metemu(i,e) eq 1 and ord(m) le (ncarbons(e)+1) )..sum((j)   $(S(i,j)      gt 0),             S(i,j)*V(j)*f(e ,m)) =e= bal2(i,e,m);
emubalanceB(i,e,m)$(not dummy(e) and not excludedmet(i) and metemu(i,e) eq 1 and ord(m) le (ncarbons(e)+1) )..bal1(i,e,m)-bal2(i,e,m)=e=0;
* Combined emu equations
$include "condensation_rxns_eqns.txt";


LABELCOMPNORMALIT(frag)..sum(n,labelcomp(frag,n))=e=1;
EMU2FRAG(frag,e,n)$(aacorr(frag,e) eq 1)..labelcomp(frag,n)*labelcompnorm(frag)=e=sum(m$(ord(m) le (ncarbons(e)+1)),gamma(frag,n,m)*f(e,m));
OBJFUNEMU..OFEMU=e=sqrt(sum((fragFit,n)$(ord(n) le (ncarbonsf(fragFit)+1)), sqr( (labelcomp(fragFit,n)-labelexp(fragFit,n))/labelstd(fragFit,n) )/labelnorm(fragFit,n) )/card(fragFit));
LABELLIM(frag,n)..sqr(labelcomp(frag,n)-labelexp(frag,n))=l=sqr(labelstd(frag,n));
OBJFUNvar(fluxFocus)..OFvar=e=VGEM(fluxFocus);
























***********UPPER AND LOWER CONSTRAINTS*********

maxflow13CU = 100;

V.lo(j)                = 0;
V.up(j)                = maxflow13CU;

V.lo(jmeas)	       = vLOWER(jmeas);
V.up(jmeas)	       = vUPPER(jmeas);

labelcompnorm.lo(frag) = labelCompNormMin;
labelcompnorm.up(frag) = 1.0;
















***********INITIAL CONDITIONS*************** 

Vini(j)                = uniform(0,1);
Vini(jmeas)            = uniform(vLOWER(jmeas),vUPPER(jmeas));
V.l(j)                 = Vini(j);
V.l(jsugg)             = vSUGG(jsugg);
OFEMU.l                = 1000;





**********FEED CONSTRAINTS*****************

$include "Source_Labeling.txt";

***********MODEL DECLARATION****************

model INIT /
STOICHIOMETRY
MAPPING
OBJFUNINIT
/;
INIT.OPTFILE = 1
	
model EMU
/
STOICHIOMETRY
MAPPING
EMUNORM
emubalance1
emubalance2
emubalanceB
$include "condensation_rxns_names.txt";
EMU2FRAG
LABELCOMPNORMALIT
OBJFUNEMU
/;
EMU.OPTFILE = 1

model EMUvar
/
STOICHIOMETRY
MAPPING
EMUNORM
emubalance1
emubalance2
emubalanceB
$include "condensation_rxns_names.txt";
EMU2FRAG
LABELCOMPNORMALIT
OBJFUNEMU
LABELLIM
OBJFUNvar
/;
EMUvar.OPTFILE = 1



















	
***********MODEL SOLUTION*********************

display randparam;
display randseed;
display V.l;

* Maximization of flux
V.fx(jsugg) = VSUGG(jsugg);
Solve emu using nlp minimizing OFEMU;

V.lo(j)                = 0;
V.up(j)                = maxflow13CU;
V.lo(jmeas)	       = vLOWER(jmeas);
V.up(jmeas)	       = vUPPER(jmeas);
*V.l(jsugg)             = vSUGG(jsugg);
*Solve emu using nlp minimizing OFEMU;
*OFEMU.up = max(1,OFEMU.l);
*OFEMU.up = OFEMU.l;

Solve EMUvar using nlp maximizing OFvar;
fluxBounds('max') = OFvar.l;
Vmax(j)              = V.l(j);
labelcompMax(frag,n) = labelcomp.l(frag,n);
OFEMUMax             = OFEMU.l;

* Minimization of flux
V.fx(jsugg) = VSUGG(jsugg);
Solve emu using nlp minimizing OFEMU;

V.lo(j)                = 0;
V.up(j)                = maxflow13CU;
V.lo(jmeas)	       = vLOWER(jmeas);
V.up(jmeas)	       = vUPPER(jmeas);
*V.l(jsugg)             = vSUGG(jsugg);
*Solve emu using nlp minimizing OFEMU;
*OFEMU.up = max(1,OFEMU.l);
*OFEMU.up = OFEMU.l;

Solve EMUvar using nlp minimizing OFvar;
fluxBounds('min') = OFvar.l;
Vmin(j)              = V.l(j);
labelcompMin(frag,n) = labelcomp.l(frag,n);
OFEMUMin             = OFEMU.l;























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
solvestat                   = emu.solvestat;
modelstat                   = emu.modelstat;
numequ                      = emu.numequ;
numnz                       = emu.numnz;
numvar                      = emu.numvar;
resusd                      = emu.resusd;

    
* Objective function and general info
file OFGenInfo/"OFGenInfo.txt"/;
put OFGenInfo
put "OF = ",OFEMU.l:10:5," ";     put /;   
put "solvestat = ",solvestat," "; put /;
put "modelstat = ",modelstat," "; put /;   
put "numequ    = ",numequ," ";    put /;
put "numnz     = ",numnz," ";     put /;
put "numvar    = ",numvar," ";    put /;
put "resusd    = ",resusd," ";    put /;    
putclose OFGenInfo;


file OFGenInfoMax/"OFGenInfoMax.txt"/;
put OFGenInfoMax
put "OF = ",OFEMUMax:10:5," ";     put /;   
put "solvestat = ",solvestat," "; put /;
put "modelstat = ",modelstat," "; put /;   
put "numequ    = ",numequ," ";    put /;
put "numnz     = ",numnz," ";     put /;
put "numvar    = ",numvar," ";    put /;
put "resusd    = ",resusd," ";    put /;    
putclose OFGenInfoMax;

file OFGenInfoMin/"OFGenInfoMin.txt"/;
put OFGenInfoMin
put "OF = ",OFEMUMin:10:5," ";     put /;   
put "solvestat = ",solvestat," "; put /;
put "modelstat = ",modelstat," "; put /;   
put "numequ    = ",numequ," ";    put /;
put "numnz     = ",numnz," ";     put /;
put "numvar    = ",numvar," ";    put /;
put "resusd    = ",resusd," ";    put /;    
putclose OFGenInfoMin;

* Fluxes output
file Voutfile/"Vout.txt"/;
put Voutfile
loop(j,
	put "v('",j.tl,"') = ",V.l(j):20:15," "; 	put /;	
     );
putclose Voutfile;

file Vmaxoutfile/"Vmaxout.txt"/;
put Vmaxoutfile
loop(j,
	put "vmax('",j.tl,"') = ",Vmax(j):20:15," "; 	put /;	
     );
putclose Vmaxoutfile;

file Vminoutfile/"Vminout.txt"/;
put Vminoutfile
loop(j,
	put "vmin('",j.tl,"') = ",Vmin(j):20:15," "; 	put /;	
     );
putclose Vminoutfile;



* EMU MDVs output
file foutfile/"fout.txt"/;
put foutfile
loop(e,
	loop(m$(ord(m) le (ncarbons(e)+1)),put "f('",e.tl:40,",",m.tl"') = ",f.l(e,m):10:5," "/);
     );
putclose foutfile;


* Simulated MS MDVs output
file labelcompoutfile/"labelcomp.txt"/;
put labelcompoutfile
loop(frag,
	loop(n,put "labelcomp('",frag.tl:40,",",n.tl"') = ",labelcomp.l(frag,n):10:5," "/);
     );
putclose labelcompoutfile;

file labelcompMaxoutfile/"labelcompMax.txt"/;
put labelcompMaxoutfile
loop(frag,
	loop(n,put "labelcompMax('",frag.tl:40,",",n.tl"') = ",labelcompMax(frag,n):10:5," "/);
     );
putclose labelcompMaxoutfile;

file labelcompMinoutfile/"labelcompMin.txt"/;
put labelcompMinoutfile
loop(frag,
	loop(n,put "labelcompMin('",frag.tl:40,",",n.tl"') = ",labelcompMin(frag,n):10:5," "/);
     );
putclose labelcompMinoutfile;

* Flux bounds 
file fluxboundsfile/"fluxBounds.txt"/;
put fluxboundsfile
*put "fluxBounds('min') = 0.0"/;
*put "fluxBounds('max') = 0.0"/;
put "fluxBounds('min') = ",fluxBounds('min'):10:5," "/;
put "fluxBounds('max') = ",fluxBounds('max'):10:5," "/;
putclose fluxboundsfile;
