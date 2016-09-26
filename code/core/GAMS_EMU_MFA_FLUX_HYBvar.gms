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
	maxflow13CL    maximum flux value for 13C problem (lower value)
	maxflow13CU    maximum flux value for 13C problem (upper value)	
;



$include "StoichNetwork.txt"
$include "vUPPER.txt"
$include "vLOWER.txt"
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
	iGEM    set all metabolites for GEM
$include "hybridmets.txt"
	jGEM    set all reactions for GEM
$include "hybridreacs.txt"
        jGEMlink(jGEM) set of reactions GEM linked to EMU transitions
$include "linkedreacs.txt"
        jGEMsugg(jGEM) set of suggested reactions
$include "jGEMsugg.txt"
;

display jGEM;


parameters
	Sbig(iGEM,jGEM)  Stochiometry matrix for GEM
	bvec(iGEM)       b vector for FBA problem (metabolite accumulation)
	lbvec(jGEM)      lower bound for flux vector
	ubvec(jGEM)      upper bound for flux vector
	cvec(jGEM)       optimization vector for FBA
	reacmap(jGEM,j)  reaction map between GEM reactions and EMU transitions
	VFBAgrowth(jGEM) Stores FBA solution for maximum growth
	Vgrowth(j)       Mapping of FBA maximum growth solution into EMU fluxes
	vFBASUGG(jGEMsugg)  suggested fluxes
	VFBAmax(jGEM)    FBA flux values for maximum of fluxfocus
	VFBAmin(jGEM)    FBA flux values for minimum of fluxfocus
	OFEMUMax         EMU objective function for max fluxfocus
	OFEMUMin         EMU objective function for max fluxfocus
	labelcompMax(frag,n)       Simulated MS data for maximum fluxfocus
	labelcompMin(frag,n)       Simulated MS data for mininum fluxfocus
	maxflowGEM        Maximum allowed flux for GEM reactions
;

$include "Sbig.txt"
$include "bvec.txt"
$include "lbvec.txt"
$include "ubvec.txt"
$include "cvec.txt"
$include "reacmap.txt"
$include "vFBASUGG.txt"


Scalar inFluxRef;
$include "inFluxRef.txt"

***************EXTRA SETS&PARAMETERS******************

sets
	fluxFocus(jGEM)  flux focused on for variability analysis
$include "fluxFocus.txt"
        bounds bounds for labeling focused on /'min','max'/
;





* Random start parameters and others
Parameters
	randparam(j)
	fluxBounds(bounds)
;







Scalar randseed;
$include "randseed.txt"
execseed = randseed;


loop(j,
	randparam(j)= uniform(0,1)
);

Scalar labelCompNormMin;
$include "labelCompNormMin.txt"




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
	OFFBA                     Initialization objective function
	OFFBAgrowth               Max growth objective function
;

Positive Variables
	VGEM(jGEM)                GEM fluxes 
;

***********EQUATIONS**************

* INIT
Equations
	STOICHIOMETRY(i)            Stoichiometric flux conservation constraint
	OBJFUNINIT                  Objective function for initialization (as close as possible to Vini meeting stoich constraints)
;



STOICHIOMETRY(i)$(i(i) - exchangemet(i) - excludedmet(i))     .. 0 =e= sum(j,S(i,j)*V(j));
OBJFUNINIT                                                    .. OFFINIT =e= sum(j, sqr(V(j)-Vini(j)) );



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
        OBJFUNvar(fluxFocus)        real objective function to be minimized
;
	
EMUNORM(e)$(not dummy(e))..sum(m$(ord(m) le (ncarbons(e)+1)),f(e,m))=e=1;
emubalance1(i,e,m)$(not dummy(e) and not excludedmet(i) and metemu(i,e) eq 1 and ord(m) le (ncarbons(e)+1) )..sum((j,ee)$(EMM(j,ee,e) gt 0), EMM(j,ee,e)*S(i,j)*V(j)*f(ee,m)) =e= bal1(i,e,m);
emubalance2(i,e,m)$(not dummy(e) and not excludedmet(i) and metemu(i,e) eq 1 and ord(m) le (ncarbons(e)+1) )..sum((j)   $(S(i,j)      gt 0),             S(i,j)*V(j)*f(e ,m)) =e= bal2(i,e,m);
emubalanceB(i,e,m)$(not dummy(e) and not excludedmet(i) and metemu(i,e) eq 1 and ord(m) le (ncarbons(e)+1) )..bal1(i,e,m)-bal2(i,e,m)=e=0;
* Combined emu equations
$include "condensation_rxns_eqns.txt";


LABELCOMPNORMALIT(frag)..sum(n$(ord(n) le (ncarbonsf(frag)+1)),labelcomp(frag,n))=e=1;
EMU2FRAG(frag,e,n)$(aacorr(frag,e) eq 1 and ord(n) le (ncarbonsf(frag)+1) )..labelcomp(frag,n)=e=sum(m$(ord(m) le (ncarbons(e)+1)),gamma(frag,n,m)*f(e,m))/labelcompnorm(frag);
OBJFUNEMU..OFEMU=e=sqrt(sum((fragFit,n)$(ord(n) le (ncarbonsf(fragFit)+1)), sqr( (labelcomp(fragFit,n)-labelexp(fragFit,n))/labelstd(fragFit,n) )/labelnorm(fragFit,n) )/card(fragFit));
LABELLIM(frag,n)..sqr(labelcomp(frag,n)-labelexp(frag,n))=l=sqr(labelstd(frag,n));
OBJFUNvar(fluxFocus)..OFvar=e=VGEM(fluxFocus);


* HYBRID
Equations
	MAPPING(jGEM)                  Mapping between transitions and GEM reactions
	STOICHIOMETRYGEM(iGEM)         Stoichiometric flux conservation constraint for GEM
	OBJFUNFBA                      Standard FBA objective function for large problem
	OBJFUNFBAgrowth                Standard FBA objective function maximizing growth for comparison
;


MAPPING(jGEM)$(jGEMlink(jGEM))     ..         VGEM(jGEM)  =e= inFluxRef*sum((j),reacmap(jGEM,j)*v(j));  
STOICHIOMETRYGEM(iGEM)             ..         bvec(iGEM)  =e= sum((jGEM),Sbig(iGEM,jGEM)*VGEM(jGEM)); 
OBJFUNFBA                          ..         OFFBA       =e= sum((jGEM),cvec(jGEM)*VGEM(jGEM));
OBJFUNFBAgrowth                    ..         OFFBAgrowth =e= VGEM('BiomassEcoli');


***********UPPER AND LOWER CONSTRAINTS*********
maxflow13CL = 10;
maxflow13CU = 100;
maxflowGEM  = 1000;

V.lo(j)                = 0;
V.up(j)                = maxflow13CU;

*        Capping bounds to numerically acceptable values
vLOWER(jmeas)	       = max(vLOWER(jmeas),-maxflow13CU);         
vUPPER(jmeas)	       = min(vUPPER(jmeas), maxflow13CU);

labelcompnorm.lo(frag) = labelCompNormMin;
labelcompnorm.up(frag) = 1.0;

lbvec(jGEM)          = max(lbvec(jGEM),-maxflowGEM);
ubvec(jGEM)          = min(ubvec(jGEM), maxflowGEM);

V.lo(jmeas)	       = vLOWER(jmeas);
V.up(jmeas)	       = vUPPER(jmeas);
VGEM.lo(jGEM)          = lbvec(jGEM);
VGEM.up(jGEM)          = ubvec(jGEM);






***********INITIAL CONDITIONS*************** 
VGEM.l(jGEM)           = uniform(0,1);

Vini(j)                = uniform(0,1);
Vini(jmeas)            = uniform(vLOWER(jmeas),max(vLOWER(jmeas),min(10,vUPPER(jmeas))));
V.l(j)                 = Vini(j);
V.l(jsugg)             = vSUGG(jsugg);
f.l(e,m)               = uniform(0,1);
labelcomp.l(frag,n)    = uniform(0,1);
OFEMU.l                = 1000;

**********FEED CONSTRAINTS*****************

$include "Source_Labeling.txt";

***********MODEL DECLARATION****************

model INIT /
STOICHIOMETRY
OBJFUNINIT
/;
model INIT2 /
STOICHIOMETRYGEM
OBJFUNINIT
MAPPING
/;
INIT2.OPTFILE = 1

model EMU
/
*STOICHIOMETRY
EMUNORM
emubalance1
emubalance2
emubalanceB
$include "condensation_rxns_names.txt";
EMU2FRAG
LABELCOMPNORMALIT
OBJFUNEMU
MAPPING
STOICHIOMETRYGEM
/;
EMU.OPTFILE = 1

model INITFBA /
STOICHIOMETRYGEM
OBJFUNFBA
MAPPING
/;
INITFBA.OPTFILE = 1

model FBAgrowth /
STOICHIOMETRYGEM
OBJFUNFBAgrowth
MAPPING
/;
FBAgrowth.OPTFILE = 1

model EMUvar
/
*STOICHIOMETRY
EMUNORM
emubalance1
emubalance2
emubalanceB
$include "condensation_rxns_names.txt";
EMU2FRAG
LABELCOMPNORMALIT
OBJFUNEMU
LABELLIM
MAPPING
STOICHIOMETRYGEM
OBJFUNvar
/;

***********MODEL SOLUTION*********************

display randparam;
display randseed;

** FBA max growth problem
* Free growth
VGEM.lo('BiomassEcoli')          = 0;
VGEM.up('BiomassEcoli')          = maxflow13CU;
Solve FBAgrowth using lp maximizing OFFBAgrowth;
VGEM.lo('BiomassEcoli')          = lbvec('BiomassEcoli');
VGEM.up('BiomassEcoli')          = ubvec('BiomassEcoli');

VFBAgrowth(jGEM) = VGEM.l(jGEM);
Vgrowth(j)       = V.l(j);

** FBA initialization
*Solve INITFBA using lp maximizing OFFBA;

** Initialization part I
*Solve INIT2 using nlp minimizing OFFINIT;
** Initilalization part II
*V.fx(j) = V.l(j);

* Maximization of flux
VGEM.fx(jGEMsugg) = vFBASUGG(jGEMsugg);
V.fx(jsugg)       = vSUGG(jsugg);
Solve emu using nlp minimizing OFEMU;

V.lo(j)                = 0;
V.up(j)                = maxflow13CU;
V.lo(jmeas)	       = vLOWER(jmeas);
V.up(jmeas)	       = vUPPER(jmeas);
VGEM.lo(jGEM)          = lbvec(jGEM);
VGEM.up(jGEM)          = ubvec(jGEM);
*OFEMU.up = max(1,OFEMU.l);
*OFEMU.up = OFEMU.l;

Solve EMUvar using nlp maximizing OFvar;
fluxBounds('max') = OFvar.l;
VFBAmax(jGEM)        = VGEM.l(jGEM);
labelcompMax(frag,n) = labelcomp.l(frag,n);
OFEMUMax             = OFEMU.l;


	
* Minimization of flux
VGEM.fx(jGEMsugg) = vFBASUGG(jGEMsugg);
V.fx(jsugg)       = vSUGG(jsugg);
Solve emu using nlp minimizing OFEMU;

V.lo(j)                = 0;
V.up(j)                = maxflow13CU;
V.lo(jmeas)	       = vLOWER(jmeas);
V.up(jmeas)	       = vUPPER(jmeas);
VGEM.lo(jGEM)          = lbvec(jGEM);
VGEM.up(jGEM)          = ubvec(jGEM);
*OFEMU.up = max(1,OFEMU.l);
*OFEMU.up = OFEMU.l;

Solve EMUvar using nlp minimizing OFvar;
fluxBounds('min') = OFvar.l;
VFBAmin(jGEM)        = VGEM.l(jGEM);
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


* EMU Fluxes output
file Voutfile/"Vout.txt"/;
put Voutfile
loop(j,
	put "v('",j.tl,"') = ",V.l(j):20:15," "; 	put /;	
     );
putclose Voutfile;


* EMU MDVs output
file foutfile/"fout.txt"/;
put foutfile
loop(e,
	loop(m,put "f('",e.tl:40,",",m.tl"') = ",f.l(e,m):10:5," "/);
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


* FBA Fluxes output
file VGEMoutfile/"VGEMout.txt"/;
put VGEMoutfile
loop(jGEM,
	put "vGEM('",jGEM.tl:30:10,"') = ",VGEM.l(jGEM):20:15," "; 	put /;	
     );
putclose VGEMoutfile;

file VFBAmaxoutfile/"VFBAmaxout.txt"/;
put VFBAmaxoutfile
loop(jGEM,
	put "vFBAmax('",jGEM.tl:30:10,"') = ",VFBAmax(jGEM):20:15," "; 	put /;	
     );
putclose VFBAmaxoutfile;

file VFBAminoutfile/"VFBAminout.txt"/;
put VFBAminoutfile
loop(jGEM,
	put "vFBAmin('",jGEM.tl:30:10,"') = ",VFBAmin(jGEM):20:15," "; 	put /;	
     );
putclose VFBAminoutfile;


* Initial FBA output
file VFBAgrowthoutfile/"VFBAgrowthout.txt"/;
put VFBAgrowthoutfile
loop(jGEM,
	put "vFBAgrowth('",jGEM.tl:30:10,"') = ",VFBAgrowth(jGEM):20:15," "; 	put /;	
     );
putclose VFBAgrowthoutfile;

* EMU Fluxes output for growth case
file Vgrowthoutfile/"Vgrowthout.txt"/;
put Vgrowthoutfile
loop(j,
	put "vgrowth('",j.tl:30:10,"') = ",Vgrowth(j):20:15," "; 	put /;	
     );
putclose Vgrowthoutfile;

* Flux bounds 
file fluxboundsfile/"fluxBounds.txt"/;
put fluxboundsfile
put "fluxBounds('min') = ",fluxBounds('min'):10:5," "/;
put "fluxBounds('max') = ",fluxBounds('max'):10:5," "/;
putclose fluxboundsfile;
