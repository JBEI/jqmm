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

***********FBA SETS&PARS***************

* Stoichiometry
sets
	i  all metabolites
$include "mets.txt"
	j  reaction indices /r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16,r17,r18,r19,r20,r21/
	a  alpha coefficient index /a1,a2,a3,a4,a5,a6,a7/
        jmeas(j) set of measured reactions /r10,r14/;       

parameters
	S(i,j)        Stoichiometry matrix
	S2(i,j)       Stoichiometry matrix for second case (pyruvate)
	alpha(a)  alpha coefficients
	/
	a1  0.32
	a2  5.20
	a3  6.20
	a4  7.07
	a5  3.11
	a6  4.27
	a7  2.68
	/
;
$include "StoichNetwork.txt";
$include "StoichNetwork2.txt";

	
* EMUs
sets
	e   set of all emus	
$include "all_emu.txt"
	dummy(e) dummy emus for condensation reacs
$include "dummy_emus.txt"
	m values of m for emu MDVs
$include "mset.txt"
;

alias(e,ee);

parameters
	ncarbons(e)        number of carbons for each em
	metemu(i,e)        1 if e is an emu of i and zero otherwise
	EMM                EMU mapping matrix
	EMM2               EMU mapping matrix 2
	vUPPER(jmeas) upper boundary for fluxes 
	vLOWER(jmeas) lower boundary for fluxes 
;
$include "emu_carbon_numbers.txt";
$include "met_emu.txt";
$include "emm_simplified.txt";
$include "emm_simplified2.txt";
$include "vUPPER.txt";
$include "vLOWER.txt";

* Ratios
parameters
	ratioGlu   Target ratio for glucose input
	ratioPyr   Target ratio for pyruvate input
;
$include "ratios.txt"

* Variables	
Positive Variables
	V(j)       fluxes for reaction j
	fA(e,m)     MDV(m) for emu e in scenario a
	fB(e,m)     MDV(m) for emu e in scenario b
	fC(e,m)     MDV(m) for emu e in scenario c
	fD(e,m)     MDV(m) for emu e in scenario d
;

Free Variables 
	OF          Objective function
	bal1A
	bal2A
	bal1B
	bal2B
	bal1C
	bal2C
	bal1D
	bal2D
;

Equations
	st1
        st2
	st3
	st4
	st5
	st6
	st7
	st8
	st9
	st11
	st12
	st13
	st15
	st16
	st17
	st18
	st19
	st20
	st21
;
	
st1  .. V('r1')  =e= 100 ;
st2  .. V('r2')  =e= V('r1')- V('r10') - V('r14') ;
st3  .. V('r3')  =e= 2*(V('r2') + 2.0/3.0*V('r11') - V('r15'));
st4  .. V('r4')  =e= V('r3') + V('r11')/3.0 - V('r16') ;
st5  .. V('r5')  =e= V('r4') - V('r9')- V('r17');
st6  .. V('r6')  =e= V('r5') - V('r18');
st7  .. V('r7')  =e= V('r6');
st8  .. V('r8')  =e= V('r6') - V('r19');
st9  .. V('r9')  =e= V('r19') + V('r20');
st11 .. V('r11') =e= V('r10') - V('r21');
st12 .. V('r12') =e= 2.0/3.0 * V('r11');
st13 .. V('r13') =e= 1.0/3.0 * V('r11');
st15 .. V('r15') =e= 0.7;
st16 .. V('r16') =e= 19.8;
st17 .. V('r17') =e= 24.7;
st18 .. V('r18') =e= 28.1;
st19 .. V('r19') =e= 10.7;
st20 .. V('r20') =e= 15.8;
st21 .. V('r21') =e= 10.4;

* EMUs
Equations
	EMUNORMA(e)                 Normalization of EMUs scenario A
	EMUNORMB(e)                 Normalization of EMUs scenario B 
	EMUNORMC(e)                 Normalization of EMUs scenario C
	EMUNORMD(e)                 Normalization of EMUs scenario D 
	emubalance1A(i,e,m)          Left side of emu balances
	emubalance1B(i,e,m)          Left side of emu balances
	emubalance1C(i,e,m)          Left side of emu balances
	emubalance1D(i,e,m)          Left side of emu balances	
        emubalance2A(i,e,m)          Right side of emu balances
        emubalance2B(i,e,m)          Right side of emu balances
        emubalance2C(i,e,m)          Right side of emu balances
        emubalance2D(i,e,m)          Right side of emu balances
	emubalanceTA(i,e,m)          Both sides of emu balances
	emubalanceTB(i,e,m)          Both sides of emu balances
	emubalanceTC(i,e,m)          Both sides of emu balances
	emubalanceTD(i,e,m)          Both sides of emu balances
**$include "condensation_rxns_names.txt";
;
	
EMUNORMA(e)$(not dummy(e))..sum(m$(ord(m) le (ncarbons(e)+1)),fA(e,m))=e=1;
emubalance1A(i,e,m)$(not dummy(e) and metemu(i,e) eq 1 and ord(m) le (ncarbons(e)+1) )..sum((j,ee)$(EMM(j,ee,e) gt 0), EMM(j,ee,e)*S(i,j)*V(j)*fA(ee,m)) =e= bal1A(i,e,m);
emubalance2A(i,e,m)$(not dummy(e) and metemu(i,e) eq 1 and ord(m) le (ncarbons(e)+1) )..sum((j)   $(S(i,j)      gt 0),             S(i,j)*V(j)*fA(e ,m)) =e= bal2A(i,e,m);
emubalanceTA(i,e,m)$(not dummy(e) and metemu(i,e) eq 1 and ord(m) le (ncarbons(e)+1) )..bal1A(i,e,m)-bal2A(i,e,m)=e=0;
** Combined emu equations
**$include "condensation_rxns_eqns.txt";
EMUNORMB(e)$(not dummy(e))..sum(m$(ord(m) le (ncarbons(e)+1)),fB(e,m))=e=1;
emubalance1B(i,e,m)$(not dummy(e) and metemu(i,e) eq 1 and ord(m) le (ncarbons(e)+1) )..sum((j,ee)$(EMM(j,ee,e) gt 0), EMM(j,ee,e)*S(i,j)*V(j)*fB(ee,m)) =e= bal1B(i,e,m);
emubalance2B(i,e,m)$(not dummy(e) and metemu(i,e) eq 1 and ord(m) le (ncarbons(e)+1) )..sum((j)   $(S(i,j)      gt 0),             S(i,j)*V(j)*fB(e ,m)) =e= bal2B(i,e,m);
emubalanceTB(i,e,m)$(not dummy(e) and metemu(i,e) eq 1 and ord(m) le (ncarbons(e)+1) )..bal1B(i,e,m)-bal2B(i,e,m)=e=0;
** Combined emu equations
**$include "condensation_rxns_eqns.txt";
EMUNORMC(e)$(not dummy(e))..sum(m$(ord(m) le (ncarbons(e)+1)),fC(e,m))=e=1;
emubalance1C(i,e,m)$(not dummy(e) and metemu(i,e) eq 1 and ord(m) le (ncarbons(e)+1) )..sum((j,ee)$(EMM2(j,ee,e) gt 0) , EMM2(j,ee,e)*S2(i,j)*V(j)*fC(ee,m)) =e= bal1C(i,e,m);
emubalance2C(i,e,m)$(not dummy(e) and metemu(i,e) eq 1 and ord(m) le (ncarbons(e)+1) )..sum((j)   $(S2(i,j)       gt 0),              S2(i,j)*V(j)*fC(e ,m)) =e= bal2C(i,e,m);
emubalanceTC(i,e,m)$(not dummy(e) and metemu(i,e) eq 1 and ord(m) le (ncarbons(e)+1) )..bal1C(i,e,m)-bal2C(i,e,m)=e=0;
** Combined emu equations
**$include "condensation_rxns_eqns.txt";
EMUNORMD(e)$(not dummy(e))..sum(m$(ord(m) le (ncarbons(e)+1)),fD(e,m))=e=1;
emubalance1D(i,e,m)$(not dummy(e) and metemu(i,e) eq 1 and ord(m) le (ncarbons(e)+1) )..sum((j,ee)$(EMM2(j,ee,e) gt 0) , EMM2(j,ee,e)*S2(i,j)*V(j)*fD(ee,m)) =e= bal1D(i,e,m);
emubalance2D(i,e,m)$(not dummy(e) and metemu(i,e) eq 1 and ord(m) le (ncarbons(e)+1) )..sum((j)   $(S2(i,j)       gt 0),              S2(i,j)*V(j)*fD(e ,m)) =e= bal2D(i,e,m);
emubalanceTD(i,e,m)$(not dummy(e) and metemu(i,e) eq 1 and ord(m) le (ncarbons(e)+1) )..bal1D(i,e,m)-bal2D(i,e,m)=e=0;
** Combined emu equations
**$include "condensation_rxns_eqns.txt";

* Objective function
Equations
	OFeq
;

*OFeq .. OF =e= V('r14');
*OFeq .. OF =e= sqr( fA('co2_c_1','1')-2.48*fB('co2_c_1','1') ) + sqr( fC('co2_c_1','1')-2.00*fD('co2_c_1','1') ) ;
OFeq .. OF =e= sqr( fA('co2_c_1','1')-ratioGlu*fB('co2_c_1','1') ) + sqr( fC('co2_c_1','1')-ratioPyr*fD('co2_c_1','1') ) ;


* Starting conditions
V.lo(j)=0.01;
V.up(j)=200.0;
	
V.lo(jmeas)	       = vLOWER(jmeas);
V.up(jmeas)	       = vUPPER(jmeas);

V.l(j) = uniform(0.01,100);


*FEED CONSTRAINTS
$include "Source_LabelingA.txt";
$include "Source_LabelingB.txt";
$include "Source_LabelingC.txt";
$include "Source_LabelingD.txt";

* Model declaration
model first /
      st1
      st2
      st3
      st4
      st5
      st6
      st7
      st8
      st9
      st11
      st12
      st13
      st15
      st16
      st17
      st18
      st19
      st20
      st21
      OFeq
      EMUNORMA
      EMUNORMB  
      EMUNORMC
      EMUNORMD      
      emubalance1A
      emubalance1B
      emubalance1C
      emubalance1D
      emubalance2A
      emubalance2B
      emubalance2C
      emubalance2D
      emubalanceTA
      emubalanceTB
      emubalanceTC
      emubalanceTD
/;

* Solve model*
Solve first using nlp minimizing OF;

* Tests
display alpha
display V.l

****OUTPUTS****
Parameters
        solvestat                       Exit status of solver
        modelstat                       Exit status of model
        numequ                          Number of equations in model 
        numnz                           Number of nonzero elements in model
        numvar                          Number of variables in model
        resusd                          Seconds used to reach final model state

;

* General info output
solvestat                   = first.solvestat;
modelstat                   = first.modelstat;
numequ                      = first.numequ;
numnz                       = first.numnz;
numvar                      = first.numvar;
resusd                      = first.resusd;

	
* Objective function and general info
file OFGenInfo/"OFGenInfo.txt"/;
put OFGenInfo
put "OF = ",OF.l:10:5," ";     put /;   
put "solvestat = ",solvestat," "; put /;
put "modelstat = ",modelstat," "; put /;   
put "numequ    = ",numequ," ";    put /;
put "numnz     = ",numnz," ";     put /;
put "numvar    = ",numvar," ";    put /;
put "resusd    = ",resusd," ";    put /;    
putclose OFGenInfo;
	
* EMU MDVs output
file foutfileA/"fAout.txt"/;
put foutfileA
loop(e,
	loop(m$(ord(m) le (ncarbons(e)+1)),put "fA('",e.tl:40,",",m.tl"') = ",fA.l(e,m):10:5," "/);
     );
putclose foutfileA;
file foutfileB/"fBout.txt"/;
put foutfileB
loop(e,
	loop(m$(ord(m) le (ncarbons(e)+1)),put "fB('",e.tl:40,",",m.tl"') = ",fB.l(e,m):10:5," "/);
     );
putclose foutfileB;
file foutfileC/"fCout.txt"/;
put foutfileC
loop(e,
	loop(m$(ord(m) le (ncarbons(e)+1)),put "fC('",e.tl:40,",",m.tl"') = ",fC.l(e,m):10:5," "/);
     );
putclose foutfileC;
file foutfileD/"fDout.txt"/;
put foutfileD
loop(e,
	loop(m$(ord(m) le (ncarbons(e)+1)),put "fD('",e.tl:40,",",m.tl"') = ",fD.l(e,m):10:5," "/);
     );
putclose foutfileD;

* Fluxes output
file Voutfile/"Vout.txt"/;
put Voutfile
loop(j,
	put "V('",j.tl,"') = ",V.l(j):20:15," "; 	put /;	
     );
putclose Voutfile;

