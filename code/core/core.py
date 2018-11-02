from __future__ import print_function
from __future__ import division
# Part of the JBEI Quantitative Metabolic Modeling Library (JQMM)
# Copyright (c) 2016, The Regents of the University of California.
# For licensing details see "license.txt" and "legal.txt".

from builtins import str
from builtins import range
from utilities import old_div
from builtins import object
import re, copy, numpy, math, string, unittest
import enhancedLists, Genes, Proteins, DB
import utilities as utils


"##################################"
" METABOLITES, REACTIONS AND EMUS"
"##################################"

# List of allowed compartment types stored in DB
compartments = DB.getCompartmentDict() 

# Used to set upper and lower bounds on fluxes.
# These values are used when the bound is declared to be present by SBML markup,
# but the markup does not specify a value (or provides an un-parseable one).
DEFAULT_LOWER_BOUND = 1000
DEFAULT_UPPER_BOUND = 1000

# Used in Reaction class, in getSignificantTranscriptomicsValue and getSignificantProteomicsValue:
TRANSCRIPTOMICS_SIGNIFICANCE_THRESHOLD = 80 # Must be 80 or above
PROTEOMICS_SIGNIFICANCE_THRESHOLD = 0       # Must be above 0


class Metabolite(object):
    """
    Class for metabolites. This class has been created for metabolites involved in 13C MFA and FBA. 
    In init:
        -) name: metabolite name as a MetaboliteName object, or a string in standard format e.g. ala-L[c]
        -) ncarbons: number of carbons in metabolite (this is redundant with formula, eliminate)
        -) excluded: whether metabolite is excluded from 13C EMU balance
        -) source: is it a source metabolite (source of labeling)?
        -) feed: labeling feed, e.g. 30% 1-C 20% U 50% UN meaning 30% labeled in first carbon, 20% fullly labeled and 50% unlabeled
        -) destination: is it a destination metabolite? e.g. is it a metabolite for which labeling is measured (destination for label)?
        -) carbonsOut: indices for the measured MDV, e.g.  [[0,1],[2,3]] for tyrL 0,1 and tyrL 2,3. tyrL 2,3 indicates that labeling for atoms 2 and 3 is measured
        -) bval: value of b vector from LP problem (sum_j S_{ij}*v_j = b_i). Used for FBA
        -) full_name: full name, e.g. L-alanine for ala-L[c]
        -) formula: composition formula, e.g C3H7NO2 for ala-L[c]
        -) extraNotes: extra notes information
    """
 
    def __init__(self, name, ncarbons=0, excluded=False, source=False, feed='none', destination=False, carbonsOut=None, bval=0, full_name=None,
                 formula='', extraNotes=None):

        self.name = name
        if isinstance(name, MetaboliteName):
            self.compartment = name.compartment
        else:
            # Attempt to extract a compartment name (or the default) if the given name is not a MetaboliteName already.
            temp_name = MetaboliteName(name)
            self.compartment = temp_name.compartment

        self.ncarbons    = ncarbons
        self.excluded    = excluded
        self.source      = source
        self.feed        = feed
        self.destination = destination
        self.carbonsOut  = carbonsOut
        self.bval        = bval
        self.extraNotes  = {}
        if extraNotes != None:
            self.extraNotes  = extraNotes
        self.full_name   = name
        if full_name is not None:
            if full_name.strip() != '':
                self.full_name = full_name
        self.formula     = formula


    def __str__(self):
        return str(self.name)


    def generateEMU(self, originalExceptions):
        """
        Generates EMU corresponding to Metabolite and given exceptions.
        Exceptions is a list looking like e.g. ['0','1','2'] indicating the carbons to be excluded from the EMU.
        For example:
            ['0'] --> f6p_1_2_3_4_5_6
            ['0','1','2'] --> f6p_3_4_5_6
            ..... etc
        """
        exceptions = copy.deepcopy(originalExceptions)
        # Eliminate 0 if present
        try:
            exceptions.remove('0')
        except ValueError:
            pass
        
        # Obtain EMU carbon numbering
        ncarbons = self.ncarbons
        if ncarbons <= 0:
            raise Exception('Metabolite %s has %s carbon atoms' % (str(self.name), str(ncarbons)))
            
        clist = list(range(1,ncarbons+1))
        for exc in exceptions:
            try:
                clist.remove(int(exc))
            except ValueError:
                raise Exception('%s not found in exceptions: %s' % (str(exc), str(exceptions)))

        sclist =[]
        for memb in clist:
            sclist.append(str(memb))
            
        emuName = self.name+'_' + '_'.join(sclist)
        return emuName



class Reactant(Metabolite):
    "reactant class, this is a metabolite in the context of a reaction (where stoichiometry makes sense)."
    def __init__(self, metabolite, stoichiometry, label='', extraLabel=[]):
        self.stoichiometry = stoichiometry
        self.originalMet   = metabolite
        if label:
            self.label = label
        self.extraLabel = extraLabel

        Metabolite.__init__(self, metabolite.name, metabolite.ncarbons, metabolite.excluded, metabolite.source, metabolite.feed, metabolite.destination, metabolite.carbonsOut)


        
class Product(Reactant):
    "product class. Same as above, a metabolite in the context of a reaction."



class Reaction(object):    
    """
    Class for reactions. This class has been created for reactions involved in 13C MFA and FBA.
    In init:
        -) name: reaction name e.g. PDH
        -) reversible: is reaction reversible?
        -) suggested: is this reaction flux part of a set of suggested initial fluxes for optimization (probably eliminate this one)
        -) measured: is this reaction flux measured (e.g. exchange flux)
        -) flux: flux for reaction in the form of an instance of the flux class
        -) FluxBounds: flux bounds this the flux associated to this reaction
        -) transitionLine: carbon transition line e.g. e.g.: ILETA	   akg + ile-L[c] <==> glu-L[c] + 3mop	ABCDE + abcdef : ABCDE + abcdef 
        -) exchange: is this an exchange reaction?
        -) reactants: reactants in reaction (instances of reactant class)
        -) products: products in reaction
        -) cval: c value in in LP problem of type max c_j*v_j subject to sum_j S_{ij}*v_j=0 (This only applies to FBA reactions)
        -) subsystem: subsystem the reaction belongs to, according to BIGG database, e.g. S_GlycolysisGluconeogenesis
        -) gene: gene(s) associated to reaction
        -) protein: protein(s) associated to reaction
        -) extraNotes: extra notes associated to reaction
        
    """
    def __init__(self, name, reversible=False, suggested=False, measured=False, flux='None', fluxBounds='None',
                 transitionLine='None', exchange=False, reactants='default', products='default', cval='default', 
                 subsystem='None', gene=None, protein=None, extraNotes=None):
        self.name = name
        self.reversible = reversible
        self.suggested  = suggested
        self.measured   = measured
        self.extraNotes = extraNotes if extraNotes != None else {}


        if flux != 'None':
            self.flux = flux
        if fluxBounds != 'None':
            self.fluxBounds = fluxBounds
        if transitionLine != 'None':
            self.transitionLine = AtomTransition(transitionLine)
        self.exchange = exchange
        if reactants != 'default':
            self.reactants = reactants  
            self.getReactDict()
        if products != 'default':
            self.products = products
            self.getProdDict()
        if cval != 'default':
            self.cval = cval
        self.subsystem = subsystem.strip()

        self.setGenes(gene)
        self.setProteins(protein)


    @classmethod
    def from_tuple(cls, tup):
        """
        Construct a new Reaction based on tuple input: (name, rev, exchange, reacts, prods)
        """
        name, rev, exchange, reacts, prods = tup        
        
        reactants = []
        products  = []
        for react in reacts:
            stoich,react,comp = react
            metabolite = Metabolite(react.strip()+'_'+comp, ncarbons=0)
            reactants.append(Reactant(metabolite, stoich))
        for prod in prods:
            stoich,prod,comp = prod
            metabolite = Metabolite(prod.strip()+'_'+comp, ncarbons=0)
            products.append(Product(metabolite,stoich))
        
        return cls(name, reversible=rev, exchange=exchange, reactants=reactants, products=products, subsystem='added')


    @classmethod
    def from_string(cls, string):
        """
        Construct a new Reaction based on string input e.g.:             
        
        ILETA	   akg + ile-L[c] <==> glu-L[c] + 3mop	ABCDE + abcdef : ABCDE + abcdef
        """
        name, rev, exchange, reacts, prods = utils.parseLine(string)
        
        return cls.from_tuple((name, rev, exchange, reacts, prods))


    def __str__(self):
        "Prints name and flux value"
        string_ = self.name
        if hasattr(self,'flux'):
            string_ = string_ + ': ' +str(self.flux)
        return string_


    def __mul__(self,other):
        "Reaction multiplication just means multiply flux and fluxbounds by factor"
        self.flux       = other*self.flux        
        self.fluxBounds = other*self.fluxBounds


    def __rmul__(self,other):
        "Reaction reverse multiplication same as multiplication (i.e. commutative)"
        result = self.__mul__(other)
        return result


    def stoichLine(self, SBML=False):
        """
        Produces stochiometry lines for the reactions, e.g.: PDH : coa_c + nad_c + pyr_c --> accoa_c  + co2_c + nadh_c
        if SBML is on, output is in SBML format
        """
        name = ReactionName(self.name) if SBML else self.name
        string_ = name + " : "
        
        for react in self.reactants:
            stoich  = '' if react.stoichiometry == 1 else str(react.stoichiometry) + '*'
            if SBML:
                string_ = string_ + stoich + MetaboliteName(react.originalMet.name) + ' + '
            else:
                string_ = string_ + stoich + react.originalMet.name + ' + '                
        string_ = string_[0:-3]
        
        # irrev/rev symbol
        if self.reversible:
            symbol = ' <==> '
        else:
            symbol = ' --> '
        string_ = string_ + symbol
        
        for prod in self.products:
            stoich  = '' if prod.stoichiometry == 1 else str(prod.stoichiometry) + '*'
            if SBML:
                string_ = string_ + stoich + MetaboliteName(prod.originalMet.name) + ' + '                
            else:
                string_ = string_ + stoich + prod.originalMet.name + ' + '
        string_ = string_[0:-3]   
    
        return string_


    def getReactDict(self):
        """
        Produces dictionary of reactants
        """
        reactDict = {}
        for react in self.reactants:
            reactDict[react.name] = react
        self.reactDict = reactDict
        return reactDict


    def getProdDict(self):
        """
        Produces dictionary of products
        """
        prodDict = {}
        for prod in self.products:
            prodDict[prod.name] = prod
        self.prodDict = prodDict
        return prodDict


    def syncTransitionLine(self,prioritySBML): 
        " Syncs information on transitionLine with information in reaction"
        
        # Check reactants and products from reaction line are in reaction
        reactsSBMLDict = self.reactDict
        prodsSBMLDict  = self.prodDict
        
        reacts     = copy.deepcopy(self.transitionLine.reactants)
        prods      = copy.deepcopy(self.transitionLine.products)
        
        for react in reacts:
            if not react.name in reactsSBMLDict:
                print(reactsSBMLDict)
                raise Exception('Reactant '+react.name+' present in transitions line but not in reaction: '+self.name+':\n'+self.stoichLine())
        for prod in prods:
            if not self.exchange:   # product are not included for exchange reactions in sbml
                if not prod.name in prodsSBMLDict:
                    print(prodsSBMLDict)
                    raise Exception('Product '+prod.name+' present in transitions line but not in reaction: '+self.name+':\n'+self.stoichLine())
        
        # Reversibilities
        if prioritySBML:
            if self.reversible:
                self.transitionLine.conv2reversible()
            else:
                self.transitionLine.conv2irreversible()
        else:
            self.reversible   = self.transitionLine.reversible
        
        # Upper and Lower bounds
        if self.fluxBounds:
            measured,flux = fluxBounds(self.fluxBounds.net.lo,self.fluxBounds.net.hi,self.reversible,self.exchange)
            self.fluxBounds = flux


    def checkTransitionLineCompleteness(self,coremets):
        """
        Checks that transition lines products and reactants are complete. 
        i.e. if they show up in core mets and in stoichiometry they should show up in transitionLine
        """
        reactsSBML = self.reactants
        prodsSBML  = self.products

        reactsSBMLDict = self.reactDict
        prodsSBMLDict  = self.prodDict

        reacts     = copy.deepcopy(self.transitionLine.reactants)
        prods      = copy.deepcopy(self.transitionLine.products)

        reactsDict = reacts.metDict
        prodsDict  = prods.metDict

        for react in reactsSBML:  
            if (react.name in reactsSBMLDict) and (react.name in coremets) and not(react.name in reactsDict): 
                print(react.name)
                print(reactsDict)
                print(self.stoichLine())
                print(self.transitionLine)
                raise Exception('Reactant '+react.name+' present in core metabolites and in stochiometry but not in transitions for reaction: '+self.name)
        for prod in prodsSBML:
            if not self.exchange:   # product are not included for exchange reactions in sbml
                if (prod.name in prodsSBMLDict) and (prod.name in coremets) and not(prod.name in prodsDict):
                    print(prod.name)
                    print(prodsDict)
                    print(self.stoichLine())
                    print(self.transitionLine)
                    raise Exception('Product '+prod.name+' present in core metabolites and in stochiometry but not in transitions for reaction: '+self.name)


    def setGenes(self, geneSets=None):
        """
           When called with an array of GeneSet objects, this function makes a deep copy of each object,
        and assigns the copies to this Reaction, saving them in self.geneSets.
        Each GeneSet contains one or more Gene objects that collectively control this Reaction.
        For example, if self.geneSets contains two GeneSets, A and B, such that:
            A contains Gene objects CCDC83, RPL11, LUC7L3,
            B contains Gene objects DDX26B, LUC7L3,
        then the Reaction is controlled by the logical rule "(CCDC83 AND RPL11 AND LUC7L3) OR (DDX26B AND LUC7L3)"
        This function also creates a single GeneSet object that acts as a master set, and contains a consolidated
        set of all the Gene objects in the objects in self.geneSets.  This is used to de-duplicate newly added sets.
           In addition, since each Gene object is derived a NamedRangedNumber object, each Gene can be given a
        measurement value (typically RPKMs).  Since the Gene objects are re-used between all GeneSets, changing the
        value of one named Gene object will effectively change it in all of them.
           For example, changing the measured value of LUC7L3 in GeneSet A above will also automatically change it
        for LUC7L3 in GeneSet B, and in the master set as self.genes.
        """
        # A consolidated master set containing all Gene objects
        self.genes = Genes.GeneSet()
        # A list of sets of genes, each set a potential cause of the reaction
        self.geneSets = []
        if geneSets is not None:
            # Make sure all the Gene objects are represented in the master set,
            # and that genes mentioned multiple times are represented by the same Gene object.
            for subSet in geneSets:
                self.geneSets.append(self.genes.recastSet(subSet))


    def setProteins(self, proteinSets=None):
        """
           When called with an array of ProteinSet objects, this function makes a deep copy of each object,
        and assigns the copies to this Reaction, saving them in self.proteinSets.
        Each ProteinSet contains one or more Protein objects that collectively control this Reaction.
        For example, if self.proteinSets contains two ProteinSets, A and B, such that:
            A contains Protein objects VEGF, TGFB1, P53,
            B contains Protein objects VMF, P53,
        then the Reaction is controlled by the logical rule "(CCDC83 AND RPL11 AND LUC7L3) OR (DDX26B AND LUC7L3)"
        This function also creates a single ProteinSet object that acts as a master set, and contains a consolidated
        set of all the Protein objects in the objects in self.proteinSets.  This is used to de-duplicate newly added sets.
           In addition, since each Protein object is derived a NamedRangedNumber object, each Protein can be given a
        measurement value.  Since the Protein objects are re-used between all ProteinSets, changing the
        value of one named Protein object will effectively change it in all of them.
           For example, changing the measured value of P53 in ProteinSet A above will also automatically change it
        for P53 in ProteinSet B, and in the master set at self.proteins.
        """
        # A consolidated master set containing all Protein objects
        self.proteins = Proteins.ProteinSet()
        # A list of sets of proteins, each set a potential agent of the reaction
        self.proteinSets = []
        if proteinSets is not None:
            # Make sure all the Protein objects are represented in the master set,
            # and that proteins mentioned multiple times are represented by the same Protein object.
            for subSet in proteinSets:
                self.proteinSets.append(self.proteins.recastSet(subSet))


    def getSignificantTranscriptomicsValue(self):
        """
        Walk through all the measurements (transcription values) in the GeneSet objects and test for significance.
        If all the Genes in any GeneSet are at or above the threshold, return true.
        """
        for geneSet in self.geneSets:
            if geneSet.testSignificance(TRANSCRIPTOMICS_SIGNIFICANCE_THRESHOLD):
                return True
        return False


    def getSignificantProteomicsValue(self):
        """
        Walk through all the measurements in the ProteinSet objects and test for significance
        If all the Proteins in any ProteinSet are above the threshold, return true.
        """
        for proteinSet in self.proteinSets:
            if proteinSet.testSignificance(PROTEOMICS_SIGNIFICANCE_THRESHOLD, equalIsntEnough=True):
                return True
        return False



class MetaboliteName(str):
    """Subclass of str for metabolite names, that converts to BiGG Universal format, and provides other formats.
    For example, m = MetaboliteName("ala-L[c]") will show as "ala_L_c" when printed, and
        m.std will contain "ala-L[c]",
        m.sbml will contain "ala_L_c" (same as when printed),
        m.sbmlb will contain "ala_DASH_L_c",
        m.compact will contain "alaL",
        m.no_compartment will contain "ala_L",
    and m.compartment will contain "c".
    Note that some conversions will render nonsensical results for some inputs, for example
    m = MetaboliteName("12ppd__S_e") is a BiGG Universal name and will print as such,
    but m.std will then contain the nonsensical "12ppd--S[e]". """
    
    def __new__(cls, name):
        # Find out what type of name it is and move to SBML format
        defaultCompartment = 'c'

        endingA = name[-3:]
        endingB = name[-2:]
        if (endingA[0]=='[' and endingA[2]==']') and endingA[1] in compartments:          # standard name (ala-L[c])
            newName   = re.sub('\[(\w)\]','_\g<1>',name).replace('-','_')
        elif (endingB[0]=='_') and endingB[1] in compartments and '_DASH_' not in name:   # SBML (ala_L_c)
            newName = name.replace('-','_')
        elif (endingB[0]=='_') and endingB[1] in compartments and '_DASH_' in name:       # SBMLb (ala_DASH_L_c)
            newName   = name.replace('_DASH_','_')  # Going directly to underscores, for SBML name.
        else:                                                                             # Compact (alaL)
            newName = name
            if (name[-1]=='L' or name[-1]=='D') and (name[-2] != '-'):
                newName = name[0:-1]+'_'+name[-1]
            newName = newName.replace('-L','_L').replace('-D','_D').replace('-','_')          
            newName = newName + '_'+str(defaultCompartment)

        return str.__new__(cls, newName)


    def __init__(self, name):
        # Standard name, stored in .std
        # TODO: This can be improved. It turns BiGG's metabolite "12ppd__S_e" into "12ppd--S[e]", for example,
        # which is a nonsensical name.
        # http://bigg.ucsd.edu/models/iY75_1357/metabolites/12ppd__S_e
        self.compartment = name[-1]
        nameNoComp = name[0:-2]
        self.no_compartment = nameNoComp
        self.std   = nameNoComp.replace('_','-')+'['+self.compartment+']'
        # SBML name, a.k.a. BiGG Universal Identifier, stored in .sbml
        self.sbml = name
        # SBMLb name, used in some old SBML documents
        self.sbmlb  = name.replace('_L','_DASH_L').replace('_D','_DASH_D')
        # Compact name, stored in .compact
        newName = name[0:-2]
        self.compact = newName[0:-2] + newName[-2:].replace('_L','L').replace('_D','D')


    def changeCompartment(self,newComp):
        "Returns another MetaboliteName instance with different compartment information"
        return MetaboliteName(self.sbml[0:-1] + newComp)



class ReactionName(str):
    """Class for reaction names, with methods to change the name to different formats"""

    def __new__(cls,name):
        # Move name to SBML format (EX_glc_e_)
        newName   = name.replace('_LPAREN_e_RPAREN_','(e)').replace('_LPAREN(e)RPAREN_','(e)').replace('(e)','_e_').replace('_DASH_','_') 

        return newName

    def getNameIn(self, nameType):
        "Produces name in desired format: std (standard, EX_glc(e)) or SBML (EX_glc_e_)"
        
        if nameType == 'std':
            newName = self.name.replace('_LPAREN_e_RPAREN_','(e)').replace('_LPAREN(e)RPAREN_','(e)').replace('_e_','(e)').replace('_DASH_','-') 
        elif nameType == 'SBML':
            newName = self.name
        else:
            raise Exception('nameType unknown: '+str(nameType))

        return newName



class EMU(object): # TODO(Hector): redo, holding the Metabolite instance?
    """
    Class for Elementary Metabolite Units (EMUs) as defined in 
    Antoniewicz MR, Kelleher JK, Stephanopoulos G: Elementary metabolite units (EMU): a novel framework for modeling isotopic distributions.
    Metab Eng 2007, 9:68-86.
    """

    def __init__(self,name,equivalent='default'):
        """
        name is the EMU name (e.g. cit_1_2_3_4_5_6)
        equivalente is an equivalent name
        """
        self.name     = name
        self.getMetName()
        if equivalent != 'default':
            self.equivalent = equivalent
            
        # Finding if it is a dummy EMU
        emuNames = re.findall('([\w_]+?(?:_\d\d?)+)',name)
        self.dummy = len(emuNames) > 1

        # Finding number of carbons
        self.ncarbons = self.findnCarbons()


    def findnCarbons(self):
        """
        Finds the number of carbons for the EMU from the name, e.g. cit_1_2_3_4_5_6 --> 6
        """
        if self.dummy:
            emuNames = re.findall('([\w_]+?(?:_\d\d?)+)',self.name)
            ncarbs   = 0
            for name in emuNames:
                emuInst = EMU(name)
                ncarbs  = ncarbs + emuInst.findnCarbons()
        else:
            carbonInds    = re.sub(self.met,'',self.name)
            ncarbs = float(carbonInds.count('_'))    
        
        return ncarbs


    def getMetName(self):
        """
        Finds name of metabolite related to EMU: e.g.  cit_1_2_3_4_5_6 --> cit
        """
        newname   = re.sub('(_\d\d?)+','',self.name)
            
        self.met  = newname
        return newname


    def getIndices(self):
        """
        Obtains indices for EMU: e.g.  cit_1_2_3_4_5_6 --> [1,2,3,4,5,6]
        """
        indList = re.findall('(_\d\d?)', self.name)
        
        inds = []
        for ind in indList:
            inds.append(int(ind.replace('_','')))
            
        return inds


    def getSortedName(self):
        """ 
        Provides names with carbon numbers sorted, 
        e.g.: dhap_c_1_2_3 instead of dhap_c_3_2_1
        """
        name = self.getMetName()
        inds = self.getIndices()
        
        # Sort indices
        indsInt =[]
        for ind in inds:
            indsInt.append(int(ind))
        indsInt = sorted(indsInt)
        indsFin =[]
        for ind in indsInt:
            indsFin.append(str(ind))
        
        sortedName = name+'_'+'_'.join(indsFin)
        
        return sortedName


    def getEmuInSBML(self):
        """
        Provides EMU name in SBML format 
        """
        name    = self.name
        metname = self.met
        ending  = name.replace(metname,'')
        newname = MetaboliteName(metname)+ending
        self.name = newname
        self.getMetName()
        return newname  



##################################
# ATOM AND EMU TRANSITIONS
##################################

class EMUTransition(object):
    """
    Class for EMU transitions that contain information on how different EMUs transform intto each other. For example:
    
        TA1_b,  TAC3_c_1_2_3 + g3p_c_1_2_3 -->  f6p_c_1_2_3_4_5_6
        
    indicating that TAC3_c_1_2_3 and g3p_c_1_2_3 combine to produce f6p_c_1_2_3_4_5_6 in transition TA1_b (backward transition of TA1), or:
        
        SSALy, (0.5) sucsal_c_4 --> (0.5) succ_c_4
        
    which indicates that the fourth atom of sucsal_c becomes the fourth atom of succ_c. The (0.5) contribution coefficient indicates that
    transition SSALy contains a symmetric molecule and two labeling correspondences are equally likely. Hence this transition only 
    contributes half the flux to the final labeling. An alternative input for these types of EMU transitions is in the form of a tuple:
    
        (0.5, 'SSALy, sucsal_c_4 --> succ_c_4')
    
    """
    
    def __init__(self,input_):
        
        if type(input_)==tuple:
            contribution, line = input_
            self.__init__(line)
            self.contribution = contribution
            
        elif type(input_)==str:   
            line = input_
            self.contribution = ''
            "Initial parsing"
            rxnSign = '-->'
            # Getting reaction name:
            rxnName,rest = line.split(',')
            self.name = rxnName
            # Getting contribution coefficients, if they exist
            parse = re.findall('\(([\.\d]+)\)',rest)
            if parse:
                # Testing that all contribution coefficients are the same
                first = float(parse[0])
                for number in parse:     
                    if float(number)!=first:
                        print('All contribution coefficients must be the same: \n'+ line)
                self.contribution = first
                # Eliminating contribution coefficients    
                rest = rest.replace('('+str(first)+')','')                
            # Obtaining reactant and product EMUs
            reacts,prod = rest.split(rxnSign)
            # Obtaining reactant list (product will always be, by construction, a single product)
            reactList = reacts.split('+')
        
            # Storing products and reactants
            self.products = [EMU(prod.strip())]
            reactants = []
            for react in reactList:
                reactants.append(EMU(react.strip()))
            self.reactants = reactants
            
            # Record whether it is a condensing or non-condensing transition
            self.cond = True if len(self.reactants)>1 else False
        else:
            raise Exception("Data type: "+str(type(input_))+" erroneous for input: "+str(input_))


    def __str__(self):
        return self.printLine()


    def printLine(self, sort = True):
        "EMU transition output in line form."
        name      = self.name      
        reactants = self.reactants 
        products  = self.products  

        # Get contribution if present
        cont = '('+str(self.contribution)+') ' if (self.contribution and self.contribution !=1) else ''
        
        if sort:
            reacLine = ' + '.join([react.name for react in sorted(reactants, key=lambda react: react.name)])
        else:
            reacLine = ' + '.join([react.name for react in reactants])
        line = name+', '+cont+ reacLine +' --> '+cont+products[0].name       
        
        return line



class AtomTransition(object):  
    """
    Class for line with transitions info, e.g.: 
    
    ILETA	   akg + ile-L[c] <==> glu-L[c] + 3mop	ABCDE + abcdef : ABCDE + abcdef
    ADSL2r   25aics[c] --> aicar[c] + fum[c]	ABCDEFGHIcdba : ABCDEFGHI + (abcd;bcda)
    PTRCTA   ptrc[c] + akg[c] --> 4abutn[c] + glu-L[c]	(bcde;edcb) + ABCDE : bcde + ABCDE

    Different 'domains' (name, metabolites and labeling info) separated by tabs ('\t')
    """
    def __init__(self,line):
        # Store line for future use
        self.fullLine = line
        
        # Parsing line and obtaining alternative line 
        name, reactLine, prodLine, rLabelLine, pLabelLine, reversible, symmetric = self.parseLine(line)  # TO DO: change this to parseTransition!!
        reacts = [x.strip() for x in reactLine.split('+')]
        prods  = [x.strip() for x in prodLine.split('+')]
        rLabel = [x.strip() for x in rLabelLine.split('+')]
        pLabel = [x.strip() for x in pLabelLine.split('+')]
      
        self.name = name
        self.reversible = reversible
        # Obtaining alternative line (see createAltLine method for explanation)
        # TODO(Hector): make this global, so fum_ps1 from r1 does not get mixed with fum_ps1 from r2 ?
        self.altLine = self.createAltLine(line)

        #### Creating products and reactants #####    
        reactants = AtomTransition.makeListOfMets(reacts,rLabel,mType='reactants')
        products  = AtomTransition.makeListOfMets(prods ,pLabel,mType='products')
        metabolites = [reactant.originalMet for reactant in reactants]
        metabolites.extend([product.originalMet for product in products])
               
        # Test if there are repeated metabolites, e.g.: accoa + accoa --> fum
        repMet = False
        for mList in [reactants,products]:
            for met in mList:
                if met.stoichiometry>1:
                    repMet = True
        
        # We can't allow for the moment atom transitions with both repeated metabolites and symmetric molecules
        if symmetric and repMet:
            raise Exception('Repeated metabolites and symmetric molecules not allowed at the same time: '+str(self.fullLine))
        
        # Final storage        
        self.reactants     = enhancedLists.MetaboliteList(reactants)
        self.products      = enhancedLists.MetaboliteList(products)
        self.metabolites   = enhancedLists.MetaboliteList(metabolites)
        self.sym           = symmetric
        self.repMet        = repMet


    def __str__(self):
        return self.printLine()  


    def parseLine(self, line):
        """
        Parses original line:
            ILETA	   akg + ile-L[c] <==> glu-L[c] + 3mop	ABCDE + abcdef : ABCDE + abcdef
        
        into constituent parts:
            - name          : ILETA
            - reactLine     : akg + ile-L[c]    (line for reactants)
            - prodLine      : glu-L[c] + 3mop   (line for products)
            - rLabelLine    : ABCDE + abcdef    (labeling info for reactants)
            - pLabelLine    : ABCDE + abcdef    (labeling info for products)
            - reversible    : True              (is reaction reversible?)
            - symmetric     : False             (does reaction involve symmetric metabolites?)
        """
        #### Basic form  ####

        splitLine = line.split('\t')
        if len(splitLine) != 3:
            raise Exception('Line must comprise three dominions separated by tabs:\n' + str(line))

        # Reaction Name

        name = splitLine[0]
        metLine = splitLine[1]
        labelLine = splitLine[2]

        # Metabolites

        # Reaction symbol and reversibility
        rxnSymbols = ['-->','<==>']   # no need to make this global yet: this is the only place it is used
        # Testing if any of the allowed reaction symbols is present
        test = [symb in metLine for symb in rxnSymbols]
        if sum(test) != 1:
            raise Exception('Line does not have a known reaction symbol (acceptable symbols are '+str(rxnSymbols)+'):\n'+metLine)
        # Checking reversibility
        rxnSymb = rxnSymbols[test.index(True)]
        reversible = False
        if rxnSymb == '<==>':
            reversible = True

        # Reactants and products

        halves = metLine.split(rxnSymb)
        if len(halves) != 2:
            raise Exception('Line must have at least one product and one reactant around the reaction symbol:\n'+metLine)
        reactLine = halves[0]
        prodLine  = halves[1]

        # Labeling

        # Obtaining reaction symbol for labeling part
        rxnSymbolsPlus = rxnSymbols
        rxnSymbolsPlus.append(':')   #labeling part reaction symbol can include : as well
        test = [symb in labelLine for symb in rxnSymbolsPlus]
        if sum(test) != 1:
            raise Exception('Line labeling does not have a known reaction symbol: '+str(rxnSymbolsPlus)+'\n'+labelLine)
        rxnSymb = rxnSymbolsPlus[test.index(True)]            
        # Reactants and products labeling
        rLabelLine, pLabelLine =  labelLine.split(rxnSymb)
        symmetric = ('(' in labelLine) and (';' in labelLine) and (')' in labelLine)   # Testing symmetric molecule labeling (e.g. (abc;cba) 
        
        return name, reactLine, prodLine, rLabelLine, pLabelLine, reversible, symmetric


    def createAltLine(self,line):
        """
        In order to take into account the effects of symmetric molecules and repeated metabolites in a unified manner, an alternative line is
        created that has unique names and labelings for each metabolites. 
        
        For example (symmetric molecule):

        r4  akg_c --> suc_c + co2_c abcde : (bcde;edcb) + a            (suc_c is a symmetric molecule)
        
        becomes: 
        
        r4  akg_c + akg_c__ps1 --> suc_c + suc_c__ps1 + co2_c + co2_c__ps1   abcde + fghij : bcde + jihg + a + f

        Where suc_c has been unfolded into suc_c and suc_c__ps1 (pseudo metabolite) 
        
        Symilarly, for repeated metabolites:
        
        GLUSx   akg_c + gln_L_c --> glu_L_c + glu_L_c   abcde + fghij : abcde + fghij
        
        becomes:
        
        GLUSx   akg_c + gln_L_c --> glu_L_c + glu_L_c__ps1  abcde + fghij : abcde + fghij
        
        
        This line can be easily used to do the EMU decomposition
        """
        # Parse line
        name, reactLine, prodLine, rLabelLine, pLabelLine, reversible, sym = self.parseLine(line)
        reacts = [x.strip() for x in reactLine.split('+')]
        prods  = [x.strip() for x in prodLine.split('+')]
        rLabel = [utils.symsplit(x.strip()) if ';' in x else x.strip() for x in rLabelLine.split('+')]
        pLabel = [utils.symsplit(x.strip()) if ';' in x else x.strip() for x in pLabelLine.split('+')]        
        
        # If there is a symmetric molecule
        if sym:
            ## Duplicate metabolites
            newReacts = []
            for react in reacts:
                newReacts.append(react)
                newReacts.append(react)
            newProds = []
            for prod in prods:
                newProds.append(prod)
                newProds.append(prod)
            newReacts = AtomTransition.convert2PS(newReacts)
            newProds  = AtomTransition.convert2PS(newProds)
            
            ## Take care of labeling
            # Find out all carbon letters (= inputKeys)
            inputKeys=[]
            for lab in rLabel:
                if type(lab) == list:
                    for comp in lab:
                        inputKeys.extend([x for x in comp])
                else:
                    inputKeys.extend([x for x in lab])
            inputKeys = sorted(set(inputKeys))
            
            # Find equivalent for each of the old carbon letters 
            all_letters = string.ascii_letters   # All possible letters
            avoid = ''.join(inputKeys)           # These are taken letters, to be avoided in choosing new equivalents
            keyDict = {}
            for key in inputKeys:
                keyDict[key]  = [x for x in all_letters if x not in avoid][0]
                avoid = avoid + keyDict[key]
        
            # Add new labeling for new metabolites        
            newRlabel = []
            for lab in rLabel:
                if type(lab) is list:    # Only works for two alternative labelings !!!
                    newRlabel.append(lab[0])
                    newRlabel.append(''.join([keyDict[x] for x in lab[1]]))
                else:
                    # Creating new labeling for new metabolites
                    newRlabel.append(lab)
                    newLab = ''.join([keyDict[x] for x in lab])
                    newRlabel.append(newLab)
            newPlabel = []
            for lab in pLabel:
                if type(lab) is list:
                    newPlabel.append(lab[0])
                    newPlabel.append(''.join([keyDict[x] for x in lab[1]]))
                else:
                    # Creating new labeling for new metabolites
                    newPlabel.append(lab)
                    newLab = ''.join([keyDict[x] for x in lab])  
                    newPlabel.append(newLab)
        else:
            # Convert reactants and products to pseudo metabolits (e.g. glu_L_c --> glu_L_c__ps1) if needed
            newReacts = AtomTransition.convert2PS(reacts)
            newProds  = AtomTransition.convert2PS(prods)
            
            newRlabel = rLabel
            newPlabel = pLabel
            
        ## Join all into new line
        rxnSymb = ' <==> ' if reversible else ' --> '
        altLine =  name 
        altLine += ' \t '+' + '.join(newReacts)+rxnSymb+' + '.join(newProds)
        altLine += ' \t '+' + '.join(newRlabel)+' : '+' + '.join(newPlabel)
        
        return altLine


    @staticmethod
    def convert2PS(metNames):
        """
        Converts a list of metabolite names to a list of nonrepeated names with pseudo metabolites, e.g.:
        ['accoa','accoa','fum']  --> ['accoa','accoa__ps1','fum']
        """
        pseudoEnd ='__ps'
        metNamesPS=[]
        metTimes = {}
        for met in metNames:
            if metNames.count(met)>1:
                if met not in metTimes:
                    metTimes[met] =  0
                else:
                    metTimes[met] += 1
                post   = '' if metTimes[met] == 0 else pseudoEnd+str(metTimes[met])
                newMet = met+post
            else:
                newMet = met
                
            metNamesPS.append(newMet)
        
        return metNamesPS


    @staticmethod
    def makeListOfMets(metNames,metLabel,mType='reactants'):
        """
        Function that takes a list of names (e.g.: ['fum', 'fum', 'ac']) and labels (e.g.: ['abcd', 'ABCD', '(ef;fe)'])
        and converts them into an instance of the reactant (or product, depending on mType) class
        """
        
        # Clean input
        metNames = [name.strip() for name in metNames]
        metLabel = [name.strip() for name in metLabel]
        
        # Main algorithm
        metabolites = []    # final output
        doneMets    = []    # metabolites taken care of
        for name in metNames:
            # Obtain name, stoichiometry and labeling patterns
            act = True
            if metNames.count(name)>1: 
                if name not in doneMets:
                    # Find indices
                    indices = [i for i, x in enumerate(metNames) if x == name]
                    # Find associated labeling
                    labels = [metLabel[i] for i in indices] # TODO: take care of symmetric molecules
                    # Check that there are no symmetric molecules
                    # TODO: extend this to take care of symmetric molecules in duplicated metabolites
                    for label in labels:
                        if ('(' in label) or (';' in label) or (')' in label):
                            outStr = ' + '.join(metNames)+' --> '+' + '.join(metLabel)
                            raise Exception('Duplicated metabolites cannot be symmetric: '+str(outStr))
                    # Add to list of done metabolites to avoid repetitions and get stochiometry from repetition number
                    doneMets.append(name)
                    stoichiometry = metNames.count(name)  
                else:
                    act = False       # act on each repeated metabolite only once
            else:
                # Get labeling and test for symmetric molecules
                label = metLabel[metNames.index(name)].strip()
                labels = utils.symsplit(label) if label[0]=='(' else [label]
                
                stoichiometry =1 

            # Introduce name, stoichiometry and labeling patterns in reactant or product structure
            if act:   # act true only for first instance of repetition
                metab    = Metabolite(name, ncarbons = len(labels[0]))
                if  mType == 'reactants':
                    metabOut    = Reactant(metab, stoichiometry, label = labels[0],extraLabel = labels[1:])  
                elif mType == 'products':
                    metabOut    = Product(metab, stoichiometry, label = labels[0],extraLabel = labels[1:])   
                else:
                    raise Exception('Type must be either "reactants" or "products" :' +str(mType))
            
                metabolites.append(metabOut)
            
        return metabolites


    def printReverse(self, altLine = False):  
        """
        Prints the line in reverse, e.g.:
        
        r4  akg_c --> suc_c + co2_c abcde : (bcde;edcb) + a
        
        becomes:
        
        r4  suc_c + co2_c<==>akg_c  (bcde;edcb) + a : abcde 
        
        """       
        # Parse line
        name, reactLine, prodLine, rLabelLine, pLabelLine, reversible, symmetric = self.parseLine(self.printLine(altLine = altLine))
        
        rxnSymb = '<==>' if self.sym else '-->'
        # Inverse reactants and products
        newLine = name+'\t'+prodLine+rxnSymb+reactLine+'\t'+pLabelLine+' : '+rLabelLine
        
        return newLine


    def printLine(self, altLine = False):
        """
        Prints line
        """
        #TODO(Hector): coordinate with __str__ and printReverse
        if not altLine:
            output = self.fullLine
        else:
            output = self.altLine
        
        return output  
        

    def getTransIn(self,outputType):
        """Provides the transition in the desired output format (e.g. SBML)"""
        if outputType == 'SBML':
            name,reac,trans = self.fullLine.split('\t')
            name = self.name
            bits = reac.split(' ')
            newbits = []
            for bit in bits:
                if not(bit == '+' or bit =='-->' or bit == '<==>' or bit == ''):
                    newbits.append(MetaboliteName(bit))
                else:
                    newbits.append(bit)                     
            newline = ReactionName(name)+"\t"+" ".join(newbits)+"\t"+trans
            self.__init__(newline)
            
            output = newline        
        else:
            raise Exception('Conversion string not recognized')
        return output


    def conv2reversible(self):
        """Converts reaction to reversible"""
        newLine = self.fullLine
        newLine = newLine.replace('-->','<==>')
        self = self.__init__(newLine)        


    def conv2irreversible(self):
        """Converts reaction to irreversible"""
        newLine = self.fullLine
        newLine = newLine.replace('<==>','-->')
        self = self.__init__(newLine)        


    def getOriginDictionary(self,emu):
        """
        Returns dictionary with names of reactants involved in the creation of EMU
        Dictionary entries are the carbons indices involved in the match in tuple form, e.g.:
            
            pyr: ( [1,2,3] , [4,5,6] )
            
        This means carbons 1, 2 and 3 in pyr correspond to carbons 4,5,6 in EMU metabolite.   
        Extra labeling patterns due to symmetric molecules are ignored because they are taken care of in getTransitionDict
        In the special case that the same metabolite maps to different parts of the EMU, for example for reaction:
        
            ACLS	pyr_c + pyr_c --> alac_S_c + co2_c	cde + fgh : fgdhe + c
            
        the dictionary includes a new entry with a pseudo name, for example (emu = alac_S_c_1_2_3_4_5 ):
                                
            pyr:       ( [2,3] , [3,5] ) 
            pyr__ps1:  ( [1,2,3] , [1,2,4] ) 
            
        """
        # Find name for metabolite corresponding to EMU 
        metName   = emu.met  # This is a string
        
        # Find all carbon mappings
        # Find product carbon string
        inds = emu.getIndices()
        altTrans = AtomTransition(self.altLine)            
        for prod in altTrans.products:    # use alternative transitions where there are no repeated metabolites
            if prod.name == metName:
                prodLabel = prod.label
               
        # Find reactants and indices    
        originDict = {}
        for react in altTrans.reactants:
            mapped, indMap = utils.indexMap(inds,prodLabel,react.label)
            if mapped:
                originDict[react.name]=indMap
        
        return originDict


    def findEMUs(self,emu):
        """
        Finds EMUs involved in the creation of EMU
        For example, for this carbon transition:
        
        ILETA	   akg + ile-L[c] <==> glu-L[c] + 3mop	ABCDE + abcdef : ABCDE + abcdef        
        
        the result is:
        
        trans.findEMUs(core.EMU('3mop_1_2_3_4_5')) = ile-L[c]_1_2_3_4_5      
        
        """
        
        # Check existence of pseudo metabolites
        originDicts = []
        altTrans = AtomTransition(self.altLine)            
        for prod in altTrans.products:
            if re.sub('(__ps\d+)','',prod.name) == emu.met:
                emuName = emu.name.replace(emu.met,prod.name)
                emuNew  = EMU(emuName)                
                originDict=self.getOriginDictionary(emuNew)
                originDicts.append(originDict) 
            
        # Output is a list of lists with EMUs for each pseudo metabolite
        output = []
        for originDict in originDicts:
            outList = []
            for name in originDict:
                cleanName = re.sub('(__ps\d+)','',name)      # Name is cleaned of pseudo name additions (see getOriginDictionary)
                emuName = cleanName+'_'+'_'.join(str(x) for x in originDict[name][1])
                emuOut = EMU(emuName)
                outList.append(emuOut)
            output.append(outList)
            
        return output


    def findEMUtransition(self,emu):
        """
        Finds EMU transition(s) involved in creation of EMU
        For example, for this carbon transition: 
        
        ILETA	   akg + ile-L[c] <==> glu-L[c] + 3mop	ABCDE + abcdef : ABCDE + abcdef        

        the result is:
        
        trans.findEMUtransition(core.EMU('3mop_1_2_3_4_5')) = ['ILETA, ile-L[c]_1_2_3_4_5 --> 3mop_1_2_3_4_5']
        
        """
        # Find EMUs involved in the creation of EMU
        emuList = self.findEMUs(emu)

        # Main algorithm
        output = []
        for emus in emuList:
            out = ''
            if emus:
                nameBit   = self.name.strip()+', '                               # create name of transition (e.g. ILETA,)
                sortedEMUs = sorted(emus, key = lambda emu: emu.name)            # sort EMUs by name for standarized output
                reactsBit  = ' + '.join(x.getSortedName() for x in sortedEMUs )  # bit with reactants
                        
                out = nameBit + reactsBit + ' --> ' + emu.getSortedName()              # put it all together
                out = re.sub('(__ps\d+)','',out)                                 # eliminate all traces of pseudometabolites
                
            output.append(out)
        
        return output     




class rangedNumber(object):
    """
    A floating-point ranged number, with the value in 'best', and the range defined by 'lo' and 'hi' attributes.
    The "best" attribute needs to be a number, but "lo" and "hi" can be either a number or "None".
    (When parsing from a string, "NULL", "NA", "N.A.", "N/A", and the empty string are also acceptable,
    but will be replaced with "None".  In fact, anything that doesn't legally parse as a float, turns into "None".)
    """

    def __init__(self, lo, best, hi):
        assert utils.is_float(best), '"best" value must be interpretable as a number'
        self.best = float(best)
        if utils.is_float(lo):
            self.lo = float(lo)
        else:
            self.lo = None
        if utils.is_float(hi):
            self.hi = float(hi)
        else:
            self.hi = None


    @classmethod
    def from_number(cls, value):
        """
        Construct a new rangedNumber based on a number supplied as an int or float.
        """
        return cls(None, value, None)


    @classmethod
    def from_string(cls, value):
        """
        Construct a new rangedNumber based on a string supplied in the same format as "__str__"
        Format is [low:best:high] where low, best, and high are all integer or floating-point numbers.
        if the string does not match the format it, it will be treated as a string representation of a float or int.
        """
        if rangedNumber.looks_like_rangedNumber(value):
            parts = re.split("\s*:\s*", value)
            lo = parts[0].lstrip("[ ")
            best = parts[1]
            hi = parts[2].rstrip("] ")
            # The regular instantiator will do the rest of the type checking.
            return cls(lo, best, hi)
        else:
            return cls(None, value, None)


    @staticmethod
    def looks_like_rangedNumber(s):
        """
        Test whether a given object can be parsed into a rangedNumber.
        Expected format is a string "[low:best:high]" where low, best, and high are all integer or floating-point numbers.
        """
        if not isinstance(s, str):
            return False
        parts = re.split("\s*:\s*", s.strip())
        if len(parts) != 3:
            return False
        if not re.search('^\[', parts[0]):
            return False
        if not re.search('\]$', parts[2]):
            return False
        if not utils.is_float(parts[1]):    # Only 'best' must be a float.
            return False
        return True


    def copy(self):
        return rangedNumber(self.lo, self.best, self.hi)


    def __str__(self):
        if self.lo == None and self.hi == None:
            return str(self.best)
        else:
            # Shorten '0.0' to '0' to save a bit of space
            lo = str(self.lo)
            if self.lo == 0:
                lo = '0'
            best = str(self.best)
            if self.best == 0:
                best = '0'
            hi = str(self.hi)
            if self.hi == 0:
                hi = '0'
            return '[%s:%s:%s]' % (lo, best, hi)


    def __add__(self, other):
        "ranged number addition. lo and hi determined through standard error propagation"
        if not isinstance(other, rangedNumber):
            other = rangedNumber.from_number(other)
        final = self.best + other.best
        return self.newWithErrorPropogation(other, final)


    def __sub__(self, other):
        "ranged number subtraction. lo and hi determined through standard error propagation"
        final = self.best - other.best
        return self.newWithErrorPropogation(other, final)


    def newWithErrorPropogation(self, other, final):
        "Helper function for __add__ and __sub__"
        if utils.is_float(self.lo) and utils.is_float(other.lo):
            loDS = abs(self.lo - self.best)     
            loDO = abs(other.lo - other.best)         
            finalLo = final - math.sqrt(loDS**2 + loDO**2)
        elif utils.is_float(self.lo):
            finalLo = final - abs(self.lo - self.best)
        elif utils.is_float(other.lo):
            finalLo = final - abs(other.lo - other.best)
        else:
            finalLo = None

        if utils.is_float(self.hi) and utils.is_float(other.hi):
            hiDS = abs(self.hi - self.best)     
            hiDO = abs(other.hi - other.best)         
            finalHi = final + math.sqrt(hiDS**2 + hiDO**2)
        elif utils.is_float(self.hi):
            finalHi = final + abs(self.hi - self.best)
        elif utils.is_float(other.hi):
            finalHi = final + abs(other.hi - other.best)
        else:
            finalHi = None

        return rangedNumber(finalLo, final, finalHi)


    # TODO: Test whether 'other' is a ranged number, and handle accordingly?
    def __mul__(self, other):
        "ranged number multiplication"
        if utils.is_float(self.lo):
            loNew = self.lo * other # Exception will automatically raise if other is None
        else:
            loNew = None
        if utils.is_float(self.hi):
            hiNew = self.hi * other
        else:
            hiNew = None
        if other >= 0:
            result = rangedNumber(loNew, self.best * other, hiNew)
        else:
            result = rangedNumber(hiNew, self.best * other, loNew)
        return result


    def __rmul__(self, other):
        "ranged number reverse multiplication"
        return self.__mul__(other)


    # TODO: Test whether 'other' is a ranged number, and handle accordingly.
    def __div__(self,other):
        "ranged number division"
        if utils.is_float(self.lo):
            loNew = old_div(self.lo, other) # Exception will automatically raise if other is None
        else:
            loNew = None
        if utils.is_float(self.hi):
            hiNew = old_div(self.hi, other)
        else:
            hiNew = None
        if other >= 0:
            result = rangedNumber(loNew, old_div(self.best, other), hiNew)
        else:
            result = rangedNumber(hiNew, old_div(self.best, other), loNew)
        return result


    def __eq__(self,other):
        """Compare two ranged numbers, or a ranged number with a float.
        (A rangedNumber will be equal to a float only when lo and hi are None.)"""
        if not isinstance(other, rangedNumber):
            other = rangedNumber.from_number(other)
        A = self.lo   == other.lo
        B = self.best == other.best
        C = self.hi   == other.hi
        return A and B and C


    def getComp(self,comp='all'):
        "Returns ranged number component (lo, hi, or best)"
        if comp == 'all':
            # Not sure what would be best here to return self or copy of self.
            return rangedNumber(self.lo, self.best, self.hi)
        else:
            return getattr(self,comp)



class StoichMatrix(object):
    "Class for stoichiometry matrices"
    
    def __init__(self,gamsParameter):
        #TODO(Hector): provide initialization with other than gams parameter
        self.gamsParameter = gamsParameter
        mets = set()
        rxns = set()
        for element in gamsParameter.elements:
            met = element[0]
            rxn = element[1]
            mets.add(met)
            rxns.add(rxn)
        self.mets = list(mets)
        self.rxns = list(rxns)
        

    def getSmatrix(self):
        """
        Produces stoichiometry matrix in array form.
        """
        # S matrix 
        mets = self.mets
        rxns = self.rxns
        metsset = set(mets)
        rxnsset = set(rxns)
        Elements = self.gamsParameter.elements
        Smat = numpy.zeros(shape=(len(metsset),len(rxnsset)))
        for element in Elements:
            # find index for met and rxn
            met = element[0]
            rxn = element[1]

            imet = self.mets.index(met)
            jrxn = self.rxns.index(rxn)

            Smat[imet][jrxn] = Elements[element]

        self.S = Smat    
        return Smat    


    def getExchangeMets(self):
        "Provides exchange metabolites based on stoichiometry info"
        mets = self.mets
        # All this should be done by naming exchange reactions as EX_something. 
        S = self.getSmatrix()
        ExchaMet = []
        for met in mets:
            imet = mets.index(met)
            inFlux  = S[imet]>0
            outFlux = S[imet]<0
            # If metabolite does not have reactions coming in AND reactions coming out
            if not (inFlux.any() and outFlux.any()):
                ExchaMet.append(met)
                
        return ExchaMet
        
    # def write
    # def read



# Flux classes
class flux(object):
    """
    Class for flux information.
    """

    def __init__(self, for_back_tup=None, net_exc_tup=None, net_coeff_tup=None):
        """
        for_back_tup = forward backwards tuple = (forward flux,backward flux)
        net_exc_tup  = net exchange tuple = (net flux, exchange flux)
        inputs can be floats, ranged numbers or strings ('NA' is the only expected one)
        """
        # TODO: Why aren't we converting all these into rangedNumbers internally,
        # upon instantiation?  getLineBest and getComp would become much simpler.
        if for_back_tup != None and net_exc_tup != None:
            forward,backward = for_back_tup
            net,exchange     = net_exc_tup
            forward  = float(forward)  if utils.is_float(forward)  else forward
            backward = float(backward) if utils.is_float(backward) else backward
            net      = float(net)      if utils.is_float(net)      else net
            exchange = float(exchange) if utils.is_float(exchange) else exchange
            
            self.forward  = forward
            self.backward = backward
            self.net      = net
            self.exchange = exchange
            
        elif for_back_tup != None:
            forward,backward = for_back_tup
            forward  = float(forward)  if utils.is_float(forward)  else forward
            backward = float(backward) if utils.is_float(backward) else backward  
            
            self.forward  = forward
            self.backward = backward
            self.net      = forward - backward        
            self.exchange = min(forward,backward)
        
        elif net_exc_tup != None:
            net,exchange  = net_exc_tup
            net      = float(net)      if utils.is_float(net)      else net
            exchange = float(exchange) if utils.is_float(exchange) else exchange

            # get value of net for comparison           
            try:
                netVal = net.best                
            except AttributeError:
                netVal = net                
            
            self.net      = net
            self.exchange = exchange
            if netVal>=0:
                self.forward  = exchange + net
                self.backward = exchange
            else:
                self.forward  = exchange
                self.backward = exchange - net

        elif net_coeff_tup != None:
            net,coeff = net_coeff_tup
            net      = float(net)      if utils.is_float(net)   else net
            coeff    = float(coeff)    if utils.is_float(coeff) else coeff            
        
            exchange  = old_div(coeff,(1-coeff)) 
            
            # get value of net for comparison           
            try:
                netVal = net.best                
            except AttributeError:
                netVal = net                 
                
            self.net      = net
            self.exchange = exchange
            if netVal>=0:
                self.forward  = net + exchange 
                self.backward = exchange
            else:
                self.forward  = exchange
                self.backward = -net + exchange    

        else:
            raise Exception ('Wrong input for flux')


    def __str__(self):
        """
        Provides string description of flux
        """
        line = "Forward: " + str(self.forward) + "\n" + "Backward: " + str(self.backward) + "\n" + "Net: " + str(self.net) + "\n" + "Exchange: " + str(self.exchange) + "\n"
        return line


    def __mul__(self,other):
        "flux multiplication"
        # Multiplying available attributes
        attribs    = ['forward','backward','net','exchange']
        attribDict = {}
        for name in attribs:
            if hasattr(self,name):
                oldAttrib = getattr(self,name)
                if isinstance(oldAttrib, str):     # if attribute is 'NA' leave as it is
                    newAttrib = oldAttrib
                else:
                    newAttrib = oldAttrib*other
                attribDict[name] = newAttrib

        # Creating result
        if len(list(attribDict.keys()))==4:
            result = flux(for_back_tup=(attribDict['forward'],attribDict['backward']),net_exc_tup=(attribDict['net'],attribDict['exchange']))
        else:   
            if 'net' in attribDict and 'exchange' in attribDict:
                result = flux(net_exc_tup=(attribDict['net'],attribDict['exchange']))
            elif 'forward' in attribDict and 'backward' in attribDict:
                result = flux(for_back_tup=(attribDict['forward'],attribDict['backward']))
            
        return result
        

    def __rmul__(self,other):
        "flux reverse multiplication"
        result = self.__mul__(other)
        return result


    def getComp(self,comp='all', rangeComp='all'): 
        """
        This methods returns flux component (net,exchange,forward or backward) and range component (lo, best or hi) as desired, 
        automating and easing a commong need.
        """

        if comp == 'all':
            attributes = ['forward','backward','net','exchange']
            attribDict = {}
            for attrib in attributes:
                initial = getattr(self,attrib)
                if utils.is_float(initial):
                    if rangeComp == 'all':
                        final = rangedNumber(initial,initial,initial)
                    else:
                        final = initial
                elif isinstance(initial, str):
                    final = initial
                elif isinstance(initial, rangedNumber):
                    final = getattr(self,attrib).getComp(rangeComp)
                else:
                    raise Exception('Unknown type for attribute '+attrib)
                attribDict[attrib] = final                        

            output = flux((attribDict['forward'],attribDict['backward']),(attribDict['net'],attribDict['exchange']))    
        else:
            value = getattr(self,comp)
            if isinstance(value, rangedNumber):
                output = value.getComp(rangeComp)
            elif utils.is_float(value):
                if rangeComp == 'all':
                    output = rangedNumber(value,value,value)
                else:
                    output = value
            else:
                if rangeComp == 'all':
                    # If the requested specific component is not a float, but we are not looking for a particular range component,
                    # then we may be expecting the raw value whether or not it's a ranged number (a string 'N/A' for example)
                    output = value
                else:
                    raise Exception('Unsupported type for value: ' + str(value))
                
        return output


    def __div__(self,other):
        "flux division by scalar"
        forward  = self.forward  if isinstance(self.forward, str)  else old_div(self.forward, other)
        backward = self.backward if isinstance(self.backward, str) else old_div(self.backward, other)
        net      = self.net      if isinstance(self.net, str)      else old_div(self.net, other)
        exchange = self.exchange if isinstance(self.exchange, str) else old_div(self.exchange, other)
        
        newflux = flux((forward,backward),(net,exchange))
        return newflux        



def fluxBounds(LB,UB,reversible,exchange=False):
    "Auxiliary function to obtain flux bounds from lower and upper bound in a standardized manner."
    maxVal = 100
    # Reversible reaction        
    measured = False
    if reversible:
        if (LB == -maxVal and UB == maxVal):
            pass
        else:
            measured = True
        if LB >= 0 and UB >= 0:
            fluxOut = flux((rangedNumber(LB,0,UB),rangedNumber(0,0,0)),(rangedNumber(LB,LB,UB),'NA'))                
        elif LB < 0 and UB >= 0:
            fluxOut = flux((rangedNumber(0,0,UB),rangedNumber(0,0,abs(LB))),(rangedNumber(LB,LB,UB),'NA'))                
        elif LB < 0  and UB < 0:
            fluxOut = flux((rangedNumber(0,0,0),rangedNumber(abs(UB),0,abs(LB))),(rangedNumber(LB,LB,UB),'NA'))                                
        else:
            raise Exception("Unacceptable values for lower and upper bounds: "+str(LB)+','+str(UB))
    # Irreversible reaction
    else:             
        if (LB == 0 and UB == maxVal):
            pass
        else:
            measured = True
        # Cap lower bound
        if not exchange:
            LB = max(LB,0)
            UB = max(UB,0)
        fluxOut = flux(('NA','NA'),(rangedNumber(LB,LB,UB),'NA'))

    return measured,fluxOut



######################################################
## AUXILIARY FUNCTIONS
######################################################


##### UNIT TESTS ########

class testMetabolites(unittest.TestCase):  
    """ Testing metabolites classes """

    def testMetabolite(self):
        """ Testing metabolite class """
        
        ala = Metabolite('ala-L', ncarbons=3, source=True, feed='100% 1-C', destination=False, formula='C3H7NO2')
        
        self.assertTrue(str(ala.generateEMU([2]))    == 'ala-L_1_3')
        self.assertTrue(str(ala.generateEMU([2,3]))  == 'ala-L_1')

    def testReactantProduct(self):
        """ Testing reactant, product classes """
        
        ala = Metabolite('ala-L', ncarbons=3, source=True, feed='100% 1-C', destination=False, formula='C3H7NO2')
        R_ala = Reactant(ala, 1, 'abc')        

        self.assertTrue(str(R_ala.generateEMU([2,3]))    == 'ala-L_1')



class testReactions(unittest.TestCase):  
    """ Testing Reaction classes """

    def setUp(self):
           
        # Create reactant metabolites
        coa_c = Metabolite('coa_c')
        nad_c = Metabolite('nad_c')
        pyr_c = Metabolite('pyr_c')

        # Convert into reactants
        Rcoa_c = Reactant(coa_c, 1.0)
        Rnad_c = Reactant(nad_c, 1.0)
        Rpyr_c = Reactant(pyr_c, 1.0)

        # Create product metabolites
        accoa_c = Metabolite('accoa_c')
        co2_c   = Metabolite('co2_c')
        nadh_c  = Metabolite('nadh_c')

        # Convert into products
        Raccoa_c = Product(accoa_c, 1.0)
        Rco2_c   = Product(co2_c, 1.0)
        Rnadh_c  = Product(nadh_c, 1.0)

        # Create reaction
        self.PDH = Reaction('PDH',reactants=[Rcoa_c,Rnad_c,Rpyr_c] , products=[Raccoa_c,Rco2_c,Rnadh_c] 
                    ,subsystem='S_GlycolysisGluconeogenesis')
        self.PDH2 = Reaction.from_string('PDH : coa_c + nad_c + pyr_c --> accoa_c  + co2_c   + nadh_c  ')

    def testPrint(self):
        """ Tests reaction printing """
        
        self.assertTrue(str(self.PDH.stoichLine()) == 'PDH : coa_c + nad_c + pyr_c --> accoa_c + co2_c + nadh_c')
        self.assertTrue(str(self.PDH2.stoichLine()) == 'PDH : coa_c + nad_c + pyr_c --> accoa_c + co2_c + nadh_c')


    def testDicts(self):
        """ Tests reaction printing """
        
        self.assertTrue(str(sorted(self.PDH.getReactDict().keys())) == "['coa_c', 'nad_c', 'pyr_c']")
        self.assertTrue(str(sorted(self.PDH.getProdDict().keys()))  == "['accoa_c', 'co2_c', 'nadh_c']")



class testEmus(unittest.TestCase):  
    """ Testing EMU classes """

    def testEMU(self):
        """ tests EMU class"""
        
        cit321= EMU('cit_3_2_1')

        self.assertTrue(cit321.findnCarbons() == 3)
        self.assertTrue(str(cit321.getMetName()) == 'cit')
        self.assertTrue(str(cit321.getIndices()) == '[3, 2, 1]')
        self.assertTrue(str(cit321.getSortedName()) == 'cit_1_2_3')
        self.assertTrue(str(cit321.getEmuInSBML()) == 'cit_c_3_2_1')

    def testEMUtransition(self):
        """ tests EMUTransition class"""
        
        emuTrans = EMUTransition('TA1_b,  TAC3_c_1_2_3 + g3p_c_1_2_3 -->  f6p_c_1_2_3_4_5_6')        
        self.assertTrue(str(emuTrans) == 'TA1_b, TAC3_c_1_2_3 + g3p_c_1_2_3 --> f6p_c_1_2_3_4_5_6')



class testAtomTransitions(unittest.TestCase):  
    """ Testing AtomTransition class """

    def setUp(self):
        self.AT = AtomTransition('AKGDH	akg --> succoa + co2	abcde : bcde + a')


    def testfindEMUtransition(self):
        emu1 = EMU('co2_1')
        self.assertTrue(str(self.AT.findEMUtransition(emu1)) == "['AKGDH, akg_1 --> co2_1']")
        emu2 = EMU('succoa_1_2_3_4')
        self.assertTrue(str(self.AT.findEMUtransition(emu2)) == "['AKGDH, akg_2_3_4_5 --> succoa_1_2_3_4']")


    def testfindEMUs(self):
        emu2 = EMU('succoa_1_2_3_4')
        self.assertTrue(self.AT.findEMUs(emu2)[0][0].name == 'akg_2_3_4_5')


    def testgetOriginDictionary(self):
        emu2 = EMU('succoa_1_2_3_4')                
        self.assertTrue(str(self.AT.getOriginDictionary(emu2)) == "{'akg': (['b', 'c', 'd', 'e'], [2, 3, 4, 5])}")



class testRangedNumbers(unittest.TestCase):  
    """ Testing ranged number classes """

    def testAll(self):
        
        A = rangedNumber(0.3,0.6,0.9)
        B = rangedNumber(0.1,0.15,0.18)

        self.assertTrue(str(A) == '[0.3:0.6:0.9]')
        self.assertTrue(str(B) == '[0.1:0.15:0.18]')
        self.assertTrue(str(A+B) == '[0.445861873485:0.75:1.05149626863]')
        self.assertTrue(str(A-B) == '[0.145861873485:0.45:0.751496268634]')
        self.assertTrue(str(2*A) == '[0.6:1.2:1.8]')
        self.assertTrue(str(old_div(B,3)) == '[0.0333333333333:0.05:0.06]')



class testFlux(unittest.TestCase):  
    """ Testing flux class """

    def testAll(self):
     
        A = rangedNumber(0.3,0.6,0.9)
        B = rangedNumber(0.1,0.15,0.18)
        
        netFlux      = A
        exchangeFlux = B
        flux1 = flux(net_exc_tup=(netFlux,exchangeFlux)) 

        self.assertTrue(str(flux1.forward) == '[0.445861873485:0.75:1.05149626863]')
        self.assertTrue(str(flux1.backward) == '[0.1:0.15:0.18]')
        self.assertTrue(str(3*flux1.forward) == '[1.33758562046:2.25:3.1544888059]')
        self.assertTrue(str(3*flux1.backward) == '[0.3:0.45:0.54]')
