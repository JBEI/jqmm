# Part of the JBEI Quantitative Metabolic Modeling Library (JQMM)
# Copyright (c) 2016, The Regents of the University of California.
# For licensing details see "license.txt" and "legal.txt".
"""
The ReactionNetworks module holds all data regarding reaction neworks,
i.e. all reactions along with metabolites and how they are interrelated.
Special classes for C13 MFA and 2S-13C MFA are provided
for their particular information requirements.
"""
from __future__ import print_function
from __future__ import division

from builtins import str
from builtins import range
from past.builtins import basestring
from builtins import object
from utilities import old_div
import os, re, numpy, random
import GAMSclasses, core, DB, enhancedLists, labeling, sbmlio
import utilities as utils


class ReactionNetwork(object):
    """ 
        A general class for reaction networks.
        Currently supports SBML file input and output, as well as a tuple defining the major pieces of a ReactionNetwork.
        Future versions will include other input types (e.g. cameo models).
    """
    def __init__(self, source=None, exceptionOnError=True):
        """ input_ can be:
              1-) A SBML document as a string,
              2-) An SBML file, which can be in the form of a filename or a file tuple=(filename,string) where string is the file contents
              3-) An information tuple=(name,Notes,MetaboliteList,ReactionList)
              4-) None, in which case all checks below fall through, and an empty ReactionNetwork is created.
        """
        self.errorCount = 0
        self.logStrings = []

        self.name    = 'undefined'
        self.notes        = {}
        self.metList      = []
        self.reactionList = []

        errorCount = 0
        logStrings = []

        sbmlFileName = None
        sbmlFileContents = None

        if isinstance(source, basestring):
            # If the string has XML brackets, it's likely not a filename.
            # (We don't accept filenames that are complex shell commands, e.g. with STDOUT routings.)
            if b'<' in source and b'>' in source:
                sbmlFileContents = source
            else:
                sbmlFileName = source
        if isinstance(source, tuple):
            if len(source) == 2:
                sbmlFileName, sbmlFileContents = source

        # If we were given a file name, but were not given file content (or were given empty content),
        # we attempt to read in a file and assign it to the content.
        if isinstance(sbmlFileName, basestring) and ((sbmlFileContents is None) or (sbmlFileContents.strip() == '')):
            with open(sbmlFileName, 'r') as file:
                sbmlFileContents = file.read()
            file.closed

        # If any of the above attempts to fill in the file contents succeeded, we parse it as SBML.
        # If not, we do one more test to see if we were given a tuple suitable for unpacking directly into ReactionNetwork parts.
        if isinstance(sbmlFileContents, basestring):
            if sbmlFileContents.strip() != '':
                sbmlImporter = sbmlio.SBMLImporter(sbmlFileContents)

                self.name         = sbmlImporter.model_name
                self.notes        = sbmlImporter.model_notes
                self.metList      = sbmlImporter.metList
                self.reactionList = sbmlImporter.reactionList

                self.errorCount   = sbmlImporter.errorCount
                self.logStrings   = sbmlImporter.logStrings
            else:
                errorCount = 1
                logStrings = ['Given SBML is empty.']
        elif isinstance(source, tuple):
            if len(source) == 4:
                name, notes, metList, reactionList = source
                self.name         = name
                self.notes        = notes
                self.metList      = metList
                self.reactionList = reactionList
            else:
                errorCount = 1
                logStrings = ['Given input tuple is wrong length to contain base ReactionNetwork parts.']
        else:
            errorCount = 1
            logStrings = ['Input does not appear to be SBML, or file does not exist.']

        self.errorCount = self.errorCount + errorCount
        self.logStrings.extend(logStrings)
        if errorCount > 0 and exceptionOnError:
            raise Exception('Errors creating ReactionNetwork:\n%s' % ('\n'.join(self.logStrings)))


    def addReaction(self,reaction):
        """
        Adds a reaction to the network.  
        Input can be a reaction instance or a string of the type:
            
            'HMGCOAR[c] : coa + mev-R + 2 nadp <==> 2 h + hmgcoa + 2 nadph'
        """
        
        # If the reaction is not a reaction object, try to convert it into one.
        if reaction.__class__.__name__ == 'Reaction':
            reactionTmp = reaction
        elif reaction.__class__.__name__ == 'str':
            line = reaction
            reactionTmp = core.Reaction.from_string(line)
        else:
            raise Exception('Input format must be string or reaction instance. Instead got: '+str(reaction))

        # Convert to SBML
        reaction = core.Reaction.from_string(reactionTmp.stoichLine(SBML=True))

        # Find metabolites already in network
        newMetabolites = []
        metList = self.reactionList.getMetList() 
        metabolite_dict = {met.name: met for met in metList}
        
        newReactants=[]                        
        for reactant in reaction.reactants:
            if reactant.name in metabolite_dict:
                reactant = core.Reactant(metabolite_dict[reactant.name], reactant.stoichiometry)
                newReactants.append(reactant)
            else:
                newReactants.append(reactant)
                newMetabolites.append(reactant.originalMet)
        newProducts=[]                        
        for product in reaction.products:
            if product.name in metabolite_dict:
                product = core.Product(metabolite_dict[product.name], product.stoichiometry)
                newProducts.append(product)
            else:
                newProducts.append(product)
                newMetabolites.append(product.originalMet)
        reaction.reactants = newReactants
        reaction.products  = newProducts
        
        # Create new reactionList
        reactionList = self.reactionList
        newReactions = []
        newReactions.append(reaction)
        newReactions.extend(reactionList.reactions)
        self.reactionList = enhancedLists.ReactionList(newReactions)

        # Add metabolites not already present
        if newMetabolites:    
            self.metList.addMetabolites(newMetabolites)
        return True


    def deleteReaction(self, reactionName):
        """
        Eliminates the reaction with the given name from the network.    
        """
        # Create new reactionlist without the reaction corresponding to reactionName
        newReactions = []
        for reaction in self.reactionList.reactions:
            if reaction.name != reactionName:
                newReactions.append(reaction)
                
        self.reactionList = enhancedLists.ReactionList(newReactions)
        self.metList      = self.reactionList.getMetList() 
        
        # Notes and name remain the same


    def changeName(self,name):
        """
        Change the name of the network. This is a trivial function now, but it is provided because it may get more complex in the future.
        """        
        self.name    = name

    # TODO: make this independent of GAMS structures?
    def getStoichMetRxnFBAFiles(self,vSuggested='default'):
        """
        Provides reaction and metabolite files needed for FBA (Flux Balance Analysis)
        """        
        
        fileSet = []

        # Genome-scale stoichiometry matrix
        stoichMatrix = self.reactionList.getStoichMatrix()
        stoichMatrixBigString = stoichMatrix.gamsParameter.write('toString')
        fileSet.append(('Sbig.txt',stoichMatrixBigString))

        # Genome-scale sets of mets and reactions
        hybridMetsSet = GAMSclasses.GAMSSet('hybridmets',set(stoichMatrix.mets))
        hybridReacsSet = GAMSclasses.GAMSSet('hybridreacs',set(stoichMatrix.rxns))
        fileSet.append(('hybridmets.txt',hybridMetsSet.write('toString')))
        fileSet.append(('hybridreacs.txt',hybridReacsSet.write('toString')))

        # Suggested fluxes (i.e. fluxes used for problem initialization)
        if vSuggested == 'default':
            jGEMSUGGset  = GAMSclasses.GAMSSet('jGEMsugg',set([]))    
            vFBASUGGpar  = GAMSclasses.GAMSParameter('vFBASUGG',[])
        else:
            assert vSuggested.__class__.__name__=='fluxGAMSpar', 'vSuggested must be a fluxGAMSpar class'   
            reactionList = vSuggested.getReactionList()
            jGEMSUGGset = GAMSclasses.GAMSSet('jsugg',set(reactionList.getReactionNameList(level=1)))
            vFBASUGGpar = vSuggested            
        fileSet.append( ('jGEMsugg.txt',jGEMSUGGset.write('toString')) )
        fileSet.append( ('vFBASUGG.txt',vFBASUGGpar.write('toString')) )
        
        # Rest of FBA parameters
        reactionList = self.reactionList
        metList      = self.metList

        UBPar   = reactionList.getUB()     # Upper bounds
        LBPar   = reactionList.getLB()     # Lower bounds
        CvecPar = reactionList.getCvec()   # C vector
        BvecPar = metList.getBvec()        # B vector

        # Check that upper and lower bounds (measured fluxes) make sense
        measFluxes = reactionList.getMeasuredFluxes(level=1)
        for name in measFluxes:
            flux = measFluxes[name]
            if flux.lo > flux.hi:
                raise Exception('Flux '+name+' has lower bound higher than upper bound: lb='+str(flux.lo)+',ub='+str(flux.hi))

        fileSet.append( ('ubvec.txt',UBPar.write('toString')) )
        fileSet.append( ('lbvec.txt',LBPar.write('toString')) )
        fileSet.append( ('cvec.txt',CvecPar.write('toString')) )
        fileSet.append( ('bvec.txt',BvecPar.write('toString')) )

        return fileSet        


    def getStoichMetRxnFVAFiles(self,studiedReactions='default'):
        """
        Provides reaction and metabolite files needed for FVA (Flux Variability Analysis)
        """  
        # Most files are the same as FBA
        FVAfiles = self.getStoichMetRxnFBAFiles()

        # Only add files to indicate which reactions are being investigated
        if studiedReactions == 'default':
            stoichMatrix = self.reactionList.getStoichMatrix()
            rxnsSet = GAMSclasses.GAMSSet('studiedfluxes',stoichMatrix.rxns)
            FVAfiles.append( ('studiedfluxes.txt', rxnsSet.write('toString')) )    
        else:
            studiedRxnsSet = GAMSclasses.GAMSSet('studiedfluxes',studiedReactions)
            FVAfiles.append( ('studiedfluxes.txt', studiedRxnsSet.write('toString')) )    

        return FVAfiles


    def getSolverOptFiles(self): 
        """
        Provides optional GAMS files
        """                                              

        conOptFileName   = 'conopt.opt'
        cplexOptFileName = 'cplex.opt'
        try: 
            dirGAMSFile =  os.environ['QUANTMODELPATH']+'/code/core/' 
        except KeyError:
            raise Exception('QUANTMODELPATH environment variable not set')        

        conOptFileName   = dirGAMSFile + conOptFileName
        cplexOptFileName = dirGAMSFile + cplexOptFileName

        return([utils.file2tuple(cplexOptFileName),utils.file2tuple(conOptFileName)])    


    def write(self,fileName):  
        """
        Output reaction network in SBML format
        """  
        if fileName == 'toString':
            return self.getSBMLString()
        else:
            sbmlExporter = sbmlio.SBMLExporter(self)
            return sbmlExporter.writeSBMLToFile(fileName)


    def getSBMLString(self):
        """
        Shortcut call
        """
        sbmlExporter = sbmlio.SBMLExporter(self)
        return sbmlExporter.getSBMLString()


    def changeFluxBounds(self, reactionName, flux, measured=False):
        """
        Changes flux bounds (upper and lower bounds) for reaction, e.g.
        
            reaction.changeFluxBounds('PDH',(0,2))
        
        changes lower bound to 0 and upper bound to 2, and:
        
            reaction.changeFluxBounds('PDH'0)
            
        changes both bounds to zero.
        Input can also be a flux class instance.
        """
        if isinstance(flux, core.flux):        # full flux instantiation provided
            reacDict = self.reactionList.getReactionDictionary()
            reaction = reacDict[reactionName]
            reaction.fluxBounds = flux 
            if measured:
                reaction.measured = True
        elif isinstance(flux, tuple):
            lb,ub = flux
            fluxNew = core.flux(((core.rangedNumber(lb,old_div((lb+ub),2),ub)),(core.rangedNumber(0,0,0))))
            self.changeFluxBounds(reactionName, fluxNew)      
        else:      # shortened version                 
            try:            
                bounds = float(flux)                      
                fluxNew = core.flux(((core.rangedNumber(bounds,bounds,bounds)),(core.rangedNumber(0,0,0))))
                self.changeFluxBounds(reactionName, fluxNew)       
            except:
                print('wrong input!!!')


    def loadFluxBounds(self, fileName, convert2SBML=True, measured=False):
        """
        Loads bounds for fluxes from file. Format of either of these types:
        
            EX_ac(e):  0.3 [==] 0.45
            EX_ac(e)   0.45
            
        convert2SBML (note: defunct) decides whether names are converted to SBML format
        measured option decides whether reactions corresponding to provided flux bounds are labeled as measured
        """
        fileName, string = utils.file2tuple(fileName)    
        lines = string.split('\n')    
        # Load data from file
        bounds = {}
        for line in lines:
            if line.strip():
                match =  re.search('(?P<name>[()_\w]+):?\s+(?P<rest>\S[-\s\d\.=\[\]]+)',line.strip())
                if match:
                    mDict = match.groupdict()
                    # Get reaction name and convert to SBML
                    rxn   = core.Reaction(core.ReactionName(mDict['name'].strip()))
                    # Get upper and lower bounds
                    rest = mDict['rest']
                    if '[==]' in rest:
                        lb    = float(rest.split('[==]')[0])
                        ub    = float(rest.split('[==]')[1])
                        bounds[rxn.name] = (lb,ub)                    
                    else:
                        lb = ub = float(rest.strip())
                else:
                    raise Exception('Flux bound input line: '+line.strip()+'\n not of type "EX_ac(e):  0.3 [==] 0.45" or "EX_ac(e)  0.45"'  )
                # Assign bounds
                bounds[rxn.name] = (lb,ub)
            
        # Put bounds in reactions
        reactionDict = self.reactionList.getReactionDictionary()    
        for name in bounds:
            lb,ub = bounds[name]
            reaction = reactionDict[name]
            self.changeFluxBounds(name, core.fluxBounds(lb, ub, reaction.reversible, reaction.exchange)[1], measured=measured )


    def fixFluxes(self, names='all'):
        "Sets flux bounds (upper and lower) to be the flux values stored in smbl file"
        self.reactionList.fixFluxes(level=2, names=names)


    def capFluxBounds(self, maximum=100):
        "Caps fluxes to a maximum value"
        self.reactionList.capFluxBounds(maximum)


    def getHighestCarbonInput(self, bounds=True):
        """
        Finds name of exchange reaction with highest carbon input.
        This reaction is typically used as a reference for drawing fluxes.
        
            - Bounds indicates whether the flux should be inferred from the lower and upper bounds
            (i.e. flux calculation has not been run) or from the actual flux value
        """
        # Obtain carbon dictionary. 
        carbonDict = utils.getCarbonDict(self.reactionList)  

        minim = 0
        name  = ''
        for reaction in self.reactionList.reactions:
            if reaction.exchange:
                if bounds:
                    ub = reaction.fluxBounds.net.hi
                    lb = reaction.fluxBounds.net.lo
                    if abs(ub-lb) < max(0.1*abs(ub),0.01): 
                        value = -abs(old_div((ub+lb),2)) 
                    else:
                        value = 0
                else:
                    value = reaction.flux.getComp('net','best')
                # Calculate input carbon for the reaction    
                value = value*carbonDict[reaction.reactants[0].name]  
                
                # Find minimum values (maximum input value)
                if value < 0 and value < minim:
                    minim = value
                    name  = reaction.name
        return name



class C13ReactionNetwork(ReactionNetwork):
    """ A derived class for reaction networks containing 13C MFA information """

    def __init__(self, input_):
        debug = False                
                
        # Proceed as for general Reaction Network
        ReactionNetwork.__init__(self, input_)
        if debug:        
            self.reactionList.writeReactionFile('C13lines.txt')  # For debugging purposes
        
        # Process carbon transitions and make necessary checks
        self.reactionList.processCarbonTransitions()        
        # Create carbon dictionary
        self.reactionList.getCarbonDict()


    def addLabeling(self,filename,dataType,STDfilename='default',minSTD=0.001):
        """
        Loads labeling data from files (filename input) of the following type:
            
            Amino acid	Fragment	Mass distribution									
                       m0	      m1	      m2	      m3	     m4	     m5	   m6	       m7	     m8	 m9
            fdp	M-0	0.381 	0.244 	0.081 	0.116	     0.041	     0.017	   0.119	  -	     -	 -
            dhap	0.584 	0.165 	0.094 	0.157     -	     -	   -	       -	     -	 -
            3pg	M-0	0.636 	0.173 	0.038 	0.153     - 	     -	   -	       -	     -	 -

        dataType indicates the kind of data: GCMS, LCMS or CEMS. 
        LCMS and CEMS are equivalent. GCMS varies in naming convention for fragments

        Labeling error come in the same format from STDfileName:
            
           Amino acid	Fragment	Mass distribution									
                       m0	      m1	      m2	      m3	     m4	     m5	   m6	       m7	     m8	 m9
            fdp	M-0	0.021 	0.006 	0.004 	0.016	     0.003	     0.002	   0.009	 -	     -	 -
            dhap	0.007 	0.003 	0.006 	0.002      -	     -	   -	       -	     -	 -
            3pg	M-0	0.005 	0.007 	0.008 	0.005      - 	     -	   -	       -	     -	 -            
            
        minSTD is the minimum error allowable (it needs to be larger than zero to avoid infinities in the objective function).
        
        """ 
        # Check dataType input
        types = ['GCMSLabelData','LCMSLabelData','CEMSLabelData']
        if dataType not in types:
            raise Exception('Type '+dataType+' not among supported types:\n'+str(types))
        
        ## Reading File
        fileName, string = utils.file2tuple(filename)    
        MDVlines = string.split('\n')[2:]    # Leaving headers out
        
        if STDfilename != 'default':   
            STDfilename, string = utils.file2tuple(STDfilename)    
            STDlines = string.split('\n')[2:]    # Leaving headers out
                    
        MSDict = {}
        for line in MDVlines:
            if line:
                STDline = STDlines.pop(0) if STDfilename != 'default' else ''
                MSdata = getattr(labeling, dataType)(lines=(line,STDline),minSTD=minSTD)
                MSDict[MSdata.abbrev] = MSdata
        
        # Modifying Notes from file tuple
        Notes           = self.notes
        Notes[dataType] = MSDict             


    def addFeed(self,filename):
        """
        Adds feed labeling information contained in filename of the type:

        0.4% Glucose: 30% 1-C 20% U 50% UN          (50% unlabeled, 30% labeled in first carbon, 20% full labeled)
        
        """
        fileName, string = utils.file2tuple(filename)        
        
        # Get dictionary of metabolites
        metabolite_dict = DB.metAliases()                     
                     
        # Create feed dictionary
        Feed = {}                     
        for line in string.split('\n'):
            if line:
                inputMet  = re.search(' (\w+):',line).group(1)
                for met in metabolite_dict[inputMet]:
                    Feed[met] = line
            
        # Adding labeling to met information
        mets = self.metList.mets   
        for met in mets:
            if met.name in Feed:
                met.source = True
                met.feed   = Feed[met.name]


    def randomizeLabeling(self):
        """
        Randomizes labeling data patterns within confidence intervals.
        This function is used to find the sensitivity of fluxes to labelin data
        """
        dataTypes = ['GCMSLabelData','CEMSLabelData','LCMSLabelData']
        for type in dataTypes:
            try:
                data = self.notes[type]
                for abbrev in data:
                    data[abbrev].randomize()
            except KeyError:
                pass


    def fragDict(self,verbose=False):
        """
        Produces dictionary of fragments to be fit
        """
        # ## Getting basic dictionary
        fragDictBasic = DB.fragDictBasic()
        
        # ## Creating GCMS, CEMS and LCMS Dictionaries 
        GCMS = self.notes['GCMSLabelData'] if 'GCMSLabelData' in self.notes else []        
        LCMS = self.notes['LCMSLabelData'] if 'LCMSLabelData' in self.notes else []        
        CEMS = self.notes['CEMSLabelData'] if 'CEMSLabelData' in self.notes else []        
        
        GCMSfits = self.notes['GCMSLabelData(fit)'] if 'GCMSLabelData(fit)' in self.notes else []        
        LCMSfits = self.notes['LCMSLabelData(fit)'] if 'LCMSLabelData(fit)' in self.notes else []        
        CEMSfits = self.notes['CEMSLabelData(fit)'] if 'CEMSLabelData(fit)' in self.notes else [] 
        
        # ## Join all data types into single dictionary
#        allData = [GCMS,LCMS,CEMS]

        fragDict = {}
        fragDict.update(GCMS)        
        fragDict.update(LCMS)        
        fragDict.update(CEMS)        

        fragDictFits = {}
        fragDictFits.update(GCMSfits)        
        fragDictFits.update(LCMSfits)        
        fragDictFits.update(CEMSfits)        
        
        # ## Getting maximum MDV length (maximum value of m)
        maxLength = 0
        for abbrev in fragDict:
            length = len(fragDict[abbrev].mdv)
            if length > maxLength:
                maxLength = length
                
        # ## Storing MDV,std and norm in fragDictOut   
        #    and limiting GCMS to the EMUs (elementary metabolite units) present in the data set 
                
        # Obtaining all EMU names in the reaction networks
        # TODO(Hector): create functions in C13RN and TSRN
        try:                             # valid for C13ReactionNetwork
            allEMUs      = self.reactionList.emuList.allEMU.set
        except AttributeError:           # valid for TSReactionNetwork
            allEMUs      = self.C13ReacNet.reactionList.emuList.allEMU.set
            
        # ## Selecting fragments that are being measured
        fragDictOut = {}
        for abbrev in fragDict:
            frag  = fragDict[abbrev]
            frag.std = fragDict[abbrev].std if hasattr(fragDict[abbrev], 'std') else 0.03 * numpy.ones(shape=(maxLength))
            abbrev = abbrev.strip()
            frag.MDVfit = fragDictFits[abbrev].mdv if abbrev in fragDictFits else numpy.array([])   
            if abbrev in fragDictBasic:
                fragRef     = fragDictBasic[abbrev]
                
                present = False
                if fragRef.emu in allEMUs:
                    fragRef.mdv    = frag.mdv
                    fragRef.MDVfit = frag.MDVfit if frag.MDVfit.any() else numpy.array([])
                    fragRef.std    = frag.std
                    fragRef.inFit  = frag.inFit      # This decides whether fragment info is included in the fit or not   
                    fragRef.nmeas  = len(frag.mdv)    
                    fragRef.norm   = numpy.ones(shape=(maxLength)) * len(frag.mdv)
                    fragDictOut[abbrev]= fragRef
                    present     = True                
                    
                if not present and verbose:
                    print("emu " + fragRef.emu + " not present in emu network")
            else:
                print("fragDictBasic:")
                print(sorted(fragDictBasic))
                raise Exception('Fragment ' + abbrev + ' not in fragment dictionary')
        
        return fragDictOut    


    # Output of GCMS data
    # TODO: This can probably be done better through a table class in GAMS
    def getLabelingDataFiles(self):    
        """
        This methods produces all GAMS files related to labeling
        """
        fragDict = self.fragDict()
        
        # keys
        maxLength = 13   # This is the maximum length of the labeling data
        nkeylist  = []
        for i in range(maxLength):
            nkeylist.append(str(i))
        # Values
        MDVs        = []
        MDVerrs     = []
        MDVnorms    = []
        fragAbbrevs = []
        for abbrev in fragDict:
            frag    = fragDict[abbrev]
            if len(frag.mdv) > maxLength:
                raise Exception('Maximum length for MDV is: '+str(maxLength)+' ; fragment '+str(abbrev)+': '+str(frag.mdv))  
            else:
                mdv     = numpy.concatenate((frag.mdv, numpy.zeros(maxLength-len(frag.mdv))))
                std     = numpy.concatenate((frag.std, numpy.zeros(maxLength-len(frag.std))))
                norm    = numpy.concatenate((frag.norm[0:len(frag.mdv)], numpy.zeros(maxLength-len(frag.mdv)))) 

                MDVs.append(mdv)
                MDVerrs.append(std)
                MDVnorms.append(norm) 
                fragAbbrevs.append(abbrev)
            
        # GCMS values
        filenameGCMS = 'GCMSout.txt'    
        name         = 'labelexp'
        xkeys        = GAMSclasses.GAMSTableKeys('n',nkeylist)
        ykeys        = GAMSclasses.GAMSTableKeys('frag',fragAbbrevs)
        tableOut     = GAMSclasses.GAMSTable(name,xkeys,ykeys,numpy.array(MDVs))
        GCMSSt = tableOut.write('toString')

        # GCMS errors
        filenameGCMSstd = 'GCMSerr.txt'    
        name            = 'labelstd'
        xkeys           = GAMSclasses.GAMSTableKeys('n',nkeylist)
        ykeys           = GAMSclasses.GAMSTableKeys('frag',fragAbbrevs)
        tableOut        = GAMSclasses.GAMSTable(name,xkeys,ykeys,numpy.array(MDVerrs))
        GCMSstdSt = tableOut.write('toString')

        # GCMS normalizations
        filenameGCMSnorm = 'GCMSnorm.txt'    
        name             = 'labelnorm'
        xkeys            = GAMSclasses.GAMSTableKeys('n',nkeylist)
        ykeys            = GAMSclasses.GAMSTableKeys('frag',fragAbbrevs)
        tableOut         = GAMSclasses.GAMSTable(name,xkeys,ykeys,numpy.array(MDVnorms))
        GCMSnormSt = tableOut.write('toString')
         
        return [(filenameGCMS,GCMSSt),(filenameGCMSstd,GCMSstdSt),(filenameGCMSnorm,GCMSnormSt)]


    # Get fragment information files
    def getFragmentInfoFiles(self, procString):
        """
        Files on fragment information: fragment names, fragments to be fit, number of carbons, correspondence between EMUs and frag names... etc
        """
        # Unpacking
        fragDictAll = self.fragDict()
        allEMUs = self.reactionList.emuList.allEMU.set

        # Creating fragDict
        fragDict={}
        fitFragNames = []
        for abbrev in fragDictAll:
            frag = fragDictAll[abbrev]
            emuSBML = frag.emu
            if emuSBML in allEMUs:
                fragDict[abbrev]=frag
                if frag.inFit:
                    fitFragNames.append(frag.abbrev)

        if not fitFragNames:
            raise Exception('No fragments found') 

        # Correspondence between EMUs and frag names
        Elements = {}
        for abbrev in fragDict:
            frag    = fragDict[abbrev]
            Elements[tuple([abbrev,frag.emu])] = 1
            
        aacorr = GAMSclasses.GAMSParameter('aacorr',Elements)    

        # number of fragment carbons 
        Elements = {}
        for abbrev in fragDict:
            frag    = fragDict[abbrev]
#            Elements[tuple([abbrev])]=frag.ncarbons
            Elements[tuple([abbrev])]=frag.nmeas-1

        ncarbonsf = GAMSclasses.GAMSParameter('ncarbonsf',Elements)    

        # fragments set
        fragmentsSet = GAMSclasses.GAMSSet('fragments',set(fragDict.keys()))
        
        # Set of fragments involved in the fit
        fitFragsSet = GAMSclasses.GAMSSet('fitFrags',set(fitFragNames))       
        
        # Gamma matrices (to take care of derivatization effects, as per Wahl et al Biotechnol Bioeng 85: 259-268 (2004))
        # m and n sets
        twelveset = ['0','1','2','3','4','5','6','7','8','9','10','11','12']
        mset   =  GAMSclasses.GAMSSet('mset',set(twelveset))
        nset   =  GAMSclasses.GAMSSet('nset',set(twelveset))

        # gamma matrix
        Elements = {}
        for abbrev in fragDict:
            # if data is processed, gamma is unity and changes nothing
            if procString == 'proc':
                gamma = labeling.createGammaMatrix(len(twelveset),'')
            # Else the gamma matrix encodes the derivatization info
            elif procString == 'raw':
                frag  = fragDict[abbrev]
                gamma = labeling.createGammaMatrix(len(twelveset),frag.comp)
            else:
                raise Exception('Unrecognized procString: '+procString)
            # GAMS par elements    
            for nind in range(len(twelveset)):
                for mind in range(len(twelveset)):
                    Elements[tuple([abbrev,str(nind),str(mind)])] = gamma[nind,mind]

        gammapar = GAMSclasses.GAMSParameter('gamma',Elements)            

        # File output
        aacorrfilename = 'aacorr.txt'
        aacorrSt = aacorr.write('toString')
        aacorrFile = aacorrfilename
        
        msetfilename = 'mset.txt'
        msetSt = mset.write('toString')
        msetFile  = msetfilename
        nsetfilename = 'nset.txt'
        nsetSt = nset.write('toString')
        nsetFile  = nsetfilename

        ncarbonsfFilename = 'frag_carbon_numbers.txt'
        ncarbonsfSt = ncarbonsf.write('toString')
        ncarbonsfFile = ncarbonsfFilename
        
        fragmentsSetFilename = 'fragmentsMS.txt'
        fragmentsSetSt = fragmentsSet.write('toString')
        fragmentsSetFile = fragmentsSetFilename
        
        fitFragsSetFilename = 'fitFrags.txt'
        fitFragsSetSt = fitFragsSet.write('toString')
        fitFragsSetFile = fitFragsSetFilename        
        
        gammafilename = 'gamma_mat.txt'
        gammaparSt = gammapar.write('toString')
        gammaFile = gammafilename
        return [(aacorrFile,aacorrSt),(msetFile,msetSt),(nsetFile,nsetSt),(ncarbonsfFile,ncarbonsfSt),
                (fragmentsSetFile,fragmentsSetSt),(fitFragsSetFile,fitFragsSetSt),(gammaFile,gammaparSt)]


    def getSourceLabelFile(self):
        """
        Files for feed labeling
        """
        # all EMUs and carbon dictionary
        rawFeedDict  = self.metList.getFeedDict()
        allEMUs      = self.reactionList.emuList.emus

        # TODO: it would be better to create the carbonDict fresh every time
        #carbonDict   = self.reactionList.getCarbonDict()
        carbonDict   = self.reactionList.carbonDict

        # Store all feeds in a feed dictionary
        feedDict = {}
        for metName in rawFeedDict:
            if rawFeedDict[metName].strip() == 'free': 
                continue
            # Eliminate first part and work with rest, e.g.: string2 = 30% 1-C 20% U 50% UN 
            string2  = rawFeedDict[metName].split(':')[1]          
            # Get feed components, e.g.: comps = ['30% 1-C', '20% U', '50% UN']
            comps    = re.findall('(\d+% [\w\,\-]+)',string2) 
            compDict = {}
            for comp in comps:
                compName = comp.split(' ')[1] # e.g.: UN or 1-C
                compPerc = old_div(float(comp.split(' ')[0].replace('%','')), 100) # e.g.: 0.30
                compDict[compName] = compPerc
            # Feed Dict , e.g. {'Glu': {'1-C': 0.3, 'UN': 0.5, 'U': 0.2}}    
            feedDict[metName] = compDict 

        # Test they all add up to 100
        for metName in feedDict:
            total = 0
            for label in feedDict[metName]:
                total = total + feedDict[metName][label]
            if total != 1.0:
                raise Exception('Sum of feeds for '+metName+' does not add up to 100%: '+str(feedDict[metName]))   
   
        # Determining MDVs for source metabolites
        SourceEMUsMDVDict = {}

        emuCount = 0
        for emu in allEMUs:
            emuCount = emuCount + 1
            metName = emu.getMetName()
            emuName = emu.name
            if metName in feedDict:
                if carbonDict[metName] == 0:
                    raise Exception('Emu number %s "%s" has 0 in carbonDict for metName "%s"' % (str(emuCount), emuName, metName))                       
                MDVtest = labeling.LabelinglString('UN', carbonDict[metName]).mdv(emuName)  # Used to get length only
                MDVmat = numpy.zeros(shape=(len(feedDict[metName]),len(MDVtest)))
                counter = 0
                for comp in feedDict[metName]:
                    s = labeling.LabelinglString(comp, carbonDict[metName])
                    MDVmat[counter] = s.mdv(emuName)*feedDict[metName][comp]
                    counter  = counter+1
                # MDVs dict, e.g.: {'accoa_c_1_2': array([ 0.75,  0.25,  0.  ])}
                SourceEMUsMDVDict[emuName] = numpy.sum(MDVmat, axis=0) 
        
        # Writing file with values
        SourceLabelingStr = ''
        for emu in sorted(SourceEMUsMDVDict):
            counter = 0
            for value in SourceEMUsMDVDict[emu]:
                line = "f.fx('" + emu + "','" + str(counter) + "')=" + str(value) + ";\n"
                SourceLabelingStr = SourceLabelingStr + line
                counter = counter+1
            SourceLabelingStr = SourceLabelingStr + '*\n'
        SourceLabelingStr = SourceLabelingStr + '*\n'
        
        SrcLabelFile = 'Source_Labeling.txt'
        return [(SrcLabelFile, SourceLabelingStr)]


    def getStoichMetRxn13CFiles(self,vSuggested='default'):
        """
        Produces stoichiometry, metabolites and reactions files for 13C MFA problem
        """

        FileOut = []
        # S matrix, mets and reacs
        #stoichMatrix  = self.reactionList.stoichMatrix
        stoichMatrix = self.reactionList.getStoichMatrix('S',2)

        fileS = 'StoichNetwork.txt'
        stS   = stoichMatrix.gamsParameter.write('toString')        
                
        metsSet = GAMSclasses.GAMSSet('mets',set(stoichMatrix.mets))
        rxnsSet = GAMSclasses.GAMSSet('rxns',set(stoichMatrix.rxns))
        filemets = 'mets.txt'
        filerxns = 'rxns.txt'
        stmets = metsSet.write('toString')
        strxns = rxnsSet.write('toString')

        FileOut.extend([(fileS,stS),(filemets,stmets),(filerxns,strxns)])

        # Exchange metabolites
        ExchaMet = stoichMatrix.getExchangeMets()
        
        filename    = 'ExchaMets.txt'
        ExchaMetSet = GAMSclasses.GAMSSet('Exchamets',set(ExchaMet))
        stMetSet = ExchaMetSet.write('toString')

        FileOut.extend([(filename,stMetSet)])
        
        # Measured fluxes
        MeasFluxes = self.reactionList.getMeasuredFluxes(level=2)   
        jmeasSet = GAMSclasses.GAMSSet('jmeas',set(MeasFluxes.keys()))
        
        
        #     UPPER and LOWER bounds
        ElementsHI = {}
        ElementsLO = {}
        for rxn in MeasFluxes:
            HI = MeasFluxes[rxn].hi
            LO = MeasFluxes[rxn].lo 
            ElementsHI[tuple([rxn])] = HI
            ElementsLO[tuple([rxn])] = LO
            # Error checking (this should be its own method)
            if (HI<LO) or ((LO<0 or HI<0) and (not 'EX' in rxn)):  # Make exception for exchange reactions
                raise Exception('Reaction '+rxn+' has unacceptable bounds: '+str(LO)+','+str(HI))
        
        vUPPERpar = GAMSclasses.GAMSParameter('vUPPER',ElementsHI)
        vLOWERpar = GAMSclasses.GAMSParameter('vLOWER',ElementsLO)
      
        jmeasFile  = 'jmeas.txt'
        jmeasSt    = jmeasSet.write('toString')
        vUPPERFile = 'vUPPER.txt'
        vUPPERSt   = vUPPERpar.write('toString')        
        vLOWERFile = 'vLOWER.txt'
        vLOWERSt   = vLOWERpar.write('toString')
        FileOut.extend([(jmeasFile,jmeasSt),(vUPPERFile,vUPPERSt),(vLOWERFile,vLOWERSt)])

        # Excluded metabolites
        excludedMetabs   = self.metList.getExcluded()
        ExcMetsSet       = GAMSclasses.GAMSSet('ExcMets',set(excludedMetabs))
        filenameout      = ExcMetsSet.name+'.txt'
        ExcMetsSetSt = ExcMetsSet.write('toString')
        ExcMetFile = filenameout
        FileOut.extend([(ExcMetFile,ExcMetsSetSt)])
        
        # Source Metabolites  
        sourceMetabs   = self.metList.getSource()
        SourceMetsSet  = GAMSclasses.GAMSSet('SourceMets',set(sourceMetabs))
        filenameout    = SourceMetsSet.name+'.txt'
        SourceMetsSetSt = SourceMetsSet.write('toString')
        SourceMetsFile = filenameout
        FileOut.extend([(SourceMetsFile,SourceMetsSetSt)])
        
        # Suggested fluxes
        filenameSet = 'jsugg.txt'
        filenamePar = 'vSUGG.txt'
        if vSuggested == 'default':
            jSUGGset  = GAMSclasses.GAMSSet('jsugg',set([]))    
            vSUGGpar  = GAMSclasses.GAMSParameter('vSUGG',[])
        else:
            assert vSuggested.__class__.__name__=='fluxGAMSpar', 'vSuggested must be a fluxGAMSpar class'   
            reactionList = vSuggested.getReactionList()
            jSUGGset = GAMSclasses.GAMSSet('jsugg',set(reactionList.getReactionNameList(level=2)))
            vSUGGpar = vSuggested

        jSUGGsetSt = jSUGGset.write('toString')
        vSUGGparSt = vSUGGpar.write('toString')
        
        FileOutSet = filenameSet
        FileOutPar = filenamePar
        FileOut.extend([(FileOutSet,jSUGGsetSt), (FileOutPar,vSUGGparSt)])
        
        return FileOut    


    # Output emu data
    def getEMUfiles(self):
        """
        EMU (elementary metabolite unit) files output
        """
        transitions = self.reactionList.transitions

        filesList = []

        equivEMUs   = transitions.equivEMUs
        dummyEMUs   = transitions.dummyEMUs
        emuList     = transitions.getEMUList(equivEMUs.elements)

        # EMU basics
        allEMU     = emuList.getAllEMUSet()
        metEMU     = emuList.getMetEMUPar()
        emuCarbons = emuList.getEMUNCarbons()

        filesList.append( ('all_emu.txt', allEMU.write('toString')) )
        filesList.append( ('met_emu.txt', metEMU.write('toString')) )
        filesList.append( ('emu_carbon_numbers.txt', emuCarbons.write('toString')) )

        # dummy and equivalent emus
        filesList.append( ('equivalent_emus.txt', equivEMUs.write('toString')) )
        filesList.append( ('dummy_emus.txt', dummyEMUs.write('toString')) )
        
        # Condensation equations
        eqns,eqnNames = transitions.getCondGAMSReacs()

        lineEqns     = []
        lineEqnNames = []  
        
        for line in eqns:
            lineEqns.append(line)
        for line in eqnNames:
            lineEqnNames.append(line)

        filesList.append( ('condensation_rxns_eqns.txt', ''.join(lineEqns)) )
        filesList.append( ('condensation_rxns_names.txt', ''.join(lineEqnNames)) )
        
        # EMM matrix
        emmPar  = transitions.getEMMmat()
        filesList.append( ('emm_simplified.txt', emmPar.write('toString')) )

        return filesList


    def getAuxiliaryFiles(self,labelCompNormMin=0.99):
        """
        Auxiliary files such as the random seed
        """
        
        # Random seed
        filename1 = 'randseed.txt'      
        line1 = 'randseed = '+str(random.randint(1,1000000000))+';'
                
        # Minimum for labelCompNorm TODO(Hector): eliminate this one. Not used anymore
        filename2 = 'labelCompNormMin.txt'
        line2 = 'labelCompNormMin = '+str(labelCompNormMin)+';'  
        
        return [(filename1,line1),(filename2,line2)]


    # Reactions map from genome-scale model to 13C model
    def getReacMapFiles(self):               
        filename1 =  'reacmap.txt'            
        filename2 =  'linkedreacs.txt'

        reacMapPar,linkedReacsSet = self.reactionList.getReactionMap()
        reacMapParSt     = reacMapPar.write('toString')
        linkedReacsSetSt = linkedReacsSet.write('toString')

        return [(filename1,reacMapParSt),(filename2,linkedReacsSetSt)]


    # Output all 13CMFA files
    def get13CMFAfiles(self, labelCompNormMin=0.99, vSuggested='default', procString='raw'):
        "Output of all needed 13C MFA files "

        if not vSuggested == 'default':
            vSuggested = vSuggested[0]
        
        filesOut = []
        filesOut.extend(self.getStoichMetRxn13CFiles(vSuggested=vSuggested))
        filesOut.extend(self.getEMUfiles())
        filesOut.extend(self.getSourceLabelFile())
        filesOut.extend(self.getAuxiliaryFiles(labelCompNormMin=labelCompNormMin))
        filesOut.extend(self.getLabelingDataFiles())
        filesOut.extend(self.getFragmentInfoFiles(procString))
        filesOut.extend(self.getReacMapFiles())    
        filesOut.extend(self.getSolverOptFiles())

        return filesOut

    def secureFluxToDest(self):
        """
        This method secures a minimum flux to destination metabolites to make sure labeling is not left unconstrained.
        I do this here instead of reactionlist since this only really pertains to 13C networks.
        """        
        
        metList      = self.metList
        reactionList = self.reactionList
        
        # Get destination metabolites
        destMets = metList.getDestinationMets()
        
        # Prop up lower bound to at least 10^-7 for reaction flowing into each destination metabolite
        for reaction in reactionList.reactions:
            newBound = 0.000001
            prods  = reaction.products
            reacts = reaction.reactants
            
            if not reaction.reversible:
                for prod in prods:
                    if prod.name in destMets:
                        lb = reaction.fluxBounds.net.lo
                        reaction.fluxBounds.net.lo = newBound if (lb <= 0) else lb    
                        reaction.fluxBounds.net.hi = max(newBound,reaction.fluxBounds.net.hi)
                        reaction.measured = True
            else:
                for prod in prods:
                    if prod.name in destMets:
                        lb = reaction.fluxBounds.forward.lo
                        reaction.fluxBounds.forward.lo = newBound if (lb <= 0) else lb    
                        reaction.fluxBounds.forward.hi = max(newBound,reaction.fluxBounds.forward.hi)                       
                        reaction.measured = True
                for react in reacts:
                    if react.name in destMets:
                        lb = reaction.fluxBounds.backward.lo
                        reaction.fluxBounds.backward.lo = newBound if lb <= 0 else lb
                        reaction.fluxBounds.backward.hi = max(newBound,reaction.fluxBounds.backward.hi)
                        reaction.measured = True


    def changeStds(self,inputData):
        """
        Changes standard deviations for labeling data to notes given in input or to the fit, if available
        See Garcia Martin et al PLOS Comp Biol. 2015 Equation 23
        """
        # TODO: This would fit better as a Notes method.
        
        # If meant to adapt to the fits
        if (inputData.__class__.__name__=='str' and inputData=='toFits'):
            fragDict = self.fragDict() 
        
            dataTypes = ['GCMSLabelData','CEMSLabelData','LCMSLabelData']
                
            for dataType in dataTypes:
                if dataType in self.notes:
                    data    = self.notes[dataType]
                    for abbrev in data:
                        newFrag    = data[abbrev]
                        frag       = fragDict[abbrev]
                        # Increasing the error if std is smaller than the difference to the fit                        
                        for i in range(len(newFrag.std)):                        
                            delta = abs(frag.mdv[i]-frag.MDVfit[i])
                            if delta > newFrag.std[i]:
                                newFrag.std[i] = 1.1*delta
                            else:
                                newFrag.std[i] = newFrag.std[i]
        
        # If meant to adapt to an input fragment dictionary
        elif inputData.__class__.__name__=='dict':   # Assuming fragdict dictionary
            fragDictOut = inputData
            dataTypes = ['GCMSLabelData','CEMSLabelData','LCMSLabelData']
                
            for dataType in dataTypes:
                if dataType in self.notes:
                    data    = self.notes[dataType]
                    for abbrev in data:
                        newFrag     = data[abbrev]
                        fragOut     = fragDictOut[abbrev]
                        newFrag.std = fragOut.std
            
        else:
            raise Exception('Unrecognized input for labeling standard deviations modification')
                            

    def excludeFits(self,mets):
        """Exclude the provided metabolites (mets) from fits. The input is a name or list of names"""
        
        if mets.__class__.__name__=='str':  # Single mets to be eliminated
            metsOut = [mets]
        elif mets.__class__.__name__=='list':
            metsOut = mets   # list of mets to be eliminated

        dataTypes = ['GCMSLabelData','CEMSLabelData','LCMSLabelData']
        for name in metsOut:
            present = False
            for dtype in dataTypes:
                if dtype in self.notes:
                    if name in self.notes[dtype]:
                        self.C13ReacNet.notes[dtype][name].inFit = False   # This eliminates met name from fit 
                        self.notes[dtype][name].inFit = False    
                        present = True
                    if not present:
                        raise Exception('Metabolite '+name+' not present among labeling measurements: \n'+str(self.notes[dtype]))



class TSReactionNetwork(C13ReactionNetwork):
    """ A derived class for Two-Scale 13C sbml files """
    
    def __init__(self, input_, inFluxRef=['EX_glc(e)','EX_glc_e_']):
        # sbml files and input tuples =(name, Notes, metList, ReactionList) contain the same information
        debug = False        

        # Proceed as for general SBML file
        ReactionNetwork.__init__(self, input_)
        if debug:
            self.write('BeginningC13linesTS.sbml')

        # Obtain reference flux uptake rate 
        self.inFluxRef = self.getFluxRefScale(inFluxRefs=inFluxRef)

        # Process carbon transitions for 13C system
        C13ReactionList = self.reactionList.getCarbonTransitionReactions(self.inFluxRef)    
        
        C13ReacNet = C13ReactionNetwork( (self.name, self.notes, C13ReactionList.getMetList(), C13ReactionList))  
        if debug:        
            C13ReacNet.reactionList.writeReactionFile('C13linesTS.txt')
            C13ReacNet.write('C13.sbml')
        self.C13ReacNet = C13ReacNet
        
        # Check transitions reaction line completeness
        metsCore = enhancedLists.MetaboliteList(set(self.C13ReacNet.metList.mets) - set(self.C13ReacNet.metList.getExcluded()))
        metsCoreDict = metsCore.metDict
        self.reactionList.checkTransitionLines(metsCoreDict)
        
        # Secure flux to destination metabolites
        self.C13ReacNet.secureFluxToDest()


    def addTransitions(self,filename,prioritySBML=False): 
        """ 
        Adds reactions into SBML file from reaction file of type:
            
        &MLIST s7p 0
        &MLIST mal-L 0

        &SOURCE glc-D[e]

        # Input Reactions
        EX_glc(e)	glc-D[e] <==> glcDEx	abcdef : abcdef
        GLCt2	glc-D[e] --> glc-D	abcdef : abcdef
        HEX1	glc-D --> g6p	abcdef : abcdef
      
      
        prioritySBML decides which file has priority in reaction reversibility. 
         """

        # Opening file
        Ctransitions      = enhancedLists.AtomTransitionList(filename)
        EMUDestMetabList  = Ctransitions.EMUDestMetabList
        EMUExcludedMetabs = Ctransitions.EMUExcludedMetabs

        # Dictionary for Destination Metabolites    
        dest_metabolite_dict = {}
        for met in EMUDestMetabList:
            if met[0] in dest_metabolite_dict:
                dest_metabolite_dict[met[0]].append(met[1:])
            else:
                dest_metabolite_dict[met[0]] = [met[1:]]

        # Get metabolite dict
        metList  = self.metList
        metsDict = metList.metDict

        # Adding destination and excluded tags
        for metName in dest_metabolite_dict:
            met = core.Metabolite(core.MetaboliteName(metName))
            metsDict[met.name].destination = True
            metsDict[met.name].carbonsOut = dest_metabolite_dict[metName]
        for metName in EMUExcludedMetabs:
            met = core.Metabolite(core.MetaboliteName(metName))
            metsDict[met.name].destination = True
            metsDict[met.name].excluded = True

        # Get reaction dict
        reactionsDict = self.reactionList.getReactionDictionary()

        # Adding transitions to reactions
        for trans in Ctransitions.transitions:
            transLine = core.AtomTransition(trans.getTransIn('SBML'))            
            
            try:
                reaction = reactionsDict[transLine.name]
            except KeyError:
                    raise Exception('Reaction '+transLine.name+' from carbon transitions file not present in SBML file')
                    
            reaction.transitionLine = transLine
            # Sincronyze info for both reaction and reaction line
            reaction.syncTransitionLine(prioritySBML)


    def getFluxRefScale(self,inFluxRefs=[],bounds=True):
        """
        Obtains reference flux information (e.g. glucose input)
        """

        reactionDict = self.reactionList.getReactionDictionary()
        
        # TODO(Hector): test for ethanol or acetate inputs
        # Keep trying all inFluxRefs inputs until one is found in the reaction dictionary. 
        if inFluxRefs:
            maxTries = 100     
            trying = True
            tries  = 0 
            while trying and (tries < maxTries):
                for inp in inFluxRefs:
                    tries        +=1
                    try:
                        inReaction   = reactionDict[inp]
                        trying       = False
                    except:
                        pass
        
            if tries >= maxTries:
                raise Exception('Could not find an input reaction in: '+str(inFluxRefs))      
        else:
            inputFluxName = self.getHighestCarbonInput(bounds=bounds)
            if not inputFluxName:
                raise Exception('Could not find reference flux:'+str(inputFluxName))
            inReaction = reactionDict[inputFluxName]
        
        # Decide wether to use bounds or net flux to give the reference
        if bounds:
            inFluxRef = abs(old_div((inReaction.fluxBounds.net.hi+inReaction.fluxBounds.net.lo),2))   
        else:
            inFluxRef  = abs(inReaction.flux.getComp('net','best'))

        return inFluxRef


    def getReacsFlowingIntoCore(self):
        """
        Provides the list of fluxes flowing into core metabolism
        """
        # Core metabolites and all reactions for genome scale model
        coreMetList      = self.C13ReacNet.metList
        coreMetList.Dict = coreMetList.metDict 
        allReactions = self.reactionList
     
        # Screen out metabolites not involved in carbon transitions (e.g. nadh)
        carbonDict = allReactions.getCarbonDict()
        newList = []
        for met in coreMetList.mets:
            if met.name in carbonDict:
                newList.append(met)
        coreMetList = enhancedLists.MetaboliteList(newList)
        coreMetList.Dict = coreMetList.metDict 

        # Main loop
        reacsFlowingIntoCore = []
        for reaction in allReactions.reactions:
            reaction.hasForwardFlow  = False
            reaction.hasBackwardFlow = False
            
            for prod in reaction.products:
                if prod.name in coreMetList.Dict:
                    reaction.hasForwardFlow = True
                    reacsFlowingIntoCore.append(reaction)
            if reaction.reversible:
                for reac in reaction.reactants:
                    if reac.name in coreMetList.Dict:
                        reaction.hasBackwardFlow = True
                        reacsFlowingIntoCore.append(reaction)
        reacsFlowingIntoCore = set(reacsFlowingIntoCore)

        # Eliminate reactions in core
        newReacs = []
        coreDict = self.C13ReacNet.reactionList.getReactionDictionary()
        
        for reaction in reacsFlowingIntoCore:
            if not (reaction.name in coreDict):
                newReacs.append(reaction)
           
        reacsFlowingIntoCore = set(newReacs) 
        
        # Make into a reactionList instance      
        return enhancedLists.ReactionList(sorted(list(reacsFlowingIntoCore)))


    # Output all 2S-13CMFA files
    def getTwoS13CMFAfiles(self, labelCompNormMin=0.99, vSuggested='default', procString='raw'):
        """Output of TS files."""

        if not vSuggested == 'default':
            vSuggestedA   = vSuggested[0]
            vFBASuggested = vSuggested[1]
        else:
            vSuggestedA   = 'default'
            vFBASuggested = 'default'

        # 13C
        filesOut = []
        filesOut.extend(self.C13ReacNet.getStoichMetRxn13CFiles(vSuggested=vSuggestedA))
        filesOut.extend(self.C13ReacNet.getEMUfiles())
        filesOut.extend(self.C13ReacNet.getSourceLabelFile())
        filesOut.extend(self.C13ReacNet.getAuxiliaryFiles(labelCompNormMin=labelCompNormMin))
        filesOut.extend(self.C13ReacNet.getLabelingDataFiles())
        filesOut.extend(self.C13ReacNet.getFragmentInfoFiles(procString))
        filesOut.extend(self.getReacMapFiles())       
        filesOut.extend(self.getSolverOptFiles())
        
        # 2S
        filesOut.extend(self.getStoichMetRxnFBAFiles(vSuggested=vFBASuggested))
        filesOut.extend(self.getInFluxRefFile())

        return filesOut        
    

    def getReacMapFiles(self):
        """
        Provides the files containing the map that connects reactions in the genome-scale model to the 13C model
        """
        filename1 = 'reacmap.txt'            
        filename2 = 'linkedreacs.txt'

        reacMapPar,linkedReacsSet = self.C13ReacNet.reactionList.getReactionMap()
        reacMapParSt     = reacMapPar.write('toString')
        linkedReacsSetSt = linkedReacsSet.write('toString')

        return [(filename1,reacMapParSt),(filename2,linkedReacsSetSt)]


    def getInFluxRefFile(self):
        """Finds influx reference flux information and puts out the corresponding file."""
        # Reference flux (e.g glucose) uptake calculation
        inFluxRef = self.getFluxRefScale(bounds=True)       
        
        if inFluxRef == 0:
            raise Exception('Flux reference scale cannot be zero') 
        
        # Output to file
        filename = 'inFluxRef.txt'
        line = 'inFluxRef = '+str(inFluxRef)+';'
        
        return [(filename,line)]


##### UNIT TESTS ########

import unittest

# Expected result from tests below
expectedResult =  '153.68'  



class testRN(unittest.TestCase):  
    """ Testing basic usage of Reaction Network """

    def setUp(self):
        
         GLCpts = core.Reaction.from_string('GLCpts: glc_dsh_D_e + pep_c --> g6p_c + pyr_c')
         PGI    = core.Reaction.from_string('PGI: g6p_c <==> f6p_c')
         PFK    = core.Reaction.from_string('PFK: atp_c + f6p_c --> adp_c + fdp_c + h_c')  
           
         reactionList = enhancedLists.ReactionList([GLCpts,PGI,PFK])
         self.ReactionNetwork = ReactionNetwork(('test',{},reactionList.getMetList(),reactionList))           


    def testRxnsMets(self):
        """ Testing expected reactions and metabolites in basic usage"""
        self.assertTrue(str(self.ReactionNetwork.reactionList) == 'GLCpts\nPFK\nPGI')
        self.assertTrue(str(self.ReactionNetwork.metList)      == 'adp_c\natp_c\nf6p_c\nf6p_c\nfdp_c\ng6p_c\ng6p_c\nglc_dsh_D_e\nh_c\npep_c\npyr_c')


    def testAddReaction(self):
        """ Testing result of adding reactions"""
        self.ReactionNetwork.addReaction('F6PA: f6p_c <==> dha_c + g3p_c')
        
        self.assertTrue(str(self.ReactionNetwork.reactionList) == 'F6PA\nGLCpts\nPFK\nPGI')
        self.assertTrue(str(self.ReactionNetwork.metList)      == 'adp_c\natp_c\ndha_c\nf6p_c\nf6p_c\nfdp_c\ng3p_c\ng6p_c\ng6p_c\nglc_dsh_D_e\nh_c\npep_c\npyr_c')    


    def testChangeFluxBounds(self):
        """ Testing changeFluxBounds in basic usage"""
        self.ReactionNetwork.changeFluxBounds('PFK',(0.1,0.5))
        rxnDict = self.ReactionNetwork.reactionList.getReactionDictionary()
        
        self.assertTrue(str(rxnDict['PFK'].fluxBounds.net) == '[0.1:0.3:0.5]')


    def tearDown(self):
        pass



class testC13RN(unittest.TestCase):  
    """ Testing usage of C13 Reaction Network """

    def setUp(self):
        
         self.dirBase = os.environ['QUANTMODELPATH']+'/data/tests/TCAtoy/' 
         REACTIONSfilename   = self.dirBase+'REACTIONStca.txt' 
           
         atomTransitions = enhancedLists.AtomTransitionList(REACTIONSfilename)
         self.reactionNetwork = atomTransitions.getReactionNetwork('E. coli wt5h 13C MFA')


    def testAddFeedLabel(self):
        """ Testing Adding feed and labeling"""
        
        FEEDfilename        = self.dirBase+'FEEDtca.txt'
        CEMSfilename        = self.dirBase+'GCMStca.txt'
        CEMSSTDfilename     = self.dirBase+'GCMSerrtca.txt'
        FLUXESFreefilename  = self.dirBase+'FLUXtca.txt'        
        
        self.reactionNetwork.addLabeling(CEMSfilename,'LCMSLabelData',CEMSSTDfilename,minSTD=0.001)
        self.reactionNetwork.addFeed(FEEDfilename)
        self.reactionNetwork.loadFluxBounds(FLUXESFreefilename)

        self.assertTrue(str(self.reactionNetwork.notes['LCMSLabelData']['Glu'].mdv)  == '[ 0.346  0.269  0.27   0.081  0.028  0.004]')
        self.assertTrue(str(self.reactionNetwork.notes['LCMSLabelData']['Glu'].std)  == '[ 0.02  0.02  0.02  0.02  0.02  0.02]')

                  
    def tearDown(self):
        pass    



class testTSRN(unittest.TestCase):  
    """ Testing usage of 2S C13 Reaction Network """

    def setUp(self):
        
        self.datadir = os.environ["QUANTMODELPATH"]+'/data/tests/Toya2010/2S/wt5h/'

        BASEfilename      = self.datadir + 'EciJR904TKs.xml'
        FLUXESfilename    = self.datadir + 'FLUXwt5h.txt'

        # Load initial SBML file
        reactionNetwork = TSReactionNetwork(BASEfilename)
    
        # Add Measured fluxes
        reactionNetwork.loadFluxBounds(FLUXESfilename)

        self.RN = reactionNetwork


    def testAddtransitions(self):
        """ Testing adding transitions to genome-scale"""

        REACTIONSfilename = self.datadir + 'REACTIONSwt5h.txt' 
        MSfilename        = self.datadir + 'GCMSwt5h.txt'
        FEEDfilename      = self.datadir + 'FEEDwt5h.txt'
        MSSTDfilename     = self.datadir + 'GCMSerrwt5h.txt'

        # Add carbon transitions
        self.RN.addTransitions(REACTIONSfilename)
        # Add measured labeling information
        self.RN.addLabeling(MSfilename,'LCMSLabelData',MSSTDfilename,minSTD=0.001)
        # Add feed labeling information
        self.RN.addFeed(FEEDfilename)

        rxnDict = self.RN.reactionList.getReactionDictionary()

        self.assertTrue(str(rxnDict['PDH'].stoichLine())  == 'PDH : coa_c + nad_c + pyr_c --> accoa_c + co2_c + nadh_c')
        self.assertTrue(str(rxnDict['PDH'].transitionLine)  == 'PDH\tpyr_c --> co2_c + accoa_c\tabc : a + bc')


    def tearDown(self):
        pass    

