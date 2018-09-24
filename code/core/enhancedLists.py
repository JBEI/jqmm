# Part of the JBEI Quantitative Metabolic Modeling Library (JQMM)
# Copyright (c) 2016, The Regents of the University of California.
# For licensing details see "license.txt" and "legal.txt".
"""
The enhancedLists module provides classes with lists that have been enhanced with methods useful for the objects
therein. So, e.g. the reactionList class has a carbonTransitionsOK method that tests if carbon transition and reaction
information are compatible.
"""
from __future__ import print_function
from __future__ import division


from builtins import zip
from builtins import str
from builtins import range
from past.utils import old_div
from builtins import object
import core, GAMSclasses, ReactionNetworks
import math, copy, re
import numpy
import utilities as utils


class MetaboliteList(object):
    "Class for list of metabolites in a reaction network. This is just a metabolite list with some useful functions atached to it."

    def __init__(self,mets):
        """
        Input is just a list of metabolites
        """
        # Sort mets by name for standardization
        self.mets = sorted(mets,key=lambda met: met.name)
        self._remakeDict()


    def __str__(self):
        lines = []
        for met in self.mets:
            lines.append(str(met))
        return "\n".join(lines)


    def __iter__(self):
        for met in self.mets:
            yield met


    def _remakeDict(self):
        self.metDict  = {}
        for met in self.mets:
            self.metDict[met.name] = met


    def addMetabolites(self, metabolites):
       """
       Adds metabolite to list.
       Input should be a list of Metabolite objects
       """
       for m in metabolites:
            if m.name not in self.metDict:
                self.mets.append(m)
       newMetabolites = self.mets
       self.mets = sorted(newMetabolites,key=lambda met: met.name)
       self._remakeDict()


    def getCarbonDict(self):
        """
        Returns carbon dictionary with number of carbons for each metabolite
        """
        carbonDict = {}        
        for met in self.mets:
            carbonDict[met.name] = met.ncarbons
        return carbonDict    


    def getCompartments(self):
        """
        Returns an array of all compartment identifiers used in this Metabolite list
        """
        compDict = {}
        for met in self.mets:
            if hasattr(met, 'compartment'):
                compDict[met.compartment] = True
            # Some Metabolites may be from old shelved instantiations and not have a component set
            elif isinstance(met.name, core.MetaboliteName):
                compDict[met.name.compartment] = True
            else:
                temp_name = core.MetaboliteName(met.name)
                compDict[temp_name.compartment] = True
        return list(compDict.keys())


    def getExcluded(self):
        """
        Returns list of excluded metabolites (metabolites not included in EMU balance)
        """
        mets     = self.mets
        excluded = []        
        for met in mets:
            if met.excluded:
                excluded.append(met.name) 
        return excluded


    def getSource(self):
        """
        Returns list of source metabolites
        """
        mets   = self.mets
        source = []        
        for met in mets:
            if met.source:
                source.append(met.name) 
        return sorted(list(set(source)))  # return unique list of metabolites (sorted)


    def getFeedDict(self):
        """
        Returns dictionary of feeds (feed = labeling for source metabolite)
        """
        mets     = self.mets
        feedDict = {}        
        for met in mets:
            if met.source:
                feedDict[met.name] = met.feed 
        self.feedDict = feedDict
        return feedDict


    def getDestinationMets(self):
        """
        Returns list of destination metabolites (measured metabolites)
        """
        mets        = self.mets
        destination = []        
        for met in mets:
            if met.destination:
                destination.append(met.name) 
        return destination
 

    def getDestinationEMUs(self):
        """
        Returns list of destination EMUs (measured EMUs)
        """
        destinationEMUs = []
        mets            = self.mets
              
        for met in mets:
            if met.destination:
                carbonsOut = met.carbonsOut
                for carb in carbonsOut:
                    destinationEMUs.append(self.metDict[met.name].generateEMU(carb))  
  
        destinationEMUs = sorted(list(set(destinationEMUs)))  # Eliminating repeats
        return destinationEMUs


    def getBvec(self): # TODO: rethink if this can be done more elegantly
        """
        Provides objective vector b for LP problem of the type sum_j S_{ij}*v_j = b_i
        """
        ElementsBvec = {}
        for met in self.mets:
            b = met.bval
            if met.name[-2:] !='_b':     # This should probably be done through exchange metabolites
                ElementsBvec[tuple([met.name])] = b
        bVecPar = GAMSclasses.GAMSParameter('bvec',ElementsBvec)
        return bVecPar


    def updateNcarbons(self,carbonDict):
        """
        Updates number of carbons in MetaboliteList with carbonDict information, IF not already present
        """
        for name in self.metDict:
            if name in carbonDict:
                if self.metDict[name].ncarbons <= 0:
                    self.metDict[name].ncarbons = carbonDict[name]



class ReactionList(object):
    """
    Class for a list of Reactions in a reaction network.
    This is just a reaction list with some useful functions attached to it.
    """

    def __init__(self, reactions, sort=True):
        """
        input is a list of reactions, sort indicates if they need to be sorted
        """
        # Sort reactions by name for standardization
        self.reactions = reactions
        if sort:
            # Sort reactions by name for standardization
            self.sort('by name')
            
        self.revReacs   = self.getRevReacs()
        self.SuggFluxes = self.getSuggFluxes()
        self.MeasFluxes = self.getMeasuredFluxes()
        self.reacNames1 = self.getReactionNameList(level=1)
        self.reacNames2 = self.getReactionNameList(level=2)


    def __str__(self):   # TODO(Hector): may want to reconsider this (low priority)
        lines = []
        for reaction in self.reactions:
            lines.append(str(reaction))
        return "\n".join(lines)


    def __iter__(self):
        for reaction in self.reactions:
            yield reaction


    def __mul__(self,other):
        "ReactionList multiplication"
        for reaction in self:
            reaction = other * reaction
        return self


    def __rmul__(self,other):
        "reaction reverse multiplication"
        result = self.__mul__(other)
        return result    


    def __div__(self,other):
        "reaction reverse divison"
        result = self.__rmul__(old_div(1,other))
        return result    
    

    def getStoichMatrix(self,parName='Sbig',level=1):    
        "Creates stoichiometry matrix from reaction information, level is the level of description: 1 -> PGI ; 2 -> PGI_f, PGI_b"
        reactions = self.reactions
                
        Elements  = {}
        for reaction in reactions:
            if level == 2 and reaction.reversible and not reaction.exchange:  # Include forward and backward reactions
                forwMark = '_f'
                backMark = '_b'
                for react in reaction.reactants:
                    Elements[(react.name,reaction.name+forwMark)] = -react.stoichiometry 
                    Elements[(react.name,reaction.name+backMark)] =  react.stoichiometry                     
                for prod in reaction.products:
                    if not reaction.exchange:                           # eliminate from big stoichiometry matrix if exchange. Could be done in a different way by changing gams code
                        Elements[(prod.name,reaction.name+forwMark)] =  prod.stoichiometry    
                        Elements[(prod.name,reaction.name+backMark)] = -prod.stoichiometry    
            else:
                for react in reaction.reactants:
                    Elements[(react.name,reaction.name)] = -react.stoichiometry 
                for prod in reaction.products:
                    if not reaction.exchange:                           # eliminate from big stoichiometry matrix if exchange. Could be done in a different way by changing gams code
                        Elements[(prod.name,reaction.name)] = prod.stoichiometry
                        
        # Stoich.matrix
        SmatPar = GAMSclasses.GAMSParameter(parName,Elements)
        Smat = core.StoichMatrix(gamsParameter=SmatPar)
        
        #self.Smat = Smat    
        return Smat


    def eliminateDupMets(self):   
        "Eliminates duplicate metabolites"
        # Create metabolite dictionary to choose one metabolite instance per metabolite name
        reactions = self.reactions
        metDict = {}
        for reaction in reactions:
            prods  = reaction.products
            reacts = reaction.reactants
            for prod in prods:
                metDict[prod.name] = prod.originalMet
            for react in reacts:
                metDict[react.name] = react.originalMet
        
        # Substitute all original mets for the chosen instance
        for reaction in reactions:
            prods  = reaction.products
            reacts = reaction.reactants
            for prod in prods:
                prod.originalMet = metDict[prod.name]
            for react in reacts:
                react.originalMet = metDict[react.name]
        
        self.reactions = reactions
        self.revReacs  = self.getRevReacs()
        self.SuggFluxes= self.getSuggFluxes()
        self.MeasFluxes= self.getMeasuredFluxes()
        self.reacNames1= self.getReactionNameList(level=1)
        self.reacNames2= self.getReactionNameList(level=2)
        
        return reactions


    def getRevReacs(self):
        "Returns list of reversible reactions"
        reactions  = self.reactions
        RevReacs   = []
        for reaction in reactions:
            if reaction.reversible:
                RevReacs.append(reaction.name)
        return RevReacs


    def getSuggFluxes(self):
        "Returns dictionary of suggested fluxes"
        reactions  = self.reactions
        suggFluxes = {}
        for reaction in reactions:
            if reaction.suggested:
                try:
                    suggFluxes[reaction.name]= reaction.flux
                except AttributeError:
                    raise Exception('Suggested flux not initialized')
        return suggFluxes        


    def getMeasuredFluxes(self,level=1):
        "Returns dictionary of measured fluxes, Level indicates the level of description: 1 -> net flux, 2 -> forward and backward."        
        reactions  = self.reactions
        measFluxes = {}
        for reaction in reactions:             
            if reaction.measured:
                try:                    
                    if level == 1:
                        measFluxes[reaction.name]= reaction.fluxBounds.net
                    elif level == 2:
                        if reaction.reversible and not reaction.exchange:
                            measFluxes[reaction.name+"_f"]= reaction.fluxBounds.forward
                            measFluxes[reaction.name+"_b"]= reaction.fluxBounds.backward
                        else:
                            measFluxes[reaction.name]= reaction.fluxBounds.net
                    else:
                        raise Exception('Wrong level value')
                except AttributeError:
                    raise Exception('Measured flux not initialized')     
        
        return measFluxes


    def getReactionNameList(self,level=1,sort=True):                
        "Returns list of reaction names, e.g.: PDH (level 1) or PDH_f and PDH_b (level 2)"
        reactions    = self.reactions
        reactionList = []
        for reaction in reactions:
            if level == 1:    
                reactionList.append(reaction.name)
            elif level == 2:
                if reaction.reversible:
                    reactionList.append(reaction.name + "_f")
                    reactionList.append(reaction.name + "_b")
                else:
                    reactionList.append(reaction.name)
            else:
                raise Exception('Wrong level')
       
        if sort:
            reactionList = sorted(reactionList)   # This sort should probably be eliminated   
        return reactionList


    def getReactionMap(self):
        "Returns reaction map between level 1 and level 2 reaction names"
        reactions  = self.reactions
        reacMap = {}
        for reaction in reactions:
            if reaction.reversible and not reaction.exchange:
                reacMap[reaction.name] = [reaction.name+'_f',reaction.name+'_b'] 
            else:
                reacMap[reaction.name] = [reaction.name,''] 
                
        # Put into parameters
        Elements = {}
        for name in reacMap:
            positive = reacMap[name][0]
            negative = reacMap[name][1]
            if positive != '':
                key = (name,positive)
                val = 1
                Elements[key] = val
            if negative != '':
                key = (name,negative)
                val = -1
                Elements[key] = val
        
        reacMapPar = GAMSclasses.GAMSParameter('reacmap',Elements)
        linkedReacsSet = GAMSclasses.GAMSSet('linkedreacs',set(reacMap.keys()))
       
        return reacMapPar,linkedReacsSet


    def zeroFluxExchange(self,names='all'):
        "zeros exchange flux for reaction name sgiven"
        for reaction in self:
            if (reaction.name in names) or names=='all':                
                net = reaction.flux.net
                if isinstance(net, core.rangedNumber):
                    exch = core.rangedNumber(0,0,0)
                elif utils.is_float(net):
                    exch = 0
                else:
                    "Unexpected type for net flux: "+str(net)                
                reaction.flux = core.flux(net_exc_tup=(net,exch))   
                

    def fixFluxes(self,level=1,names='all',epsilon=0.01,maximum= math.pow(10,9)):
        "Fix upper and lower fluxes to the values in stored flux"
        reactions = self.reactions
        # maximum should be the same as in GAMS_EMU_MFA_FLUX_HYB.lst, lbvec(jFBA) = max(lbvec(jFBA),-1000);
        # This can be solved by putting bounds on the MOMA results, as well
        
        for reaction in reactions:  
            if (reaction.name in names) or names == 'all':
                if hasattr(reaction,'flux'):
                    if isinstance(reaction.flux.net, core.rangedNumber):
                        netVal = float(reaction.flux.net.best)
                    else:
                        netVal = float(reaction.flux.net)
                    Delta = min(epsilon,abs(netVal)*epsilon)
                    Delta = max(0.00001,Delta)   # This is required so conopt does not collapse lower and upper bounds into one and drives itself nuts
                    
                    hi = netVal+Delta
                    lo = netVal-Delta

                    # maximum here to avoid problems when you fix fluxes
                    if abs(hi)<maximum and abs(netVal)<maximum and abs(lo)<maximum:    
                        reaction.fluxBounds.net.hi   = hi
                        reaction.fluxBounds.net.best = netVal                    
                        reaction.fluxBounds.net.lo   = lo
                    
                    # marking as measured
                    reaction.measured = True                    
                    
                    if level == 1:                    
                        pass
                    
                    elif level == 2:
                        if reaction.reversible:
                            if isinstance(reaction.flux.forward, core.rangedNumber):
                                forwardVal  = float(reaction.flux.forward.best)
                            else:
                                forwardVal = float(reaction.flux.forward)
                            Delta = min(epsilon,abs(forwardVal)*epsilon)

                            hi =     forwardVal+Delta
                            lo = max(forwardVal-Delta,0)

                            if abs(hi)<maximum and abs(netVal)<maximum and abs(lo)<maximum:
                                reaction.fluxBounds.forward = core.rangedNumber(lo, forwardVal, hi)

                            if isinstance(reaction.flux.backward, core.rangedNumber):
                                backwardVal = float(reaction.flux.backward.best)
                            else:
                                backwardVal = float(reaction.flux.backward)
                            Delta = min(epsilon,abs(backwardVal)*epsilon)

                            hi =     backwardVal+Delta
                            lo = max(backwardVal-Delta,0)

                            if abs(hi) < maximum and abs(netVal) < maximum and abs(lo) < maximum:
                                reaction.fluxBounds.backward = core.rangedNumber(lo,backwardVal,hi)
                        else:
                            reaction.fluxBounds.net.lo = max(0,reaction.fluxBounds.net.lo) 

                    else:
                        raise Exception('Wrong level')


    def getFluxDictionary(self, level=1, norm=1, fluxComp='all', rangeComp='best', deepCopy='True'):
        """Gets flux dictionary. 
        rangeComp decides whether to give out a single number or a ranged number. 
        norm can be a number, or the name of a flux in the network, to normalize everything with (e.g. GLCpts, or EX_glc_e_).
        This should be done by defining __rdiv__ properly
        fluxComp is the flux component to be output (e.g. forward,backward,net,exchange) in the level=1 case.
        copy decides wether the dictionary is a deepcopy or not. (The False option should be further developed.)
        """
        reactions = self.reactions

        # Dealing with normalization
        newNorm = norm
        if not utils.is_float(norm):
            for reaction in reactions:
                if reaction.name == norm:
                    newNorm = reaction.flux.getComp('net','best')
                    break
            if not utils.is_float(newNorm):
                # If we fail to match a reaction name from which to draw a flux,
                # or the flux component returns a non-numeric value like "N/A", fall back to a value of 1.
                newNorm = 1
        norm = abs(newNorm)
        if norm == 0:
            norm = 1    # Even if it's unlikely, we should make sure to avoid division by 0...

        fluxDict  = {}

        if deepCopy:
            # Main loop to get fluxes        
            for reaction in reactions:
                if level == 1:
                    r = reaction.flux.getComp(fluxComp,rangeComp)
                    if not isinstance(r, str):  # Might be "NA" or similar, in which case normalization doesn't apply
                        r = old_div(r, norm)
                    fluxDict[reaction.name] = r
                else:
                    if reaction.reversible and not reaction.exchange:
                        fluxDict[reaction.name+"_f"] = old_div(reaction.flux.getComp('forward' ,rangeComp), norm)
                        fluxDict[reaction.name+"_b"] = old_div(reaction.flux.getComp('backward',rangeComp), norm)
                    else:
                        fluxDict[reaction.name] = old_div(reaction.flux.getComp('net' ,rangeComp), norm)
        else:
        # TODO: fluxcomp, rangecomp options should be worked out for THIS ONE!!!        
            for reaction in reactions:
                fluxDict[reaction.name] = old_div(reaction.flux,norm)
        return fluxDict


    def addFluxes(self, fluxDictionary):
        fluxDict = self.findCommonFluxes(fluxDictionary)
        reactionDict = self.getReactionDictionary()
        for name in fluxDict:
            reactionDict[name].flux = fluxDict[name]


    def addFluxesWithDiffCheck(self, fluxDictionary):
        """
        Does the equivalent of addFluxes above, using a more direct method that does
        not combine forward and backward fluxes, then prints a series of debug
        statements if the results differ from the addFluxes method.
        """
        fluxDictA = {} ### out
        for reaction in self.reactions:
            if reaction.name in fluxDictionary:
                reaction.flux = fluxDictionary[reaction.name] 
                fluxDictA[reaction.name] = fluxDictionary[reaction.name] #### out                    

        fluxDictB = self.findCommonFluxes(fluxDictionary)

        #### Comparison out
        diff = set(fluxDictA.keys()) - set(fluxDictB.keys())      
        if diff:
            print(diff)
        else:
            for key in fluxDictA:
                if str(fluxDictA[key]) != str(fluxDictB[key]):
                    print(fluxDictA[key])
                    print(fluxDictB[key])  


    def findCommonFluxes(self, fluxDictionary):
        """
        Helper function for addFluxes and addFluxesWithDiffCheck.
        Returns a dictionary of flux objects for all the fluxes in fluxDictionary that are also present in
        this ReactionList's dictionary.  Combines forward and backward fluxes into single fluxes along the way.
        """
        reactionDict = self.getReactionDictionary()
        fluxDict = {}
        for name in fluxDictionary:
            fluxOrig = fluxDictionary[name]
            # Consolidate forward (X_f) and backward (X_b) fluxes into a single flux with name X.
            if name[-2:] == '_f':   #  level 2 input
                newName  = name[0:-2]
                # Deliberately not catching an exception here because backward reaction should always be present
                back     = fluxDictionary[newName+'_b']
                fluxOut  = core.flux(for_back_tup=(fluxOrig,back))
            else:                #  level 1 input
                newName  = name
                if isinstance(fluxOrig, core.flux):
                    fluxOut = fluxOrig
                else:
                    fluxOut = core.flux(net_exc_tup=(fluxOrig,'NA'))  
            # Add to reactions if present            
            if newName in reactionDict:
                fluxDict[newName] = fluxOut
        return fluxDict


    def addFluxBounds(self,fluxDictionary):   # not really used anywhere, eliminate?
        "Adds flux bounds, mirrors addFluxes"
        reactionDict = self.getReactionDictionary()
        for name in fluxDictionary:
            if name[-2:] == '_f':   #  level 2 input
                newName  = name[0:-2]
                forw     = fluxDictionary[name]
                # Deliberately not catching an exception here because backward reaction should always be present
                back     = fluxDictionary[newName+'_b']
                fluxOut  = core.flux(for_back_tup=(forw,back))
            else:                #  level 1 input
                newName  = name
                fluxOrig = fluxDictionary[newName]
                if isinstance(fluxOrig, core.flux):
                    fluxOut = fluxOrig
                else:
                    fluxOut = core.flux(net_exc_tup=(fluxOrig,'NA'))  
            # Add to reactions if present            
            if newName in reactionDict:
                reactionDict[newName].fluxBounds = fluxOut


    def getFluxBoundsDictionary(self):
        "Provides a dicionary of flux bounds"
        reactions = self.reactions
        fluxDict  = {}
        for reaction in reactions:
            fluxDict[reaction.name]=reaction.fluxBounds
        return fluxDict


    def capFluxBounds(self,maximum=100):
        # caps the maximum of flux bounds. ASSUMES NO DEEPCOPY!!!
        maximum = abs(maximum)
        boundDict = self.getFluxBoundsDictionary()
        for name in boundDict:
            bound = boundDict[name]
            bound.net.hi =  maximum if bound.net.hi>  maximum else bound.net.hi
            bound.net.lo = -maximum if bound.net.lo< -maximum else bound.net.lo    
        
        
    def getReactionDictionary(self):     # Do we really need both of this and getFluxDictionary?
        """
        Creates dictionary of reactions
        """
        reactions = self.reactions
        reactDict  = {}
        for reaction in reactions:
            reactDict[reaction.name]=reaction
        return reactDict


    def getTransitionLines(self):
        "Provides carbon transition lines"
        reactions = self.reactions
        transitionLines = []
        for reaction in reactions:
            try:
                transitionLines.append(reaction.transitionLine.fullLine)
            except AttributeError:
                pass    
        return transitionLines
    
    
    def writeReactionFile(self,filename):
        "Writes reaction file = reaction lines + options"
        # ## Options
        EMUNumberingStart = 1
        metList = self.getMetList()

        newlines = []

        #   NUMBERSTART
        newlines.append("&NUMBERSTART " + str(EMUNumberingStart))

        #   MLIST
        carbonDict = self.getCarbonDict()
        for EMUname in metList.getDestinationEMUs():
            oneEmu     = core.EMU(EMUname)
            numbers    = oneEmu.getIndices()
            metName    = oneEmu.getMetName()
            ncarbons   = carbonDict[metName]

            for n in range(ncarbons+1):
                if n not in numbers:
                    newlines.append('&MLIST '+metName+' '+str(-int(n)))                            

        #   EXCLUDE    
        for met in metList.getExcluded():
            newlines.append('&EXCLUDE '+met)    

        #   SOURCE
        for met in metList.getSource():
            newlines.append('&SOURCE '+met)

        #   VERBOSE    
        newlines.append('&VERBOSE')
        
        # ## Reaction lines
        transitionLines = self.getTransitionLines()
        newlines.extend(transitionLines)
        
        # ## Output
        file = open(filename,'w')
        for line in newlines:
            file.write(line+'\n')
        file.close()


    def getCarbonDict(self):  
        "Provides dictionary of number of carbons for each metabolite (from reaction line information)"
        metList = self.getMetList()
        carbonDict = metList.getCarbonDict()
        self.carbonDict = carbonDict

        return carbonDict


    def carbonTransitionsOK(self):
        "Checks if carbon transition information is ok (i.e. all reactions have the carbon transition info)"        
        reactions = self.reactions
        OK = True
        for reaction in reactions:   
            if not hasattr(reaction, 'transitionLine'):
                OK = False
        return OK
        # Need to check product and reaction names with reaction line names!!
        

    def processCarbonTransitions(self, metList='default'):    # TODO(Hector): pack better!!
        """Calculates all relevant information for the EMU reaction system
            New version which does not use Suthers' code. In progress"""
        debug = True        

        # ### Default MetList is the one included in ReactionList
        if metList == 'default':
            metList = self.getMetList()
        
        # ### Check that transitions make sense (all info is included)
        if not self.carbonTransitionsOK():
            raise Exception("Not all reactions have carbon transition information")
               
        # ### Main processing
        self.EMUSourceMetabs   = metList.getSource()
        self.EMUExcludedMetabs = metList.getExcluded()
        self.rawEMUReactions   = self.getTransitionLines()

        # Get carbon transition lines
        tLines = [core.AtomTransition(line) for line in self.rawEMUReactions]
        cTrans = AtomTransitionList(tLines)
        carbonDict = cTrans.getCarbonDict()
        
        # Integrate carbonDict information from atom transitions into metList if non-existent
        metList = self.getMetList()
        metList.updateNcarbons(carbonDict)
        
        # Get destination EMUs
        destEMUs          = metList.getDestinationEMUs()
        
        # Get all EMU transitions. TODO: make method in AtomTransitionList to get all EMU transitions
        transitions = []
        for emuName in destEMUs:
            transitions.extend(cTrans.findEMUTransitions2Source(core.EMU(emuName),self.EMUSourceMetabs,cleanCache=True))                
                
        transitions = sorted(list(set(transitions)),key=lambda trans: trans[1])
        EMUtransitions = EMUTransitionList(transitions)        
        
        if debug:
            EMUtransitions.write('VerboseFullEMUNetwork.txt')
        
        # Transition lines preprocessing
        EMUtransitions.preprocess(self.EMUExcludedMetabs) # TODO: rewrite preprocess
        
        self.transitions = EMUtransitions   

        # Equivalent emus and dummy emus data
        equivEMUs = EMUtransitions.equivEMUs
           
        # EMU list, names, metabolites and carbons
        emuList             = EMUtransitions.getEMUList(equivEMUs.elements)
        emuList.allEMU      = emuList.getAllEMUSet()
        emuList.metEMU      = emuList.getMetEMUPar()
        emuList.emuCarbons  = emuList.getEMUNCarbons()
        
        self.emuList = emuList
        

    def getCarbonTransitionReactions(self,inFluxRef):  
        "Provides only the reactions with carbon transition information and limits metabolites to those involved in carbon exchange"
        attribs = ['forward','backward','net','exchange']        
        inFluxRef = inFluxRef if inFluxRef !=0 else 1.0               
        self.inFluxRef = inFluxRef
        reactions = self.reactions

        # General metabolite dictionary to avoid repeated metabolites
        newReactions = []
        newMetaboliteList = copy.deepcopy(self.getMetList())
        newMetDict = {met.name: met for met in newMetaboliteList}
                
        for reaction in reactions:
            if hasattr(reaction, 'transitionLine'):               
                newReaction = copy.deepcopy(reaction)

                metsAll      = [x.name for x in reaction.transitionLine.products]
                metsAll.extend([x.name for x in reaction.transitionLine.reactants])
                
                newProducts = []
                newReacts   = []
                # Eliminate metabs from prods and reacs not in line
                for prod in newReaction.products:
                    if prod.name in metsAll: 
                        prod.originalMet = newMetDict[prod.name]
                        newProducts.append(prod)
                for react in newReaction.reactants:
                    if react.name in metsAll:
                        react.originalMet = newMetDict[react.name]
                        newReacts.append(react)
                
                newReaction.products  = newProducts
                newReaction.reactants = newReacts
                # Scale flux bounds by inFluxRef
                for attrib in attribs:
                    entry = getattr(newReaction.fluxBounds,attrib)
                    if isinstance(entry, core.rangedNumber):
                        entry.lo   = old_div(entry.lo, abs(inFluxRef))
                        entry.best = old_div(entry.best, abs(inFluxRef))
                        entry.hi   = old_div(entry.hi, abs(inFluxRef))
                # Add to new set of reactions
                newReactions.append(newReaction)
                
        newReactionsList = ReactionList(sorted(newReactions))
        
        return newReactionsList


    def getMetList(self):
        "Provides list of metabolites associated with these reactions"
        reactions = self.reactions
        allMets = []
        for reaction in reactions:
            for prod in reaction.products:
                allMets.append(prod.originalMet)   # Adding original metabolite, not product or reactant
            for react in reaction.reactants:
                allMets.append(react.originalMet)

        allMets = list(set(allMets))
        newMetList = MetaboliteList(allMets)
        return newMetList


    def getUB(self):
        "Provides upper bounds for reactions"
        ElementsUB = {}
        for reaction in self.reactions:
            ub = reaction.fluxBounds.net.hi
            ElementsUB[tuple([reaction.name])] = ub

        UBPar = GAMSclasses.GAMSParameter('ubvec',ElementsUB)
        return UBPar


    def getLB(self):
        "Provides lower bounds for reactions"
        ElementsLB = {}
        for reaction in self.reactions:
            lb = reaction.fluxBounds.net.lo
            ElementsLB[tuple([reaction.name])] = lb

        LBPar = GAMSclasses.GAMSParameter('lbvec',ElementsLB)
        return LBPar


    def getCvec(self):
        """
        Provides objective vector for linear problem of the type (FBA):
            max c_j*v_j subject to sum_j S_{ij}*v_j=0          
        """
        ElementsCvec = {}
        for reaction in self.reactions:
            if hasattr(reaction,'cval'):
                c = reaction.cval
                ElementsCvec[tuple([reaction.name])] = c
        cVecPar = GAMSclasses.GAMSParameter('cvec',ElementsCvec)
        return cVecPar


    def checkTransitionLines(self,metsCore):
       "Checks that reaction lines products and reactants are complete." 
       
       for reaction in self.reactions:
           if hasattr(reaction, 'transitionLine'):
               reaction.checkTransitionLineCompleteness(metsCore)


    def getReactionSubSet(self,rxnNames='default',cofactor='default',function='default'):
        "Gets a subset of reactions given by the rxn names"
        newReactionList = []            
        if isinstance(rxnNames,list):  #list of names    
            reactionDict = self.getReactionDictionary()
            for name in rxnNames:
                newReactionList.append(reactionDict[name])
        elif isinstance(rxnNames,str):  # string
            if rxnNames == 'default': # consider other input
                if cofactor == 'default':
                    if function != 'default':
                        for reaction in self:
                            if function(reaction):
                                newReactionList.append(reaction)
                    else:
                        raise Exception('Wrong input for ReactionList.getReactionSubSet')
                else:
                    for reaction in self:
                        metDict = reaction.getProdDict()
                        metDict.update(reaction.getReactDict())
                        if cofactor in metDict:
                            newReactionList.append(reaction)
            elif rxnNames == 'all':   # Give all reactions (duplication)
                for reaction in self:
                    newReactionList.append(reaction)    
            elif rxnNames == 'exchange':
                for reaction in self:
                    if reaction.exchange:
                        newReactionList.append(reaction)
            else:     # name of a subsystem
                for reaction in self:
                    if reaction.subsystem == rxnNames:
                        newReactionList.append(reaction)
        
        output = ReactionList(newReactionList,sort=False)
        
        return output


    def sort(self,function='default'):
        "Sorts reactions according to function"
        if function == 'default' or function=='by name':
            function = lambda reaction: reaction.name     # default is based on name
            revVal   = False
        elif function == 'by net flux':
            function = lambda reaction: reaction.flux.net.best
            revVal   = True
        elif function == 'by absolute net flux':
            function = lambda reaction: abs(reaction.flux.net.best)
            revVal   = True
        elif function == 'by absolute low flux':
            function = lambda reaction: abs(reaction.flux.net.lo)
            revVal   = True
        elif function == 'by absolute high flux':
            function = lambda reaction: abs(reaction.flux.net.hi)
            revVal   = True

        self.reactions = sorted(self.reactions,key=function,reverse=revVal)


    def extend(self,other):
        "Extends the list with the elements from other ReactionList"
        reactions      = self.reactions
        otherReactions = other.reactions
        
        reactions.extend(otherReactions)
        self = self.__init__(reactions,sort=False)        


    def euclidDist(self,other):
        "Finds mean euclidean distance with other ReactionList"
        # Check that lists are of the same lenght
        otherReactionDict = other.getReactionDictionary()
        if len(self.reactions) != len(other.reactions):
            raise Exception('Reactions lists are not the same length')
        # Find distance
        ED=0
        for reaction in self:
            val      = reaction.flux.getComp('net','best')
            otherVal = otherReactionDict[reaction.name].flux.getComp('net','best')
            delta    = 1
            ED       = ED + (old_div((val-otherVal),delta))**2
        ED = old_div(math.sqrt(ED),math.sqrt(len(self.reactions)))
        
        return ED


    def printFluxes(self,fileName='None',names='default',fluxDict='default',nonZero=False,brief=True):
        """
        Prints flux values.        
        nonZero decides whether to print fluxes with nonzero values
        """

        if fluxDict == 'default':
            Fluxes = self.getFluxDictionary(rangeComp='all')
        else:
            Fluxes = fluxDict    
            
        lines = []
        
        # Fluxes to print
        if names == 'default':
            names = list(Fluxes.keys())
        if names == 'exchange':
            tolerance = 0.0001
            names = []
            for name in list(Fluxes.keys()):
                if (('EX_' in name) or ('iomass' in name)) and (abs(Fluxes[name].net.best) > tolerance):
                    names.append(name)            
        
        # sort by quantity
        for flux in sorted(names, key=lambda name: abs(Fluxes[name].net.best),reverse=True):         
            if brief:
                lines.append(flux+": \t"+str(Fluxes[flux].net.best))
            else:                
                lines.append(flux+':\n')
                lines.append(Fluxes[flux].__str__())
                lines.append('\n')

        # Output
        if fileName == 'None':    # Print to screen
            for line in lines:
               print(line)
        elif fileName == 'toString':  # Print to string and return
            out = '\n'.join(lines)
            return out             
        else:                   # Print to file
            file = open(fileName,'w')
            for line in lines:
                file.write(line)
            file.close()    


    def getBiomassReaction(self):
        "Obtains appropriate biomass reaction"
        
        # Get all candidates        
        bmRxns = []
        for reaction in self.reactions:
            if 'biomass' in reaction.name or 'Biomass' in reaction.name:
                bmRxns.append(reaction)
        
        # Choose most appropriate 
        output = []
        if bmRxns:
            for rxn in bmRxns:
                if "core" in rxn.name or "Core" in rxn.name:
                    output = rxn
            if not output:
                output = bmRxns[0]   # Choose the first
        
        return output


    def barMaxMinFluxes(self,maxPlot=10,titleFig='',save='default'):
        "Produces a bar graph showing the minimum (lo) and maximum (hi) of fluxes"

        import matplotlib.pyplot as plt
        from pylab import figure, title, savefig

        # Produce data for bar
        names = []
        lower = []
        upper = []
        total = []
        for reaction in self:
            if len(names) < maxPlot:
                names.append(reaction.name.replace('EX_','',1).replace('_e_','',1).replace('(e)','',1).replace('_','-'))
                lower.append(reaction.flux.net.lo)
                upper.append(reaction.flux.net.hi-reaction.flux.net.lo)
                total.append(reaction.flux.net.hi)
        
        # bar data
        N = len(names)
        ind = numpy.arange(N)+0.5    # the x locations for the groups
        width = 0.45                 # the width of the bars
        
        fig = figure()
        axes = fig.add_subplot(111) 
        axes.axis([0,N+width,min(1.1*min(min(lower),0),min(total)),max(total)*1.1])

        p1 = plt.bar(ind, lower,  width, color='r', alpha=0.1)
        p2 = plt.bar(ind, upper, width, color='r', bottom=lower,alpha=0.8)

        plt.ylabel('Flux (mmol/gdw/h)')
        plt.title('Minimum and maximum values for exchange fluxes')
        plt.xticks(ind+old_div(width,2.),tuple(names),rotation='vertical')
        plt.legend( (p1[0], p2[0]), ('Minimum', 'Maximum') )

        title(titleFig)
        # Change x tick label size
        for label in axes.get_xticklabels():
                label.set_fontsize(16)
        # Inclination for xlabels
        fig.autofmt_xdate()
        # Save fig        
        if save != 'default':
            savefig(save) 


    def plotFluxes(self,color='k',fmt='o',markersize=10,names='default',axes='default',
                   xLabel='',yLabel='',titleFig='',label='label',plotType='default',axesRanges=[-30,30]):
        "Plots fluxes with confidence intervals in square axes"  # TODO(Hector): This can be improved

        import matplotlib.pyplot as plt
        from pylab import figure, title, xlabel, ylabel
        from matplotlib.patches import Polygon

        axisMax = axesRanges[1]  # = 10
        axisMin = axesRanges[0]  # = -5
       
        # default for names is all names
        if names == 'default':
            names = self.getReactionNameList(level=1,sort=False)
            
        # Do the plot
        reactionDict = self.getReactionDictionary()
        x  = []  
        y  = [] 
        DyU= []
        DyL= []
        counter = 1
        for name in names:
            reaction = reactionDict[name]
            y.append(reaction.flux.getComp('net','best'))
            DyU.append(reaction.flux.getComp('net','hi')-reaction.flux.getComp('net','best'))
            DyL.append(-(reaction.flux.getComp('net','lo')-reaction.flux.getComp('net','best')))
            x.append(counter)
            counter+=1
        x   = numpy.array(x)            
        y   = numpy.array(y)            
        DyU = numpy.array(DyU)            
        DyL = numpy.array(DyL)            

        if axes =='default':
            axes = figure(figsize=(8,8)).add_subplot(111)
        axes.axis([0,len(x)+1,max(axisMin,1.1*min(y-DyL)),min(axisMax,max(y+DyU)*1.1)])
        axes.axis([0,len(x)+1,axisMin,axisMax])

        if (len(x)!=0) and (len(y)!=0):
            if plotType == 'default':
                hand = plt.plot(x,y,color+fmt,markersize=markersize,label=label)    
            if plotType == 'errorbar':
                hand = axes.errorbar(x,y,yerr=[DyL,DyU],fmt=fmt,color=color,markersize=markersize,label=label)    
            elif plotType == 'transparent':
                hand = plt.plot(x,y,color+fmt+':',markersize=markersize,label=label,alpha=0.6)    
            elif plotType == 'transparent std':
                # Make polygon
                verts = list(zip(x,y+DyU)) + [i for i in reversed(list(zip(x,y-DyL)))] 
                poly  = Polygon(verts, facecolor=color, edgecolor=color,alpha=0.2)
                hand  = axes.add_patch(poly)

            xlabel(xLabel)
            ylabel(yLabel)
    
            title(titleFig)    
        else:
            print("No data in plot")
            print(x)
            print(y)        

        # return handle
        return hand


class EMUList(object):
    """
    Class for list of emus in a reaction network. Like in the case for MetaboliteList and ReactionList, this is just an emu list with some useful functions atached to it.
    """
    
    def __init__(self,emus):
        self.emus = emus

    def getMetEMUPar(self):
        "Produces metabolite emu relationship in the form of a GAMS parameter"
        metEMUDict = {}
        for emu in self.emus:
            if not emu.dummy:
                metEMUDict[tuple([emu.met,emu.name])] = 1.0
        metEMUPar = GAMSclasses.GAMSParameter('metemu',metEMUDict)   
        return metEMUPar


    def getEMUNCarbons(self):
        "Produces GAMS parameter holding the number of carbons for each emu" 
        nCarbonsDict = {}
        for emu in self.emus:
            nCarbonsDict[emu.name] = emu.ncarbons
        EMUNCarbonsPar = GAMSclasses.GAMSParameter('ncarbons',nCarbonsDict)
        return EMUNCarbonsPar


    def getAllEMUSet(self):
        "Get set of all emus in the form of a GAMS parameter" 
        emuSet = set([])
        for emu in self.emus:
            emuSet.add(emu.name)
        return GAMSclasses.GAMSSet('allemu',emuSet)        



class EMUTransitionList(object):
    """
    Class for list of EMU transitions in a reaction network. Like in the case for MetaboliteList and ReactionList, 
    this is just an EMU transition list with some useful functions atached to it.
    """

    def __init__(self,input_):
        if type(input_)==str: # Assume filename
            filename = input_
            # Open file and load transitions
            transfile   = open(filename)
            transitions = []
            for line in transfile:
                transitions.append(core.EMUTransition(line))      
            transfile.close()    
            # Store transitions
            self.transitions = transitions
        elif type(input_)==list: # Assume list of transitions
            # Test these are really EMU transitions            
            newInput = []
            for element in input_:
                if element.__class__.__name__=='EMUTransition':
                    newInput.append(element)
                elif element.__class__.__name__=='tuple':
                    newInput.append(core.EMUTransition(element))
                else:
                    raise Exception("List element: "+str(element)+" is not an EMU transition") 
            # Store transitions    
            self.transitions = newInput    
        else:
            raise Exception('Invalid input: '+str(input_))

        cond    = []
        nonCond = []
        self.condensation    = cond
        self.nonCondensation = nonCond


    def getCondNonCond(self):
        """
        Divides EMU transitions into condensation reactions (e.g.):

        PSCVT, pep_c_2_3 + skm5p_c_1_2_3_4_5_6_7 --> 3psme_c_1_2_3_4_5_6_7_8_9'
        
        and non-condensation reactions:
        
        'HEX1, glc_D_c_1_2_3_4_5_6 --> g6p_c_1_2_3_4_5_6'        
        
        """        
        transitions  = self.transitions
        cond         = []
        nonCond      = []
        for trans in transitions:
            if trans.cond:
                cond.append(trans)
            else:
                nonCond.append(trans)
        return cond,nonCond


    def preprocess(self,excMetSet):
        """
        Does preprocessing of condensation reactions. 
        """
        # TODO(Hector): check if still necessary
        # Creates new transition networks to account for condensation reactions
        transitions = self.transitions

        # Reg exps
        newTransitions     = []
        equivEMUsDict = {}
        dummyEMUsSet  = set([])		
        for trans in transitions:
            # If not a condensation reaction leave alone
            if not trans.cond:
                newTransitions.append(trans) 
            # Change things if it is a condensation reaction    
            else:
                # Keep old transition
                newTransitions.append(trans) 
                # Create new trans for dummy emus from the condensation reaction
                productEMU  = trans.products[0]
                # Exclude output if output metabolite is in the excluded metabolite set
                if productEMU.met not in excMetSet:
                    combined_reactant = ''.join([x.name for x in trans.reactants])
                    newTrans = core.EMUTransition(trans.name + ", " + combined_reactant + " --> " + productEMU.name + "\n")
                    newTransitions.append(newTrans)
                    equivEMUsDict[(combined_reactant,productEMU.name)] = 1
                    dummyEMUsSet.add(combined_reactant)
            
        eqEMUsPar = GAMSclasses.GAMSParameter('eqemu',equivEMUsDict)              # TODO: do this somewhere else
        dummyEMUsGAMSSet      = GAMSclasses.GAMSSet('dummyEMUs', dummyEMUsSet)
        self.transitions      = newTransitions
        self.equivEMUs        = eqEMUsPar
        self.dummyEMUs        = dummyEMUsGAMSSet            
        
        # Getting condensation and non-condensation reactions
        cond,nonCond = self.getCondNonCond()
        self.condensation    = cond
        self.nonCondensation = nonCond 
                

    def getEMUList(self,equivEMUsDict):
        """
        Obtains EMU list from transitions
        """
        transitions = self.transitions        
        
        emuSet =  set([])
        for trans in transitions:
            for react in trans.reactants:
                emuSet.add(react.name)
            for prod in trans.products:
                emuSet.add(prod.name)
                    
        # Storing in emuList class
        emus = []
        for emuName in emuSet:
            if tuple([emuName]) in equivEMUsDict:
                equivalent = equivEMUsDict[tuple([emuName])]
            else:
                equivalent = 'default' 
            emu = core.EMU(emuName, equivalent)
            emus.append(emu)
        return EMUList(emus)


    def getEMMmat(self):
        """
        Produces EMM matrix (see Suthers et al) in GAMS format for future use
        """
        transitions = self.transitions

        emmDict   = {}
        for trans in transitions:
            cont = trans.contribution if trans.contribution else 1.0        
            if not trans.cond:       # Do this only for non-condensing EMU transitions (others taken care of in getCondGAMSReacs)
                rxn = trans.name
                reactant = trans.reactants[0].name
                product  = trans.products[0].name
                emmDict[tuple([rxn,reactant,product])] = cont
                emmDict[tuple([rxn,product,reactant])] = -cont

        return GAMSclasses.GAMSParameter('EMM',emmDict)   #TODO: eventually do this somewhere else to decouple from GAMS


    def getCondGAMSReacs(self):  # TODO(Hector): this can use some further organization
        """
        Creates the condensation reactions for the GAMS code of the type:
        
        eqn_condensation_rxn_1at2..f('accoa_c_1_2oac_c_2','2') =e= f('accoa_c_1_2','1')*f('oac_c_2','1')+f('accoa_c_1_2','2')*f('oac_c_2','0');
        
        from condensation EMU reactions of the type:
        
        r1, accoa_c_1_2 + oac_c_2_3_4 --> cit_c_1_2_3_4_5
        
        IT ASSUMES ONLY TWO REACTANTS!!
        """

        ## Auxiliary functions
        def findCombinationTuples(k,maxLs):
            # Need to change only this if more than two reactants are needed
            output = []
            for L1 in range(maxLs[0]+1):
                L2 = k-L1
                if L2 >=0 and L2 <= maxLs[1]:
                    output.append((L1,L2))
            
            return output


        def construcTerm(nameValues): 
            return '*'.join(["f('"+name+"','"+str(value)+"')" for name,value in nameValues])


        def getNameValuesList(Kvals,reactants): 
            reacNames = [x.name for x in reactants]
            Nreacs = len(reacNames)
            nameValuesList = []
            for Kval in Kvals:
                nameValues = []
                for ind in range(Nreacs):
                    nameValues.append((reacNames[ind],Kval[ind]))
                nameValuesList.append(nameValues)
        
            return nameValuesList
        
        # Main Loop
        count    = 1
        eqns     = []
        eqnNames = []
        for trans in self.condensation:

            # Find parameters
            reactants = trans.reactants
            
            combinedReactant = ''.join([x.name for x in reactants])
            maxK  = sum([int(react.ncarbons) for react in reactants])
            maxLs = [int(react.ncarbons) for react in reactants]

    
            # Main loop
            for k in range(maxK+1):
                # Find convolution combinations
                Kvals = findCombinationTuples(k,maxLs)
                # List of reactant names interspersed with convolution combinations
                nameValuesList = getNameValuesList(Kvals,reactants)    
    
                # Final equation
                eqnName=  "eqn_condensation_rxn_"+str(count)+"at"+str(k)
                eqn    =  eqnName+"..f('"+combinedReactant+"','"+str(k)+"') =e= "
                eqn   += '+'.join([construcTerm(nameValue) for nameValue in nameValuesList])+';\n'
        
                eqns.append(eqn)
                eqnNames.append(eqnName+"\n")
            count +=1

        return eqns,eqnNames
        

    def eliminateDuplicates(self):
        "Eliminates duplicated EMU transitions"
        self.transitions = list(set(self.transitions))


    def compare(self,other): # TODO: overload '=' for this comparison?
        """
        Compares with othe EMUTransitionList to see if they are exaclty the same
        """
        # Write output in the same sorted order
        local = []
        for transition in self.transitions:
            #transition = EMUTransition(line)
            local.append(transition.printLine())
        local = sorted(local)
        
        new = []
        for transition in other.transitions:
            #transition = EMUTransition(line)
            new.append(transition.printLine())
        new = sorted(new)

        # Check lengths
        if len(local) != len(new):
            raise Exception('Different number of transitions: '+str(len(local))+' vs '+str(len(new)) )
        
        # Make comparison
        for i in range(len(local)):
            if local[i] != new[i]:
                errmsg = str(local[i]) + '\n vs \n' + str(new[i])
                raise Exception('Different Lines:\n'+errmsg)

        print("Both transition lists are the same")


    def write(self,output='default'):
        """
        Transition output
        """
        if output == 'default':   # Spit as method output
            out = self.transitions
        else:                     # write to file name = output
            file_ = open(output,'w')
            for trans in sorted(self.transitions):
                file_.write(str(trans).strip()+'\n')
            file_.close()
            out =  '' 
        
        return out



class AtomTransitionList(object):
    "Class for a list of AtomTransition objects"

    def __init__(self, input_):
        """
        Input can be in the form of:
            1-) A list of AtomTransition objects
            
            2-) A atomTransitions file of the type:
        
        '                
        &MLIST s7p 0
        &MLIST mal-L 0

        &SOURCE glc-D[e]

        # Input Reactions
        EX_glc(e)	glc-D[e] <==> glcDEx	abcdef : abcdef
        GLCt2	glc-D[e] --> glc-D	abcdef : abcdef
        HEX1	glc-D --> g6p	abcdef : abcdef
        .....
        '

        Here &MLIST determines the metabolite whose labeling data (MDV) is being fit.
        The number indicates the number of the carbon left out in the measured EMU.
        For example: 
            A-) s7p 0      indicates that labeling for the EMU s7p_1_2_3_4_5_6_7 is to be fit               
            B-) s7p -1 -2  indicates that labeling for the EMU s7p_3_4_5_6_7 is to be fit  
            
        &SOURCE determines the feed metabolite.
        
        """
        
        if type(input_) == str:  # Assuming file Name 
            # Opening file
            fileName = input_
            transitionLines, destMets, sources, labelExcluded = self.readTransitionFiles(fileName)
            self.EMUDestMetabList  = destMets
            self.EMUSourceMetabs   = sources
            self.EMUExcludedMetabs = labelExcluded

            # Obtaining transitions
            transitions = []
            for line in transitionLines:
                transition = core.AtomTransition(line)
                transitions.append(transition)
            self.transitions = transitions
            
        elif type(input_) == list:
            transitions = input_
            # Test they are really AtomTransition objects
            for trans in transitions:
                assert trans.__class__.__name__ == 'AtomTransition'
            self.transitions = transitions
        else:
            raise Exception('Input type must be a list of AtomTransition objects or an atomTransitions file: '+str(input_))

        # Storing madeIn dictionary for future use
        self.madeIn = self.getMadeInDictionary(altLine=True)  


    def readTransitionFiles(self,fileName):
        """
        Parses carbon transitions files of the type detailed above in __init__ (case 2)
        
        & Indicates options (as opposed to transitions)
        
        # Indicates comments
        
        """
        
        options         = []
        transitionLines = []
    
        # Read input lines from file
        fileName, string = utils.file2tuple(fileName)        
        inputLines = string.split('\n')
    
        # Separate transition lines from options input
        for line in inputLines:
            line = line.strip()
            if line:      # eliminate blanks
                if line[:1] not in '#&':
                    transitionLines.append(line.strip())
                elif line[:1] in '&':
                    options.append(line.strip('&').strip())
                # Lines starting with # are considered commments and ignored
        
        # Separate options into different outputs
        destMets      = []
        sources       = []  
        labelExcluded = []
        for line in options:
            if 'MLIST' in line.upper():  
                parts = line.replace('MLIST','').split()
                name    = parts[0]
                numbers = [abs(int(x)) for x in parts[1:]]
    
                out     = [name]
                out.extend([str(x) for x in numbers])
    
                destMets.append(out)
                
            elif 'SOURCE' in line.upper():
                sources.append(line.replace('SOURCE','').strip())
            elif 'EXCLUDE' in line.upper():
                labelExcluded.append(line.replace('EXCLUDED','').strip())
            elif 'NUMBERSTART' in line.upper():
                pass   # Deprecated option. Left here only for backwards compatibility
            elif 'VERBOSE' in line.upper():
                pass   # Deprecated option. Left here only for backwards compatibility
            else:
                poss = ['MLIST','SOURCE','NUMBERSTART','VERBOSE','EXCLUDE']
                raise Exception('Option: '+str(line) + ' unknown. Possibilities are '+','.join(poss))
    
        return transitionLines, destMets, sources, labelExcluded


    def getReactionNetwork(self, name='no Name'):
        """
        Produces 13C reaction network out of atom transitions
        Only available if options data (SOURCES and MLIST) are available (typically obtained from a
        data file input instead of a list of transitions)
        """

        # Test that options data are present
        errMsg = 'Destination and source metabolite info not present\n'
        Dir    =  dir(self)
        assert ('EMUDestMetabList' in Dir) and ('EMUSourceMetabs' in Dir) and ('EMUExcludedMetabs' in Dir), errMsg
        
        # Dictionary for Destination Metabolites
        DestMetDict = {}
        for metPack in self.EMUDestMetabList:
            metName  = metPack[0]
            destInfo = metPack[1:]
            met = core.Metabolite(metName)            
            if met.name in DestMetDict:
                DestMetDict[met.name].append(destInfo)
            else:
                DestMetDict[met.name]=[destInfo]   

        # Creating reactions
        maxVal = 100.00
        allReactions = []
        for transition in self.transitions:
            name       = transition.name
            reversible = transition.reversible
            if reversible:
                LB = -maxVal
                UB = maxVal
            else:
                LB = 0
                UB = maxVal
            flux_   = core.flux(('NA','NA'),(core.rangedNumber(LB,LB,UB),'NA'))    

            exchange = False   # This should be revisited

            # Metabolites
            reactants = transition.reactants
            products  = transition.products

            rxn = core.Reaction(name, reversible=reversible, flux='None', fluxBounds=flux_, transitionLine=transition.fullLine, exchange=exchange, reactants=reactants, products=products)
            allReactions.append(rxn)

        reactionList = ReactionList(allReactions)
        reactionList.eliminateDupMets()  

        # Deriving metabolites and source, excluded and destination information
        metList = reactionList.getMetList()
        for met in metList.mets:
            if met.name in DestMetDict:                # Destination 
                met.destination = True
                met.carbonsOut  = DestMetDict[met.name]
            if met.name in self.EMUSourceMetabs:            # Source 
                met.source = True
            if met.name in self.EMUExcludedMetabs:          # Excluded 
                met.excluded = True
                
        #  Storing carbon transitions type as model note
        notes = {}

       # Returning reaction network
        reactionNetwork = ReactionNetworks.C13ReactionNetwork( (name, notes, metList, reactionList) )
        return reactionNetwork
    
    
    def findEMUTransitions2Source(self,emu,sources,cleanCache=False):
        """
        Finds EMUtransitions all the way to the sources as a list in the format:
        (contribution, emuTransition)
        Where contribution is the flux contribution (e.g. 0.5 if symmetric molecules are present).
        For example:
        (0.5, 'ACGS, glu_L_c_2_3 --> acglu_c_2_3' )
        """
        madeIn = self.madeIn   
        if cleanCache:
            self.cache = []              

        # Main loop         
        EMUTransitions = []
        for transition in sorted(madeIn[emu.met], key=lambda x: str(x)): 
             # Unique index used to avoid processing the same transitions
            index = emu.name+' : '+'_'.join(str(transition))
            if index not in self.cache:
                self.cache.append(index)                       # Add index to cache to avoid processing previously processed transitions 
                EMUtransIn = transition.findEMUtransition(emu) # Add EMU transitions for this transition line  
                if EMUtransIn:                                 # Add only if transition exists
                    # Eliminate duplicates and add contributions       
                    tmpEMUTrans = list(set(EMUtransIn))  
                    nEMUTrans = len(tmpEMUTrans)
                    newEMUTrans = [(old_div(1.,nEMUTrans),EMUTrans) for EMUTrans in tmpEMUTrans]    
                    EMUTransitions.extend(newEMUTrans)
                        
                # Find origin EMUs and do the same for them
                EMUsNowAlts = transition.findEMUs(emu)      # This takes into account alternative transitions arising from symmetric molecules ...etc
                for EMUsNow in EMUsNowAlts:
                    for EMUNow in EMUsNow:
                        if (EMUNow.met not in sources):  # check met is not a carbon source
                            newTrans = self.findEMUTransitions2Source(EMUNow,sources,cleanCache=False)
                            EMUTransitions.extend(newTrans) 

        return EMUTransitions


    def getTransitionDict(self, level=1, altLine=False):  # TODO: level?
        """
        Obtains a dictionary of carbon transitions. 
        Level = 1 does not separate them into forward and backward versions, level = 2 does.
        If altLine = True, an alternate line notation is provided to deal with symmetric and repeated molecules molecules (see AtomTransition class)
        """

        transDict= {}
        for transition in copy.deepcopy(self.transitions):   # Making a copy of everything!! This probably slows things down 
            if 'EX_' not in transition.name:  # Eliminating exchange reactions. TODO(Hector): Do this more elegantly 
                if level==1:
                    transDict[transition.name] = core.AtomTransition(transition.printLine(altLine = altLine))  
                elif level==2:
                    if transition.reversible:
                        origName = transition.name
                        transition.name = origName + '_f'
                        transDict[transition.name] = core.AtomTransition(transition.printLine(altLine = altLine).replace('<==>','-->').replace(origName,transition.name))
                        transition.name = origName + '_b'
                        transDict[transition.name] = core.AtomTransition(transition.printReverse(altLine = altLine).replace('<==>','-->').replace(origName,transition.name))
                    else:
                        transDict[transition.name] = core.AtomTransition(transition.printLine(altLine = altLine))                      
                else:
                    raise Exception('Incorrect level: '+str(level))
        
        return transDict


    def getMadeInDictionary(self, altLine=False):
        """
        Obtains madeIn Dictionary with transitions where metabolite key is produced, e.g.:
            madeIn[metName] = transition        
        """
        madeIn = {}
        # Obtain dictionary of transitions
        transDict = self.getTransitionDict(level=2, altLine=altLine)

        # Add transitions for each product name
        for name in transDict:
            transition = transDict[name]

            for prod in transition.products:
                prodName = re.sub('(__ps\d+)','',prod.name)
                if prodName not in madeIn:
                    madeIn[prodName] = []   
                madeIn[prodName].append(transition)

        # Eliminate repetitions
        for name in madeIn:
            madeIn[name] = list(set(madeIn[name]))            

        return madeIn


    def getCarbonDict(self):
        """
        Obtains carbon dictionary for metabolites as reflected in the transitions
        """
        transDict = self.getTransitionDict(level=2)

        carbonDict = {}
        for name in transDict:
            transition = transDict[name]
            mets = []
            mets.extend(transition.products)
            mets.extend(transition.reactants)
            for met in mets:
                if met.name in carbonDict:
                    # Test if carbon numbers are consistent                    
                    newN = met.ncarbons
                    oldN = carbonDict[met.name]
                    if oldN != newN:
                        raise Exception('Metabolite %s has inconsistent carbon numbers: %s vs %s' % (name, str(newN), str(oldN)))  
                else:
                    carbonDict[met.name] = met.ncarbons
        return carbonDict

