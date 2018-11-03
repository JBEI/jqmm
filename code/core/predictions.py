# Part of the JBEI Quantitative Metabolic Modeling Library (JQMM)
# Copyright (c) 2016, The Regents of the University of California.
# For licensing details see "license.txt" and "legal.txt".

"""
The predictions module provides a couple of small classes and a function to predict fluxes from a reference flux set using 
Minimization of Metabolic Adjustment (MoMA, Segre et al PNAS 2002) or Regulatory On/Off Minimization (ROOM, Shlomi et al PNAS 2005) 
"""
from __future__ import division

import os
import FluxModels, GAMSclasses, core, ReactionNetworks
import utilities as utils

# Find directory of GAMS files
try: 
    dirGAMSFile =  os.environ['QUANTMODELPATH']+'/code/core/' 
except KeyError:
    raise Exception('QUANTMODELPATH environment variable not set')


class MOMAmodel(FluxModels.FluxModel):
    """Class for MOMA calculations."""
    
    def __init__(self,reactionNetwork,fluxPar):  
        """
        reactionNetwork: reaction network harboring the upper and lower bounds for the final (predicted) flux profile
        fluxPar: reference fluxes in the form of a GAMS parameter
        """
        
        # GAMS file information
        GAMSfileName     = "momaGAMS.gms"
        self.gamsFile    = dirGAMSFile + GAMSfileName
        
       # Storing SBML file
        self.reactionNetwork  = reactionNetwork   
       
        # Storing input files 
        GAMSInputFiles = self.reactionNetwork.getStoichMetRxnFBAFiles()
        vGEMini = fluxPar  
        
        vGEMiniFileName = 'vGEMini.txt'  
        vGEMiniString = vGEMini.write('toString')
        GAMSInputFiles.extend([(vGEMiniFileName,vGEMiniString)])

        self.GAMSInputFiles = GAMSInputFiles        

        
        # Loading output functions
        outputFuncs = {'Info':(FluxModels.infoOut,['OFGenInfo.txt']),
                       'Vout':(GAMSclasses.GAMSParameterFromResult,['MOMAout.txt','Vout']), 
                        }
        self.outputFuncs = outputFuncs
        
    def findFluxes(self,dirFinal='default', erase=True, probName='MOMAprob'):
        """Finds fluxes through MoMA
        
        dirFinal: directory for dumping files
        erase: erase all files when finished?
        probName: name of problem to be solved, MoMA or ROOM        
        """
        
        # Default dir for dumping files
        if dirFinal == 'default':
            dirFinal = os.getcwd()
        
        # Create Directory
        origDir     = os.getcwd()
        tmpDir = os.getcwd()+'/tmpDir/'
        try:
            os.mkdir(tmpDir)
        except OSError as xxx_todo_changeme:
            (errno,strerror) = xxx_todo_changeme.args
            if (errno ==17) and (strerror == 'File exists'):
                pass   # Directory already there
            else:
                raise 
                    
        # Get problem
        gamsFile = self.gamsFile
        reactionNetwork = self.reactionNetwork
        outputFuncs = self.outputFuncs
        Name = 'MOMAprob'
                
        GAMSprob       = GAMSclasses.GAMSProblem(Name,gamsFile,self.GAMSInputFiles,outputFuncs,execType='serial',erase=True)
        GAMSprob.writeAll()
        
        # Run problem(s)
        GAMSprob.run()

        # Check if they are done
        GAMSprob.waitTilDone()

        # Collect data
        GAMSprob.collect()
        
        # Output data 
        resultsDict = GAMSprob.Results 
        results     = MOMAResults(resultsDict,reactionNetwork)        
            
        # Return to original directory
        os.chdir(origDir)
        
        # Erase tmpdir
        if erase:
            utils.eraseFilesIn(tmpDir)
            GAMSprob.eraseFiles()

        return results
        
class ROOMmodel(MOMAmodel):
    "Class for ROOM calculations."
    
    def __init__(self,reactionNetwork,fluxPar):    
        """
        reactionNetwork: reaction network harboring the upper and lower bounds for the final (predicted) flux profile
        fluxPar: reference fluxes in the form of a GAMS parameter
        """
        # Reuse MoMA init
        MOMAmodel.__init__(self,reactionNetwork,fluxPar)        
        
        # Change GAMS file information
        GAMSfileName     = "roomGAMS.gms"
        self.gamsFile    = dirGAMSFile + GAMSfileName
        
        # Reload output functions
        outputFuncs = {'Info':(FluxModels.infoOut,['OFGenInfo.txt']),
                       'Vout':(GAMSclasses.GAMSParameterFromResult,['ROOMout.txt','Vout']), 
                        }
        self.outputFuncs = outputFuncs        
        
            
    def findFluxes(self,dirFinal='default',erase=True):
        
        return MOMAmodel.findFluxes(self,dirFinal=dirFinal,erase=erase,probName='ROOMprob')        
        
        
class MOMAResults(FluxModels.Results):
    """Class for holding MOMA results """    
    
    def __init__(self,resultsDict,reactionNetwork):
        """
        resultsDict: dictionary with flux results
        reactionNetwork: reaction network harboring the upper and lower bounds for the final (predicted) flux profile
        """
        
        self.reactionNetwork = reactionNetwork
        self.resultsDict = resultsDict
        self.processResults()
        
    def processResults(self):
        Vout             = self.resultsDict['Vout']
        self.Vout        = GAMSclasses.fluxGAMSpar(Vout.name,Vout.elements)
        
        # Ranged flux dict (PROVISIONAL: it should be changed)
        fluxDict = self.Vout.getReactionList().getFluxDictionary()
        rangedFluxDictAll = {}
        for name in fluxDict:
            flux = fluxDict[name]
            rangedForw = core.rangedNumber(flux.forward ,flux.forward ,flux.forward) 
            rangedBack = core.rangedNumber(flux.backward,flux.backward,flux.backward) 
            rangedNet  = core.rangedNumber(flux.net     ,flux.net     ,flux.net) 
            rangedExch = core.rangedNumber(flux.exchange,flux.exchange,flux.exchange) 
            
            rangedFlux = core.flux((rangedForw,rangedBack),(rangedNet,rangedExch))
        
            rangedFluxDictAll[name] = rangedFlux
        
        self.rangedFluxDict = rangedFluxDictAll
  
  
# Function which implements MOMA and ROOM
def predict(results,geneName,method='MOMA',sbmlFileName='default'):
    """
    This function takes a reference flux profile and produces the flux prediction for a gene knocked out
    
    results : results instance containing the reference flux profiles
    geneName: name of gene to be knocked out
    method  : method to be used, MoMA or ROOM
    sbmlFileName: sbml input harboring the upper and lower bounds for the final (predicted) flux profile    
    """
    
    # Input of final network as sbml
    reactionNetwork = ReactionNetworks.TSReactionNetwork(sbmlFileName)        
    # Knocking genes out
    if geneName.__class__.__name__ == 'str':
        reactionNetwork.changeFluxBounds(geneName     ,core.fluxBounds( 0, 0  ,False)[1])     
    elif geneName.__class__.__name__ == 'list':   # list of genes to be knocked out
        for gene in geneName:
            reactionNetwork.changeFluxBounds(gene     ,core.fluxBounds( 0, 0  ,False)[1])     
                
    # Obtain reference fluxes
    predResults = []
    if results.__class__.__name__ == 'FBAResults':      # if reference fluxes are FBA fluxes
        fluxList = [results.VFBA.getReactionList().getFluxDictionary()]           
    else:                                                # else assume TS reference fluxes
        fluxList = results.fluxDictEnsemble.fluxDictList

    # Make predictions
    for targetFlux in fluxList:
        fluxPar = GAMSclasses.fluxGAMSpar('vGEMini',[],fluxDict = targetFlux)
        if method == 'MOMA':
            model   = MOMAmodel(reactionNetwork,fluxPar)              
        elif method == 'ROOM':
            model   = ROOMmodel(reactionNetwork,fluxPar)                          
        predResults.append(model.findFluxes(erase=True))
        
    # Store predicted fluxes in dictionary
    fluxDictList = []
    for result in predResults:
        fluxDict = result.Vout.getReactionList().getFluxDictionary()
        fluxDictList.append(fluxDict)
    fluxNames   = list(predResults[0].Vout.getReactionList().getFluxDictionary().keys()) 
    fluxDictRef = fluxDictList[0]
    
    fluxDictEnsemble = FluxModels.fluxDictEns(fluxDictList,fluxNames,fluxDictRef=fluxDictList[0])
    attribs = ['forward','backward','net','exchange']
    fluxesVarHi,fluxesVarLo  = fluxDictEnsemble.getFluxDictStdRef(attribs)
    
    ranged = {}
    rangedFluxDictAll  = {}
    for flux in sorted(fluxDictRef):    
        for attrib in attribs:
            ranged[attrib] = core.rangedNumber(fluxesVarLo[attrib][flux], getattr(fluxDictRef[flux],attrib) , fluxesVarHi[attrib][flux])
        rangedFlux = core.flux((ranged['forward'],ranged['backward']),(ranged['net'],ranged['exchange']))
        rangedFluxDictAll[flux] = rangedFlux     


    # Creating model for prediction output
    TSModel = FluxModels.TwoSC13Model(sbmlFileName)   
    
    TSModel.name = method
    TSModel.reactionNetwork.reactionList.addFluxes(rangedFluxDictAll)

    # Prediction for forward and backward fluxes is to keep exhange flux as it was 
    inFluxRef = TSModel.reactionNetwork.getFluxRefScale(bounds=False, inFluxRefs= ['EX_glc(e)','EX_glc_e_'])
    
    if hasattr(results.reactionNetwork,'C13ReacNet'): 
        fluxDict = results.reactionNetwork.C13ReacNet.reactionList.getFluxDictionary(rangeComp='all')    
    else:
        fluxDict = results.reactionNetwork.reactionList.getFluxDictionary(rangeComp='all')  
    
    for reaction in TSModel.reactionNetwork.C13ReacNet.reactionList:
        net           = rangedFluxDictAll[reaction.name].net*(utils.old_div(1,inFluxRef))
        exch          = fluxDict[reaction.name].exchange
        reaction.flux = core.flux(net_exc_tup=(net,exch))  
         
    return TSModel
    
    