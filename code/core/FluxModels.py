# Part of the JBEI Quantitative Metabolic Modeling Library (JQMM)
# Copyright (c) 2016, The Regents of the University of California.
# For licensing details see "license.txt" and "legal.txt".
"""
The FluxModels module provides the classes needed to hold the three types of flux models used in jQMM: Flux Balance Analysis (FBA), 
13C Metabolic Flux Analysis (13C MFA) and two-scale 13C Metabolic Flux Analysis (2S-13C MFA). 
"""

import os, time, copy, numpy, re, matplotlib, unittest
import GAMSclasses, ReactionNetworks, core, DB, enhancedLists, labeling
from pylab import figure, xlabel, ylabel, title, savefig
from FluxMaps import FluxMap
import utilities as utils

# Specifying font size for plots
# TODO: move all plotting into separate classes or functions
font = {'size'   : 20}
matplotlib.rc('font', **font)



class FluxModel:
    """Base class to derive all others model classes involved in flux analysis """

    def __init__(self, reactionNetwork, GAMSfileName, outputFuncs):
        """
        ReactNet:       reaction nework used for modeling
        GAMSfileName:   the name of the appropriate GAMS file
        outputFuncs:    Functions used to process the raw output (files from GAMS, etc) into the processed output (e.g. fluxes and labeling patterns)
        """
        # Obtain GAMS file
        try: 
            self.dirGAMSFile =  os.environ['QUANTMODELPATH']+'/code/core/' 
        except KeyError:
            raise Exception('QUANTMODELPATH environment variable not set')
        gamsFile        = self.dirGAMSFile + GAMSfileName
        self.gamsFile   = gamsFile
        
        # Storing reactionNetwork
        self.reactionNetwork = reactionNetwork

        # Loading output functions
        self.outputFuncs = outputFuncs


    # TODO: eliminate by passing it to ReactionList methods
    def recordFluxes(self,fluxDict):
        """Records fluxes in the flux Dict in the reactions"""

        for reaction in self.reactionNetwork.reactionList:
            if reaction.name in fluxDict:
                reaction.flux = fluxDict[reaction.name]
           
       
        
    def runGAMSBatch(self,problemBatch):
        "General function for running a GAMS batch"
        
        with problemBatch as batch:
            # Run problem(s)
            batch.run()

            # Check if they are done
            batch.waitTilDone()
        
            # Collect data
            batch.collect()
        
            # Put results in local format 
            resultsDict = batch.getAllResults()
        return resultsDict          


    def runGAMSProblem(self,problem):
        "General function for solving a GAMS problem"
        
        # Using the 'with' statement to provide automatic setup and wrapup
        with problem as prob:
            # Run problem
            prob.run()

            # Check if it is done
            prob.waitTilDone()

            # Collect data
            prob.collect()
        
            # Output data 
            resultsDict = prob.Results            

        return resultsDict        
        
        
    def solveProblem(self,problem):
        """ Solves problem stored in self.problem. This self.problem can be a single GAMS problem or a GAMS problem batch"""
        # Determine what kind of problem we have and solve it
        self.problem = problem
        if problem.__class__.__name__== 'GAMSProblem':
            resultsDict = self.runGAMSProblem(problem)
        elif problem.__class__.__name__== 'GAMSProblemBatch':
            resultsDict = self.runGAMSBatch(problem)
            
        return resultsDict


class FBAModel(FluxModel):
    "Class for FBA calculations, sbmlFileName can be name of a file or file in string format."

    # TODO: redo init by leveraging FluxModel init as much as possible:
    #        FluxModel.__init__(self,reactionNetwork,GAMSfileName,outputFuncs)
    def __init__(self,sbmlFileName):
        "sbmlFileName can be the name of the sbml file or a tuple of the form (fileName,string)"
        # Obtain GAMS file information
        try: 
            dirGAMSFile =  os.environ['QUANTMODELPATH']+'/code/core/' 
        except KeyError:
            raise Exception('QUANTMODELPATH environment variable not set')
        GAMSfileName     = "GAMS_EMU_MFA_FLUX_FBA.gms"
        GAMSFVAfileName  = "GAMS_EMU_MFA_FLUX_FVA.gms"
        gamsFile         = dirGAMSFile + GAMSfileName
        gamsFileFVA      = dirGAMSFile + GAMSFVAfileName
        self.gamsFile    = gamsFile
        self.gamsFileFVA = gamsFileFVA

        # Reading SBML file
        self.reactionNetwork = ReactionNetworks.ReactionNetwork(sbmlFileName)

        # Loading output functions
        outputFuncs = {'Info':(infoOut,['OFGenInfo.txt']),
                       'Vout':(GAMSclasses.GAMSParameterFromResult,['VFBAout.txt','Vout']), 
                       'Successful':(getSuccessful,[GAMSfileName.replace('.gms','.lst',1)]),
                        }
        self.outputFuncs = outputFuncs
        outputFuncsFVA = {'Info':(infoOut,['OFGenInfo.txt']),
                          'Vout':(GAMSclasses.GAMSParameterFromResult,['VFBAout.txt','Vout']), 
                          'Vmax':(GAMSclasses.GAMSParameterFromResult,['Vmaxout.txt','Vmax']), 
                          'Vmin':(GAMSclasses.GAMSParameterFromResult,['Vminout.txt','Vmin']), 
                        }
        self.outputFuncsFVA = outputFuncsFVA


    # TODO Convert to getGAMSInputFilesFBA as soon as changeFluxBoundsFAST is gone
    def getGAMSInputFiles(self):
        """Provides GAMS files for FBA"""
        return self.reactionNetwork.getStoichMetRxnFBAFiles()


    def findFluxes(self,erase=True,warnFail=True): 
        """Finds fluxes through FBA"""
        # Get files for optimization problem
        gamsFile            = self.gamsFile
        outputFuncs         = self.outputFuncs
        self.GAMSInputFiles = self.getGAMSInputFiles()  
        Name                = 'FBAprob' 
               
        # Get GAMS optimization problem
        GAMSprob       = GAMSclasses.GAMSProblem(Name,gamsFile,self.GAMSInputFiles,outputFuncs,execType='serial',erase=erase)        
        
        # Solve GAMS problem
        resultsDict = self.solveProblem(GAMSprob)
        
        # Process raw results
        results     = FBAResults(resultsDict,self.reactionNetwork)   
        
        # Warn if LP failed
        if warnFail:
            results.printSuccess()
            
        return results


    # TODO: it would be nice to have a parallel version of this at some point
    def FVA(self,reactions='default',erase=True):
        """Finds fluxes through Flux Variability Analysis"""
                  
        # Get files for optimization problem
        gamsFile        = self.gamsFileFVA
        reactionNetwork = self.reactionNetwork
        outputFuncs     = self.outputFuncsFVA
        Name            = 'FBAprob'
        reactions       = set(reactions)

        # Get GAMS optimization problem
        GAMSInputFilesFVA = reactionNetwork.getStoichMetRxnFVAFiles(studiedReactions = reactions)
        GAMSprob       = GAMSclasses.GAMSProblem(Name,gamsFile,GAMSInputFilesFVA,outputFuncs,erase=erase)

        # Solve GAMS problem
        resultsDict = self.solveProblem(GAMSprob)
            
        # Process raw results
        results     = FVAResults(resultsDict,self.reactionNetwork)        

        return results


    # TODO(Hector): these two methods may not be needed anymore. Confirm and eliminate.
    def changeFluxBoundsFAST(self, reactionName, flux):
        "Changes flux bounds for problem without having to change whole sbml file (takes too long). This could be done in a more elegant way."
        
        fileNameUB   = 'ubvec.txt'
        fileNameLB   = 'lbvec.txt'
      
        reactionList = self.reactionNetwork.reactionList
        reactDict = reactionList.getReactionDictionary()
        reactDict[reactionName].fluxBounds = flux

        UBPar = reactionList.getUB()
        LBPar = reactionList.getLB()
              
        stUB = UBPar.write('toString')
        stLB = LBPar.write('toString')
        
        self.GAMSInputFiles[3] = (fileNameUB,stUB)
        self.GAMSInputFiles[4] = (fileNameLB,stLB)
        
        # Keeping changes in a list for later use with recordFASTchanges
        if hasattr(self,'FASTchanges'):
            self.FASTchanges[reactionName] = flux
        else:
            FASTchanges = {reactionName : flux}
            self.FASTchanges = FASTchanges


    def recordFASTchanges(self):
        "Permanently records in the sbml file and tuple the changes from changeFluxBoundsFAST"
        # TODO(Hector): This may be slowing things down.  See if it can be eliminated.
        for reactionName in self.FASTchanges:    # This for loop should go
            # Change Tuple
            reacDict = self.reactionNetwork.reactionList.getReactionDictionary()
            reaction = reacDict[reactionName]
            reaction.fluxBounds = self.FASTchanges[reactionName]
            
        # Updating tuple
        self.reactionNetwork = ReactionNetworks.ReactionNetwork( (self.reactionNetwork.name, self.reactionNetwork.notes, self.reactionNetwork.metList, self.reactionNetwork.reactionList))


    def changeObjective(self, reactionName, cval):
        "Changes objective value."
        reactionList = self.reactionNetwork.reactionList
        reactDict = reactionList.getReactionDictionary()
        reactDict[reactionName].cval = cval



class C13Model(FluxModel):
    "Class for 13C MFA calculations"

    def __init__(self,sbmlFileName):
        "sbmlFileName can be the name of the sbml file or a tuple of the form (fileName,string)"
        # GAMS file name
        GAMSfileName    = "GAMS_EMU_MFA_FLUX.gms"           
        infeasFileName  = GAMSfileName.replace('gms','lst')     
        
        # sbml file
        reactionNetwork = ReactionNetworks.C13ReactionNetwork(sbmlFileName)


        # Loading output functions
        outputFuncs = {'Info'           :(infoOut           ,['OFGenInfo.txt']),
                       'Vout'           :(fluxParOut        ,['Vout.txt','Vout']), 
                       'fout'           :(GAMSclasses.GAMSParameterFromResult,['fout.txt','fout']),
                       'labelcomp'      :(labelingParOut    ,['labelcomp.txt','labelcomp']),
                       'infeasibilities':(getInfeasibilities,[infeasFileName])    
                        }
                        
        # Keeping names for needed functions
        self.getFilesFunctionName    = 'get13CMFAfiles'
        self.resultsFunctionNameRand = 'C13Results'
        self.resultsFunctionNameVar  = 'C13FVAResults'
        
        FluxModel.__init__(self,reactionNetwork,GAMSfileName,outputFuncs)


    def getGAMSInputFiles(self,reactionNetwork,labelMin=0.99,vSuggested='default',procString='raw'):
        """ 
        Creates input Files for the GAMS problem.
        I separated this from createBatch for flexibility.
        """

        GAMSInputFiles = reactionNetwork.get13CMFAfiles(vSuggested=vSuggested, labelCompNormMin=labelMin, procString=procString)       
        
        return GAMSInputFiles      


    # TODO(Hector): Take maxFlux13C and labelMin to GAMS file class?
    def createBatchRand(self,Nrep,Nrand,labelMin,maxFlux13C,erase,procString='raw'):
        "Create Batch of gams problems"
        
        reactionNetwork = self.reactionNetwork
        reactionNetworkOrig = copy.deepcopy(reactionNetwork)
        gamsFile = self.getGAMSFileRand(maxFlux13C)
        outputFuncs = self.getOutputFuncsRand()                
        self.erase = erase

        # Creating problem batch
        problems = []
        for i in range(Nrand+1):
            for j in range(1,Nrep+1):
                Name    = '/Run'+str(i)+'_'+str(j)
                dirName = os.getcwd()+Name+'-'+GAMSclasses.GAMSProblem.getGAMSTempFolderName()
                GAMSclasses.GAMSProblem.getGAMSTempFolderName()
                GAMSInputFiles = self.getGAMSInputFiles(reactionNetwork,labelMin,procString=procString) # Don't take this out of loop to preserve unit tests
                GAMSprob       = GAMSclasses.GAMSProblem(Name,gamsFile,GAMSInputFiles,outputFuncs,directory = dirName,erase=erase)
                
                problems.append(GAMSprob)
                    
            reactionNetwork = copy.deepcopy(reactionNetworkOrig)
            reactionNetwork.randomizeLabeling()

        batch = GAMSclasses.GAMSProblemBatch(problems,erase=erase)
    
        return batch


    def createBatchVar(self,Nrep,fluxNames,labelMin,maxFlux13C,erase,procString='raw'):
        """Create Batch of gams problems """
        self.erase = erase

        # Find best solution through findFluxesStds and save as suggested fluxes
        vSuggested,resultsOne = self.getSuggestedFluxes(Nrep, erase, procString)

        # Record fits 
        self.recordFits(resultsOne.EMUlabel)
        
        # Change MDV stds to match best fit
        self.reactionNetwork.changeStds('toFits')     
        
        # GAMS file
        gamsFileVar = self.getGAMSFileVar(maxFlux13C)  

        # Loading output functions
        outputFuncsVar = self.getOutputFuncsVar()        
        
        # Creating gams problems and storing them in a batch class    
        problems = []
        GAMSInputFilesCommon = self.getGAMSInputFiles(self.reactionNetwork,labelMin,vSuggested=vSuggested,procString=procString)            
        for fluxName in fluxNames:
            # Changing reaction name to sbml
            tmpRxn  = core.Reaction(fluxName)
            fluxNameNew = tmpRxn.name            
            Name    = '/Run_'+fluxNameNew
            dirName = os.getcwd()+Name
                                      
            # Add file for flux to focus on
            fluxFocus         = GAMSclasses.GAMSSet('fluxFocus',set([fluxName]))
            fluxFocusSt       = fluxFocus.write('toString')
            fluxFocusFileName = 'fluxFocus.txt' 
            
            GAMSInputFiles = copy.deepcopy(GAMSInputFilesCommon)
            GAMSInputFiles.extend([(fluxFocusFileName,fluxFocusSt)])
        
            # Creating GAMS problem
            GAMSprob  = GAMSclasses.GAMSProblem(Name,gamsFileVar,GAMSInputFiles,outputFuncsVar,erase=erase)
            GAMSprob.dir   = dirName
            problems.append(GAMSprob)         
        
        batch = GAMSclasses.GAMSProblemBatch(problems)
        
        return batch,resultsOne


    def getSuggestedFluxes(self,Nrep,erase,procString):
        """Obtains suggested fluxes for createBatchVar"""

       # Find best solution through findFluxesStds and save as suggested fluxes
        resultsOne        = self.findFluxesStds(Nrep=Nrep,Nrand=0,erase=erase,procString=procString)    
            
        # Adding previously calculated flux profile starting point for range procedure         
        self.reactionNetwork.reactionList.addFluxesWithDiffCheck(resultsOne.rangedFluxDict)
        fitFluxes    = GAMSclasses.fluxGAMSpar('vSUGG',[],fluxDict=self.reactionNetwork.reactionList.getFluxDictionary(level=2,rangeComp='best'))
        
        vSuggested = tuple([fitFluxes])
        
        return vSuggested,resultsOne 


    def GAMSFileModifications(self,gamsFile,maxFlux13C,maxTime):
        "Does standard modifications for GAMS file"
        
        biomassReaction = self.reactionNetwork.reactionList.getBiomassReaction()
        if biomassReaction:            
            changes = [['BiomassEcoli',biomassReaction.name]]   # Change biomass reaction
        else: 
            changes = []
        
        if maxFlux13C:
            #changes.append(['maxflow13CU',maxFlux13C])
            changes.append(['maxflow13CU = (\d+)','maxflow13CU = '+str(maxFlux13C)])
            
        if maxTime:
            changes.append(['reslim = (\d+)','reslim = '+str(maxTime)])

       
        gamsFile = utils.substituteInFile(gamsFile,changes) 
        
        return gamsFile


    def getGAMSFileRand(self,maxFlux13C,maxTime=[]):
        "Returns GAMS file for the ranged fluxes approach"
        gamsFileRand  = utils.file2tuple(self.dirGAMSFile+'GAMS_EMU_MFA_FLUX.gms')        
        gamsFileRand  = self.GAMSFileModifications(gamsFileRand,maxFlux13C,maxTime)      
        
        return gamsFileRand


    def getGAMSFileVar(self,maxFlux13C,maxTime=[]):
        "Returns GAMS file for the ranged fluxes approach"
        gamsFileVar  = utils.file2tuple(self.dirGAMSFile+'GAMS_EMU_MFA_FLUXvar.gms')        
        gamsFileVar  = self.GAMSFileModifications(gamsFileVar,maxFlux13C,maxTime)      
        
        return gamsFileVar


    def getOutputFuncsRand(self):
        """Returns output functions for randomization  approach"""
        
        # GAMS file name
        GAMSfileName    = self.getGAMSFileRand([])[0]
        infeasFileName  = GAMSfileName.replace('gms','lst')     

        # Loading output functions
        outputFuncs = {'Info'           :(infoOut           ,['OFGenInfo.txt']),
                       'Vout'           :(fluxParOut        ,['Vout.txt','Vout']), 
                       'fout'           :(GAMSclasses.GAMSParameterFromResult,['fout.txt','fout']),
                       'labelcomp'      :(labelingParOut    ,['labelcomp.txt','labelcomp']),
                       'infeasibilities':(getInfeasibilities,[infeasFileName])    
                        }
                                
        return outputFuncs        


    def getOutputFuncsVar(self):
        """Returns output functions for ranged fluxes approach"""
        
        # Loading output functions
        GAMSfileName    = self.getGAMSFileVar([])[0]
        infeasFileName  = GAMSfileName.replace('gms','lst')   
        outputFuncsVar = {
                       'Info'           :(infoOut           ,['OFGenInfo.txt']),      # Temporary. This should be set in init
                       'InfoMax'        :(infoOut           ,['OFGenInfoMax.txt']),      
                       'InfoMin'        :(infoOut           ,['OFGenInfoMin.txt']),  
                       'Vout'           :(fluxParOut        ,['Vout.txt','Vout']), 
                       'Vmaxout'        :(fluxParOut        ,['Vmaxout.txt','Vmaxout']), 
                       'Vminout'        :(fluxParOut        ,['Vminout.txt','Vminout']), 
                       'fout'           :(GAMSclasses.GAMSParameterFromResult,['fout.txt','fout']),
                       'labelcomp'      :(labelingParOut    ,['labelcomp.txt','labelcomp']),
                       'labelcompMax'   :(labelingParOut    ,['labelcompMax.txt','labelcompMax']),
                       'labelcompMin'   :(labelingParOut    ,['labelcompMin.txt','labelcompMin']),
                       'infeasibilities':(getInfeasibilities,[infeasFileName]),
                       'fluxBounds':(GAMSclasses.GAMSParameterFromResult,['fluxBounds.txt','fluxBounds'])        
                        }
                        
        return outputFuncsVar


    def getResultsStds(self,resultsDict):
        "Produces results out of the results dictionary"        
        results     = C13Results(resultsDict,self.reactionNetwork)           

        return results


    def getResultsRanges(self,resultsDict,resultsOne):
        "Produces results out of the results dictionary"        
        resultsVar     = C13FVAResults(resultsDict,self.reactionNetwork)
        rangedFluxDict = resultsOne.rangedFluxDict
        
        # Use the resultsOne solution (best fit) as a basis
        results = copy.deepcopy(resultsVar)
        results.rangedFluxDict = resultsOne.rangedFluxDict
        results.EMUlabel = resultsOne.EMUlabel
        results.OF = resultsOne.OF

        # Use the results of the 13C FVA to obtain confidence intervals
        for name in rangedFluxDict:
            if name in resultsVar.rangedFluxDict:
                rangedFluxDict[name].net.hi   = resultsVar.rangedFluxDict[name][1] 
                rangedFluxDict[name].net.lo   = resultsVar.rangedFluxDict[name][0]
                
                rangedFluxDict[name].forward  = 'NA'
                rangedFluxDict[name].backward = 'NA'
                rangedFluxDict[name].exchange = 'NA'
         
        # TODO(Hector): there is probably no reason to store this here anymore. Eliminate! 
        results.rangedFluxDict = rangedFluxDict
        
        # Storing result in reaction list
        results.reactionNetwork.reactionList.addFluxesWithDiffCheck(results.rangedFluxDict)       

        return results


    def findFluxesStds(self,Nrep=20,Nrand=5,erase=True,labelMin=0.99,maxFlux13C=[],procString='raw'):
        "Calculates fluxes and calculates confidence intervals through the monte carlo approach"        
        # Create directory and move into it
        self.origDir     = os.getcwd()
        self.baseDirName = time.strftime("%Y%m%dT%H%M%SGAMSstd")
        os.mkdir(self.baseDirName)
        os.chdir(self.baseDirName)
        
        # Get batch
        batch = self.createBatchRand(Nrep,Nrand,labelMin,maxFlux13C,erase,procString=procString)

        # Solve batch problem
        resultsDict = self.solveProblem(batch)
        
        # Process raw results into final results
        results     = self.getResultsStds(resultsDict)
                
        # Delete Directory. There should be no files left inside
        os.chdir('..')
        if erase:
            os.rmdir(self.baseDirName)

        return results            


    def findFluxesRanges(self,Nrep,fluxNames,erase=True,labelMin=0.99,maxFlux13C=[],procString='raw'):
        "Finds fluxes and ranges of fluxes within allowable Objective Function"
        baseDirName = time.strftime("%Y%m%dT%H%M%SGAMSrng")
        os.mkdir(baseDirName)
        os.chdir(baseDirName)
        
        oldReactionNetwork = self.reactionNetwork 
         # Get batch
        batch,resultsOne = self.createBatchVar(Nrep,fluxNames,labelMin,maxFlux13C,erase,procString=procString)

        # Solve batch problem
        resultsDict = self.solveProblem(batch)
        
        # Process raw results into final results
        results  = self.getResultsRanges(resultsDict,resultsOne)
        
        # Going back to initial reactionNetwork (with original stds)
        self.reactionNetwork = oldReactionNetwork
        #self.reactionNetwork.changeStds(oldReactionNetwork.fragDict())   # this should be the final version

        # Delete Directory. There should be no files left inside
        os.chdir('..')
        if erase:
            os.rmdir(baseDirName)
            
        return results       


    def recordFits(self,EMUlabel):
        "Records fits to experimental data in SBML file"
        Notes        = self.reactionNetwork.notes       
        
        regexp = re.compile('(\D+)(\d+)')

        fragDictBasic = DB.fragDictBasic()
        GCMS = {}
        LCMS = {}
        CEMS = {}
        for abbrev in EMUlabel:
            if abbrev in fragDictBasic:
                mdv       = EMUlabel[abbrev]
                fragment  = fragDictBasic[abbrev]
                className = fragment.__class__.__name__

                if className == 'GCMSfragment':                  
                    match = regexp.match(abbrev)
                    fragmentName = 'M-'.join(match.groups())
                    GCMS[abbrev] = core.GCMSLabelData(dataTup=(fragmentName, mdv, []))                    
                if className == 'LCMSfragment':
                    fragmentName = fragment.abbrev
                    LCMS[abbrev] = labeling.LCMSLabelData(dataTup=(fragmentName, mdv, []))
                if className == 'CEMSfragment':
                    fragmentName = fragment.abbrev
                    CEMS[abbrev] = labeling.CEMSLabelData(dataTup=(fragmentName, mdv, []))
            else:
                raise Exception("Fragment "+abbrev+" not found in fragDictBasic")

        if GCMS:
            Notes['GCMSLabelData(fit)'] = GCMS
        if LCMS:
            Notes['LCMSLabelData(fit)'] = LCMS            
        if GCMS:
            Notes['CEMSLabelData(fit)'] = CEMS



class TwoSC13Model(C13Model):
    "Class for 2S-13CMFA calculations"
    
    def __init__(self,sbmlFileName):
        "sbmlFileName can be the name of the sbml file or a tuple of the form (fileName,string)"
        # GAMS file name
        GAMSfileName    = "GAMS_EMU_MFA_FLUX_HYB.gms"
        infeasFileName  = GAMSfileName.replace('gms','lst')     

        # sbml file
        reactionNetwork = ReactionNetworks.TSReactionNetwork(sbmlFileName)

        
        # Loading output functions
        outputFuncs = {'Info'           :(infoOut,['OFGenInfo.txt']),
                       'Vout'           :(fluxParOut,['Vout.txt','Vout']), 
                       'fout'           :(GAMSclasses.GAMSParameterFromResult,['fout.txt','fout']),
                       'labelcomp'      :(labelingParOut    ,['labelcomp.txt','labelcomp']),
                       'infeasibilities':(getInfeasibilities,[infeasFileName]),    
                       'VGEMout'        :(fluxParOut,['VGEMout.txt','VGEMout']), 
                       'VFBAgrowthout'  :(fluxParOut,['VFBAgrowthout.txt','VFBAgrowthout']), 
                       'Vgrowthout'     :(fluxParOut,['Vgrowthout.txt','Vgrowthout']), 
                        }
        
        self.getFilesFunctionName    = 'getTwoS13CMFAfiles'
        self.resultsFunctionNameRand = 'TSResults'
        self.resultsFunctionNameVar  = 'TSFVAResults'

        FluxModel.__init__(self,reactionNetwork,GAMSfileName,outputFuncs)


    def getGAMSInputFiles(self,reactionNetwork,labelMin=0.99,vSuggested='default',procString='raw'):
        """
        Get Input Files for the GAMS problem. 
        I separated this from createBatch for flexibility
        """
        
        GAMSInputFiles = reactionNetwork.getTwoS13CMFAfiles(vSuggested=vSuggested, labelCompNormMin=labelMin, procString=procString)       
        
        return GAMSInputFiles 


    def getSuggestedFluxes(self,Nrep,erase,procString):
        """Obtains suggested fluxes for createBatchVar"""

        # Find best solution through findFluxesStds and save as suggested fluxes
        resultsOne        = self.findFluxesStds(Nrep=Nrep,Nrand=0,erase=erase,limitFlux2Core=False,procString=procString)    
        
        # Adding previously calculated flux profile starting point for range procedure         
        self.reactionNetwork.C13ReacNet.reactionList.addFluxesWithDiffCheck(resultsOne.rangedFluxDict)
        fitFluxes    = GAMSclasses.fluxGAMSpar('vSUGG',[],fluxDict=self.reactionNetwork.C13ReacNet.reactionList.getFluxDictionary(level=2,rangeComp='best'))
        self.reactionNetwork.reactionList.addFluxesWithDiffCheck(resultsOne.rangedFluxDictAll)       
        fitFluxesFBA = GAMSclasses.fluxGAMSpar('vFBASUGG',[],fluxDict=self.reactionNetwork.reactionList.getFluxDictionary(level=1,rangeComp='best'))
        
        vSuggested = (fitFluxes,fitFluxesFBA)
        
        return vSuggested,resultsOne 


    def getGAMSFileRand(self,maxFlux13C,maxTime=[]):
        "Returns GAMS file for the ranged fluxes approach"
        gamsFileRand  = utils.file2tuple(self.dirGAMSFile+'GAMS_EMU_MFA_FLUX_HYB.gms')        
        gamsFileRand  = self.GAMSFileModifications(gamsFileRand,maxFlux13C,maxTime)      
        
        return gamsFileRand


    def getGAMSFileVar(self,maxFlux13C,maxTime=[]):
        "Returns GAMS file for the ranged fluxes approach"
        gamsFileVar  = utils.file2tuple(self.dirGAMSFile+'GAMS_EMU_MFA_FLUX_HYBvar.gms')        
        gamsFileVar  = self.GAMSFileModifications(gamsFileVar,maxFlux13C,maxTime)      
        
        return gamsFileVar

    
    def getOutputFuncsRand(self):
        """Returns output functions for randomization  approach"""
        
        # GAMS file name
        GAMSfileName    = self.getGAMSFileRand([])[0]
        infeasFileName  = GAMSfileName.replace('gms','lst')     

        # Start with the ones for 13C MFA
        outputFuncsRand = C13Model.getOutputFuncsRand(self)   
        
        # Add functions relevant to 2S only
        GAMSfileName    = "GAMS_EMU_MFA_FLUX_HYB.gms"          # TODO(Hector): clean this up by using getGAMSFileVar above
        infeasFileName  = GAMSfileName.replace('gms','lst')   

        outputFuncsRand['infeasibilities'] = (getInfeasibilities,[infeasFileName])
        outputFuncsRand['VGEMout']         = (fluxParOut        ,['VGEMout.txt','VGEMout'])
        outputFuncsRand['VFBAgrowthout']   = (fluxParOut        ,['VFBAgrowthout.txt','VFBAgrowthout'])
        outputFuncsRand['Vgrowthout']      = (fluxParOut        ,['Vgrowthout.txt','Vgrowthout'])

        return outputFuncsRand
            
    
    def getOutputFuncsVar(self):
        """Returns output functions for ranged fluxes approach"""
        
        # Start with the ones for 13C MFA
        outputFuncsVar = C13Model.getOutputFuncsVar(self)       
        
        # Add functions relevant to 2S only
        GAMSfileName    = "GAMS_EMU_MFA_FLUX_HYBvar.gms"          # TODO(Hector): clean this up by using getGAMSFileVar above
        infeasFileName  = GAMSfileName.replace('gms','lst')   

        outputFuncsVar['infeasibilities'] = (getInfeasibilities,[infeasFileName])
        outputFuncsVar['VGEMout']         = (fluxParOut        ,['Vout.txt','Vout'])
        outputFuncsVar['VFBAmaxout']      = (fluxParOut        ,['VFBAmaxout.txt','VFBAmaxout'])
        outputFuncsVar['VFBAminout']      = (fluxParOut        ,['VFBAminout.txt','VFBAminout'])
        outputFuncsVar.pop('Vmaxout')
        outputFuncsVar.pop('Vminout')

        return outputFuncsVar


    def getResultsStds(self,resultsDict):
        "Produces results out of the results dictionary"        
        results     = TSResults(resultsDict,self.reactionNetwork)           

        results.reactionNetwork.C13ReacNet.reactionList.addFluxesWithDiffCheck(results.rangedFluxDict)
        results.reactionNetwork.reactionList.addFluxesWithDiffCheck(results.rangedFluxDictAll)

        return results


    def getResultsRanges(self,resultsDict,resultsOne):
        "Produces final results out of the raw results dictionary"    
        # TODO(Hector): compact this with the C13 model version
        resultsVar        = TSFVAResults(resultsDict,self.reactionNetwork)           
        rangedFluxDictAll = resultsOne.rangedFluxDictAll
        rangedFluxDict    = resultsOne.rangedFluxDict


        # Use the resultsOne solution (best fit) as a basis
        results = copy.deepcopy(resultsVar)
        results.rangedFluxDict = resultsOne.rangedFluxDict
        results.rangedFluxDictAll = resultsOne.rangedFluxDictAll
        results.EMUlabel = resultsOne.EMUlabel
        results.OF = resultsOne.OF   

        # Use the results of the 13C FVA to obtain confidence intervals
        for name in rangedFluxDictAll:
            if name in resultsVar.rangedFluxDict:
                rangedFluxDictAll[name].net.hi   = resultsVar.rangedFluxDict[name][1] 
                rangedFluxDictAll[name].net.lo   = resultsVar.rangedFluxDict[name][0]
                
                rangedFluxDictAll[name].forward  = 'NA'
                rangedFluxDictAll[name].backward = 'NA'
                rangedFluxDictAll[name].exchange = 'NA'

        # These are not required, but we keep them around for debugging purposes.
        results.rangedFluxDictAll = rangedFluxDictAll
        results.rangedFluxDict    = rangedFluxDict
        
        # Storing result in reaction list
        results.reactionNetwork.C13ReacNet.reactionList.addFluxesWithDiffCheck(results.rangedFluxDict)
        results.reactionNetwork.reactionList.addFluxesWithDiffCheck(results.rangedFluxDictAll)             
        
        return results


    def testCoreRxns(self):
        "Test that there are core reactions included"
        C13ReactionList = self.reactionNetwork.reactionList.getCarbonTransitionReactions(self.reactionNetwork.getFluxRefScale())   
        assert len(C13ReactionList.reactions) > 0, "No core reactions included"        


    def findFluxesStds(self,Nrep=20,Nrand=5,erase=True,labelMin=0.99,fluxLimit=[0,0.05,0.2],limitFlux2Core=True,maxFlux13C=[],procString='raw'):
        "Finds fluxes and standard deviations"

        # Test that there are core reactions included
        self.testCoreRxns()
        
        if limitFlux2Core:
            self.limitFlux2Core(fluxLimit)
                
        return C13Model.findFluxesStds(self,Nrep,Nrand,erase,labelMin=labelMin,maxFlux13C=maxFlux13C,procString=procString)


    def findFluxesRanges(self,Nrep,fluxNames,erase=True,labelMin=0.99,fluxLimit=[0,0.05,0.2],limitFlux2Core=True,maxFlux13C=[],procString='raw'):
        "Finds fluxes and ranges of fluxes within allowable Objective Function"     

        # Test that there are core reactions included
        self.testCoreRxns()    
        
        if limitFlux2Core:
            self.limitFlux2Core(fluxLimit)                                   
           
        return C13Model.findFluxesRanges(self,Nrep,fluxNames,erase,labelMin=labelMin,maxFlux13C=maxFlux13C,procString=procString)
              
    
    def limitFlux2Core(self,limits=[0,0.05]):
        # TODO: simplify!
        "Limits flux into core metabolism to either 0 or limit input"
        FBAFile = self.reactionNetwork
        inFluxRef = self.reactionNetwork.inFluxRef
        Verbose = False
        warnFail= True if Verbose else False
        
        # Initial notes with (e.g) frag info on whether to include for fit
        initNotes    = self.reactionNetwork.notes
        C13initNotes = self.reactionNetwork.C13ReacNet.notes        
                
        # Raising error if biomass lower bound is zero, because it will create problems with SecureFluxToDest
        bmrxn = self.reactionNetwork.reactionList.getBiomassReaction()       
        if bmrxn.fluxBounds.net.lo <= 0:
            raise Exception('Lower bound of biomass must be bigger than zero.')
                
        # Get fluxes flowing into core metabolites
        reacsFlowingIntoCore = self.reactionNetwork.getReacsFlowingIntoCore()
       
        # Testing which fluxes can be set to zero        
        fluxBoundsDict = self.reactionNetwork.reactionList.getFluxBoundsDictionary()
        fbaModel = FBAModel(FBAFile.getSBMLString())

        if Verbose:        
            FBAFile.write('KEIOlimfileBefore.sbml')  
            print "limits:"
            print limits
        
        resultsFBA = fbaModel.findFluxes(warnFail=warnFail)
        if Verbose:
            print "successful?(limitFlux2Core before)"
            print resultsFBA.successful
        if not resultsFBA.successful:
            raise Exception('Extracellular fluxes infeasible')

        for reaction in reacsFlowingIntoCore.reactions:
            if Verbose:
                print "limiting flux to "+reaction.name
                    
            fluxBound = fluxBoundsDict[reaction.name]
            ub = fluxBound.net.hi
            lb = fluxBound.net.lo
            
            limitIter  = limits.__iter__()
            limit      = limitIter.next()
            while True:
                # Try setting to limit
                if reaction.hasForwardFlow:
                    fluxBound.net.hi = min(ub, abs(limit*inFluxRef)) 
                elif reaction.hasBackwardFlow:
                    fluxBound.net.lo = max(lb,-abs(limit*inFluxRef))                    
                else:
                    raise Exception("Reaction %s should not be included among those flowing into core" % reaction.name)
            
                fbaModel.changeFluxBoundsFAST(reaction.name,fluxBound)
                results = fbaModel.findFluxes(warnFail=warnFail)                         

                # Break if there is no next level or FBA was successful
                limitOld = limit
                try: 
                    limit = limitIter.next()
                except StopIteration:
                    break
                if results.successful:
                    break

            if Verbose:
                print "reaction %s set to level = %s" % (reaction.name, str(limitOld))

        resultsFBA = fbaModel.findFluxes(warnFail=warnFail)
        if not resultsFBA.successful:
            raise Exception('Flux limits to core too stringent!!')
        
        if Verbose:
            print "successful?(limitFlux2Core after)"
            print resultsFBA.successful
        
        # Storing changes in reactionNetwork
        if Verbose:
            fbaModel.reactionNetwork.write('KEIOlimfile.sbml')
        fbaModel.recordFASTchanges()
        self.reactionNetwork = ReactionNetworks.TSReactionNetwork(fbaModel.reactionNetwork.write('toString'))
        self.reactionNetwork.notes = initNotes
        self.reactionNetwork.C13ReacNet.notes = C13initNotes


    # External Labeling Variability Analysis (ELVA)
    def getELVAfile(self,results):
        "Function to obtain complementary reactions, sources and feeds for ELVA"

        def complement13C(Smat,Fluxes,inFluxRef):
            "Find complementary reactions, sources and feeds for ELVA" 
            carbonstr = 'abcdefghijklmnopqrst';
            mets = Smat.mets
            rxns = Smat.rxns
            S    = numpy.matrix(Smat.getSmatrix())
            v    = numpy.zeros((len(rxns),1))
            metList = self.reactionNetwork.C13ReacNet.metList
    
        
            for flux in Fluxes:
                v[rxns.index(flux)] = Fluxes[flux].net.best
              
            #   Finding metabolites involved in core reactions (and not excluded)
            metsCore = enhancedLists.MetaboliteList(set(self.reactionNetwork.C13ReacNet.metList.mets) - set(self.reactionNetwork.C13ReacNet.metList.getExcluded()))
            metsCoreDict = metsCore.metDict
 
        
            #   Finding reactions not involved in core 
            reacsCore   = set(self.reactionNetwork.C13ReacNet.reactionList.getReactionNameList(level=1)) 
            
            self.reactionNetwork.C13ReacNet.write('C13test.sbml')
            
            indReacsCore = []
            for reac in reacsCore:
                indReacsCore.append(rxns.index(reac))
            reacsNoCore = set(rxns) - set(self.reactionNetwork.C13ReacNet.reactionList.getReactionNameList(level=1))            
            indReacsNoCore = []
            for reac in reacsNoCore:
                indReacsNoCore.append(rxns.index(reac))
        
            #   Finding core mets with nonzero net flux coming in
            inout    = S[:,indReacsNoCore]*v[indReacsNoCore]
            coreFlux = S[:,indReacsCore]*v[indReacsCore]
            metsNonZeroTuple   = []
            metNonZeroFluxDict = {}
            for index in range(inout.shape[0]):
                if mets[index] in metsCoreDict and inout[index] != 0:
                    metsNonZeroTuple.append((mets[index],inout[index][0,0],coreFlux[index][0,0]))
                    
                    # Extra info for report. This could be written more elegantly.
                    fluxStDictCore    = {}
                    fluxStDictNonCore = {}                    
                    for Jind in range(S.shape[1]):
                        fluxSt = S[index,Jind]*v[Jind]
                        if fluxSt != 0:
                            if rxns[Jind] in reacsCore:
                                fluxStDictCore[rxns[Jind]] = fluxSt
                            else:
                                fluxStDictNonCore[rxns[Jind]] = fluxSt
                    if fluxStDictCore or fluxStDictNonCore:
                        metNonZeroFluxDict[mets[index]] = (fluxStDictCore,fluxStDictNonCore) 


            # Initialization for report output
            report = {}
        
            #   Final packing of newmets and newreactions
            newmets      = []
            newreactions = []
            carbondict = self.reactionNetwork.C13ReacNet.reactionList.getCarbonDict()
            metDict    = metList.metDict
            epsilon    = 0.00001

            for met,fluxinout,coreFlux in metsNonZeroTuple:
                ncarb = carbondict[met]            
                value = abs(fluxinout/inFluxRef)
                if fluxinout < 0:
                    name  = "RO" + met
                    line  = name + "\t" + met + " --> Out"+met+  "\t"  +carbonstr[0:ncarb]+" : "+carbonstr[0:ncarb]
                
                    react = core.Reactant(metDict[met],1)
                    prod  = core.Product(core.Metabolite('OUT'+met,ncarb),1)
                    flux  = core.flux(('NA','NA'),(core.rangedNumber(value-epsilon,value,value+epsilon),'NA')) 
                    newreactions.append(core.Reaction(name,reversible=False,suggested=False,measured=False,flux='None',fluxBounds=flux,transitionLine=line,exchange=False,reactants=[react],products=[prod]))
                    newmets.append(prod)
                    report[name] = value
                else:
                    name  = "RI" + met
                    line  = name + "\t" + "OUT"+ met +" --> "+met+  "\t"  +carbonstr[0:ncarb]+" : "+carbonstr[0:ncarb]
            
                    react = core.Reactant(core.Metabolite('OUT'+met,ncarb,source=True,feed='free'),1)
                    prod  = core.Product(metDict[met],1)
                    flux = core.flux(('NA','NA'),(core.rangedNumber(value-epsilon,value,value+epsilon),'NA')) 
                    newreactions.append(core.Reaction(name,reversible=False,suggested=False,measured=False,flux='None',fluxBounds=flux,transitionLine=line,exchange=False,reactants=[react],products=[prod]))
                    newmets.append(react)   
                    report[name] = value
        
            # Output to file
            outputFileName = 'ELVAfluxReport.txt'
            outputFile = open(outputFileName,'w')
            for name in sorted(report):
                stringOut = '{0:5.4f}'.format(report[name])
                outputFile.write(str(name)+"("+stringOut+"):\n")
                fluxStDictCore,fluxStDictNonCore = metNonZeroFluxDict[name.replace('RO','',1).replace('RI','',1)]

                outputFile.write("      ----> Non Core:\n")
                for fluxName in sorted(fluxStDictNonCore):
                    outputFile.write("      "+str(fluxName)+"   "+str(fluxStDictNonCore[fluxName])+"\n")  
    
                outputFile.write("      ----> Core:\n")
                for fluxName in sorted(fluxStDictCore):
                    outputFile.write("      "+str(fluxName)+"   "+str(fluxStDictCore[fluxName])+"\n")
  
            return newmets,newreactions


        # Find complementary reactions, sources and feeds for ELVA 
        inFluxRef = self.reactionNetwork.inFluxRef
        Smat = self.reactionNetwork.reactionList.getStoichMatrix()
        Fluxes  = results.rangedFluxDictAll 

        newmets,newreactions = complement13C(Smat, Fluxes, inFluxRef)

        # Old model
        oldC13Model = C13Model(self.reactionNetwork.C13ReacNet.getSBMLString())

        # Record fluxes and fits for fixing fluxes
        oldC13Model.recordFluxes(results.rangedFluxDict)
        oldC13Model.recordFits(results.EMUlabel)
        
        # Add new mets and reactions
        reactionList = oldC13Model.reactionNetwork.reactionList
        metList      = oldC13Model.reactionNetwork.metList
        
        newReactionsAll = copy.deepcopy(reactionList.reactions)
        newMetsAll      = copy.deepcopy(metList.mets)
        
        newReactionsAll.extend(newreactions)
        newMetsAll.extend(newmets)
        
        newReactionList = enhancedLists.ReactionList(newReactionsAll)
        newMetList      = enhancedLists.MetaboliteList(newMetsAll)
        
        # New model
        name  = 'ELVAmodel'
        Notes      = oldC13Model.reactionNetwork.notes
        # TODO: eliminate from Notes fragments not included in the fit?
                
        newReactNet = ReactionNetworks.C13ReactionNetwork( (name,Notes,newMetList,newReactionList))
        
        return newReactNet.write('toString')


    def ELVA(self,results,erase=True,labelMin=0.99,procString='raw'):
        "External Metabolite Variability Analysis as explained in Garcia Martin 2015, eqns 9-15"
        
        newSBMLFileName   = 'ELVAFileB.sbml'        
        newSBMLFileString = self.getELVAfile(results)
                
        ExpC13Model = ELVAModel(newSBMLFileString)

        ExpC13Model.reactionNetwork.write(newSBMLFileName)
        
        #newResults = ExpC13Model.findFluxesStds(dirFinal=os.getcwd()+'/',erase=erase,labelMin=labelMin)    
        newResults = ExpC13Model.findFluxesStds(erase=erase,labelMin=labelMin,procString=procString)

        return newResults


    def addFluxes(self,fluxDictionary):
        "Adds fluxes from dictionary to all reactions lists in model"
        fluxDict = copy.deepcopy(fluxDictionary)
        #    Genome-scale
        self.reactionNetwork.reactionList.addFluxes(fluxDict)
        inFluxRef = self.reactionNetwork.inFluxRef 
        #    core reactions
        for name in fluxDict.keys():
            fluxDict[name] = fluxDict[name]/inFluxRef             
        self.reactionNetwork.C13ReacNet.reactionList.addFluxes(fluxDict)
        self.reactionNetwork.C13ReacNet.reactionList.zeroFluxExchange()



class ELVAModel(FluxModel,C13Model):
    "Class for ELVA = External Labeling Variability Analysis calculations."
    # TODO(Hector): check the effect of fragments not to be fit

    def __init__(self,sbmlFileName):
        #TODO(Hector): recheck this class

        # GAMS file name
        GAMSfileName    = "GAMS_EMU_MFA_GCMSvar.gms"
        
        # Reaction network
        reactionNetwork = ReactionNetworks.C13ReactionNetwork(sbmlFileName)

        # Fix fluxes    (This should be done through method in sbml class)
        reactionList  = copy.deepcopy(reactionNetwork.reactionList)
        
        reactionList.fixFluxes(level=2)
        
        #reactionList.fixFluxes(level=1)        
        reactionNetwork = ReactionNetworks.C13ReactionNetwork( (reactionNetwork.name,reactionNetwork.notes,reactionNetwork.metList,reactionList))
        reactionNetwork.write('ELVAFile2.sbml')

        # Loading output functions
        outputFuncs = {'Info':(infoOut,['OFGenInfo.txt']),
                       'Vout':(fluxParOut,['Vout.txt','Vout']), 
                       'fout':(GAMSclasses.GAMSParameterFromResult,['fout.txt','fout']),
                       'labelcomp':(labelingParOut,['labelcomp.txt','labelcomp']),
                       #'infeasibilities':(getInfeasibilities,[infeasFileName])                           
                       'bounds':(GAMSclasses.GAMSParameterFromResult,['bounds.txt','bounds'])    
                        }
                        
        self.getFilesFunctionName    = 'get13CMFAfiles'
        self.resultsFunctionNameRand = 'ELVAResults'

        
        FluxModel.__init__(self,reactionNetwork,GAMSfileName,outputFuncs)


    def findFluxesStds(self,erase=True,labelMin=0.99,procString='raw'):
        "Find fluxes and standard deviations through monte carlo procedure"
        
        results = C13Model.findFluxesStds(self,Nrep=0,Nrand=0,erase=erase,labelMin=labelMin,procString=procString)

        return results


    def getGAMSFileRand(self,maxFlux13C,maxTime=[]):
        "Returns GAMS file for the ranged fluxes approach"
        gamsFileRand  = utils.file2tuple(self.dirGAMSFile+'GAMS_EMU_MFA_GCMSvar.gms')        
        gamsFileRand  = self.GAMSFileModifications(gamsFileRand,maxFlux13C,maxTime)      
        
        return gamsFileRand      


    def getOutputFuncsRand(self):
        """Returns output functions for randomization  approach"""

        # Loading output functions
        outputFuncs = {'Info':(infoOut,['OFGenInfo.txt']),
                       'Vout':(fluxParOut,['Vout.txt','Vout']), 
                       'fout':(GAMSclasses.GAMSParameterFromResult,['fout.txt','fout']),
                       'labelcomp':(labelingParOut,['labelcomp.txt','labelcomp']),
                       #'infeasibilities':(getInfeasibilities,[infeasFileName])                           
                       'bounds':(GAMSclasses.GAMSParameterFromResult,['bounds.txt','bounds'])    
                        }
                                
        return outputFuncs        


    def getResultsStds(self,resultsDict):
        "Produces results out of the results dictionary"        
        results     =   ELVAResults(resultsDict,self.reactionNetwork)           

        return results


    def createBatchRand(self,Nrep,Nrand,labelMin,maxFlux13C,erase,procString='raw'):
        "Creates Batch of gams problems"        

        reactionNetwork = self.reactionNetwork
        gamsFile = self.getGAMSFileRand(maxFlux13C)   
        outputFuncs = self.getOutputFuncsRand()                        
        fragDict     = reactionNetwork.fragDict()
        # Nrep and Nrad are ignored
        
        problems = []
        for abbrev in fragDict:
            frag = fragDict[abbrev]
            for m in range(len(frag.mdv)):
                Name    = '/Run'+abbrev+'M'+str(m)
                dirName = os.getcwd()+Name+'-'+GAMSclasses.GAMSProblem.getGAMSTempFolderName()
                GAMSInputFiles = self.getGAMSInputFiles(reactionNetwork,labelMin,procString=procString)
                
                # Add files for fragment and m to focus on
                fragFocus = GAMSclasses.GAMSSet('fragfocus',set([str(frag.abbrev)]))
                nFocus    = GAMSclasses.GAMSSet('nfocus',set(str(m)))
                
                fragFileName = 'fragfocus.txt'
                nFileName    = 'nfocus.txt'
                
                fragFocusSt = fragFocus.write('toString')
                nFocusSt    = nFocus.write('toString')
                
                GAMSInputFiles.extend([(fragFileName,fragFocusSt),(nFileName,nFocusSt)])
                # Writing problem
                GAMSprob       = GAMSclasses.GAMSProblem(Name,gamsFile,GAMSInputFiles,outputFuncs,directory = dirName,erase=erase)
                problems.append(GAMSprob)
    
        batch = GAMSclasses.GAMSProblemBatch(problems,erase=erase)
        
        return batch    



# Classes for results in GAMSProblem
class Results:
    "Base class for results"

    def __init__(self,resultsDict,reactionNetwork):
        self.reactionNetwork = copy.deepcopy(reactionNetwork)
        self.resultsDict = resultsDict
        self.fragDict = reactionNetwork.fragDict()
        self.processResults()
        delattr(self,'resultsDict')     # These files are big, better not to keep them
        
    
    def outputData(self,bestResults):  
        "Method"
        
        # flux fit
        FITFluxes   =  bestResults[0]['Vout'].getReactionList().getFluxDictionary()
        # Std calculation
        attribs = ['forward','backward','net','exchange']
        fluxesList = []
        for randomization in bestResults:
            if randomization != 0: 
                fluxes   =  bestResults[randomization]['Vout'].getReactionList().getFluxDictionary()
                fluxesList.append(fluxes)
        fluxNames = bestResults[0]['Vout'].getReactionList().getFluxDictionary().keys()
        fluxDictEnsemble = fluxDictEns(fluxesList,fluxNames,fluxDictRef=bestResults[0]['Vout'].getReactionList().getFluxDictionary())
        FITFluxesVar     = fluxDictEnsemble.getFluxDictStd(attribs)
        FITFluxesVarHi,FITFluxesVarLo  = fluxDictEnsemble.getFluxDictStdRef(attribs)

        # Ranged fluxes calculation
        ranged = {}
        rangedFluxDict  = {}
        for flux in sorted(FITFluxesVar[attribs[0]]):
            for attrib in attribs:
                ranged[attrib] = core.rangedNumber(FITFluxesVarLo[attrib][flux], getattr(FITFluxes[flux],attrib) , FITFluxesVarHi[attrib][flux])
        
            rangedFlux = core.flux((ranged['forward'],ranged['backward']),(ranged['net'],ranged['exchange']))
            rangedFluxDict[flux] = rangedFlux
        
          
        # Computational labeling
        EMULabelDict = bestResults[0]['labelcomp'].getLabelingDict(self.fragDict)
    
        # Objective function value
        OF           = bestResults[0]['Info'].OF
        
        return rangedFluxDict,EMULabelDict,OF    



    def findxyDxDy(self,fluxDictA,fluxDictB):
        "Finds values of x, y, DxU, DxL, DyU and DyL for plot"
        
        x    = []
        y    = []
        DxU  = []
        DxL  = []
        DyU  = []
        DyL  = []
        
        fluxComp = {}
        
        for flux in fluxDictA:
            if flux in fluxDictB:
                try:
                    xval  = fluxDictA[flux].net.best
                    xvalU = fluxDictA[flux].net.hi-fluxDictA[flux].net.best
                    xvalL = fluxDictA[flux].net.best-fluxDictA[flux].net.lo
                except AttributeError:
                    try:
                        xval  = fluxDictA[flux].net
                    except AttributeError:
                        xval  = fluxDictA[flux]
                    xvalU = 0
                    xvalL = 0
                    
                try:
                    yval  = fluxDictB[flux].net.best
                    yvalU = fluxDictB[flux].net.hi-fluxDictB[flux].net.best
                    yvalL = fluxDictB[flux].net.best-fluxDictB[flux].net.lo
                except AttributeError:
                    try:
                        yval  = fluxDictB[flux].net
                    except AttributeError:
                        yval  = fluxDictB[flux]
                    yvalU = 0
                    yvalL = 0
                
                
                x.append(xval)
                DxU.append(xvalU)
                DxL.append(xvalL)
                 
                y.append(yval)
                DyU.append(yvalU)
                DyL.append(yvalL)
                
                fluxComp[flux] = [xval,xval-xvalL,xval+xvalU,yval,yval-yvalL,yval+yvalU]
        
        x    = numpy.array(x)    
        dxU  = numpy.array(DxU)  
        dxL  = numpy.array(DxL)  
              
        y    = numpy.array(y)    
        dyU  = numpy.array(DyU)  
        dyL  = numpy.array(DyL)
        
        return x,dxU,dxL,y,dyU,dyL,fluxComp


    def createFigure(self,titleFig='',axes='default',Xlabel='External fluxes',Ylabel='Internal fluxes',xMax=2.2,yMax=2.2):
        "Creates figure"
        if axes == 'default':
            axes = figure(figsize=(8,8)).add_subplot(111)
        axes.plot([0,xMax],[0,yMax],'k')
        axes.axis([0,xMax,0,yMax])

        xlabel(Xlabel)
        ylabel(Ylabel)

        title(titleFig)     
    
        return axes
        

    def compareFluxes(self, FluxDict, titleFig='',outputFileName="fluxComparison.txt",axes='default',Xlabel='default',Ylabel='default',Xmax='default',Ymax='default'):        
        
        self.compareFluxesBase(self.rangedFluxDict, FluxDict, titleFig=titleFig, outputFileName=outputFileName, axes=axes, Xlabel=Xlabel, Ylabel=Ylabel, Xmax=Xmax, Ymax=Ymax)


    def compareFluxesBase(self,fluxDictA,fluxDictB,titleFig='',outputFileName="fluxComparison.txt",axes='default',Xlabel='default',Ylabel='default',Xmax='default',Ymax='default'):
        "Core function"
        
        # Obtaining x and y data
        x,dxU,dxL,y,dyU,dyL,fluxComp = self.findxyDxDy(fluxDictA,fluxDictB)

        # Plotting data
        Xmax = 1.05*max(x) if Xmax=='default' else Xmax
        Ymax = 1.05*max(x) if Ymax=='default' else Ymax
        axes = self.createFigure(titleFig=titleFig,axes=axes,Xlabel=Xlabel,Ylabel=Ylabel,xMax=Xmax,yMax=Ymax)
        axes.errorbar(x,y,yerr=[dyL,dyU],xerr=[dxL,dxU],fmt='o')

        # Output to file
        outputFile = open(outputFileName,'w')
        outputFile.write('Local fluxes (best,lo,hi), Outside fluxes (best,lo,hi):\n')
        for flux in fluxComp:
            outputFile.write(str(flux)+": "+str(fluxComp[flux])+"\n")



    def drawFluxesBase(self, reactionNetwork, normFluxDict, fileName, svgInFileName='Ecolireactionslevel_1J.svg', title='', oldVersion=False):        
        
       if not oldVersion: 
           # Only add the default path if the given path is not absolute
           if svgInFileName[:1] != '/':
               qmodeldir   = os.environ['QUANTMODELPATH']
               svgDir      = qmodeldir+'/code/core/'
               svgin = svgDir+svgInFileName
           else:
                svgin = svgInFileName

           print "svgin:"
           print svgin     
        
           # instantiate FluxMap object
           fmap = FluxMap(inputFilename=svgin)

           # change all fluxes
           fmap.decorateMap(reactionNetwork, normFluxDict)

           # write the svg out
           fmap.writeSVG(fileName)        
       
       else:   # TODO(Hector): deprecate this feature as soon as all maps are actualized
           import FluxMaps
           
           qmodeldir   = os.environ['QUANTMODELPATH']
           svgDir      = qmodeldir+'/code/core/' 

           # set parameters for flux map
           svgin=svgDir+svgInFileName
           print "svgin:"
           print svgin     
        
           # instantiate FluxMap object
           fmap = FluxMaps.FluxMapOLD(svgin)
        
           # change title
           fmap.changetitle(title)

           # change all fluxes
           fmap.changeallfluxes(normFluxDict)

           # write the svg out
           fmap.writesvg(fileName)    


    def drawFluxes(self,fileName,svgInFileName='Ecolireactionslevel_1J.svg',norm='',titleTag=' (13C FVA)',oldVersion=False):        
    
        # Normalizing Flux dictionary
        norm = self.reactionNetwork.getHighestCarbonInput() if not norm else norm  
        
        normFluxDict = copy.deepcopy(self.reactionNetwork.reactionList.getFluxDictionary(level=1,fluxComp='net',rangeComp='all',norm=norm))
  
        # Passing it to base function
        title = self.reactionNetwork.name+self.titleTag
        self.drawFluxesBase(self.reactionNetwork, normFluxDict, fileName, svgInFileName=svgInFileName, title=title, oldVersion=oldVersion)         


    def plotMetaboliteFlux(self,metName, minFluxGroup = 0.05, maxFluxGroup=1.5,titleFig='',save='default'):
        "Sankey plot for fluxes producing and consuming a metabolite. minFluxGroup is the minimum amount of flux so it does not get grouped with others."
        # TODO: make more succint
        # Fluxes for genome-scale model
        try:
            fluxes  = self.rangedFluxDictAll
        except:
            fluxes  = self.rangedFluxDict            
                
        
        # Stoichiometry matrix
        Smat = self.reactionNetwork.reactionList.getStoichMatrix()
        mets    = Smat.mets
        rxns    = Smat.rxns
        S    = Smat.getSmatrix()
        
        # flows and labels
        flows   = []
        labels  = []
        #tempflux = 0
        #temprxns = ''
        units = []
        imet = mets.index(metName)             
        newFluxes = {}
        for rxn in rxns:
             jrxn = rxns.index(rxn)   
                     
             if S[imet][jrxn] != 0 and fluxes[rxn].net.best != 0:      
                flows.append(round(fluxes[rxn].net.best*S[imet][jrxn],2))
                stringOut = rxn       
                labels.append(stringOut)
                newFluxes[rxn] = S[imet][jrxn]*fluxes[rxn]
                units.append((newFluxes[rxn].net,rxn))       
        
        # Group flow and labels
        groups = []
        groupPlus  = (core.rangedNumber(0,0,0),'')
        groupMinus = (core.rangedNumber(0,0,0),'')
        for unit in units:
            flux,label = unit
            if abs(flux.best) >= minFluxGroup:
                groups.append(unit)
            else:
                if flux.best >= 0:
                    gPflux,gPlabel = groupPlus
                    gPflux  = gPflux + flux
                    gPlabel = gPlabel+'\n'+label if gPlabel != '' else label
                    groupPlus = (gPflux,gPlabel)
                    if abs(gPflux.best) >= maxFluxGroup:
                        groups.append((gPflux,gPlabel))
                        groupPlus  = (core.rangedNumber(0,0,0),'')
                else:
                    gMflux,gMlabel = groupMinus
                    gMflux  = gMflux + flux
                    gMlabel = gMlabel+'\n'+label if gMlabel != '' else label
                    groupMinus = (gMflux,gMlabel)
                    if abs(gMflux.best) >= maxFluxGroup:
                        groups.append((gMflux,gMlabel))
                        groupMinus  = (core.rangedNumber(0,0,0),'')
        # Take care of last group
        if abs(groupPlus[0].best) < maxFluxGroup and groupPlus[1] != '':
            groups.append((gPflux,gPlabel))
        if abs(groupMinus[0].best) < maxFluxGroup and groupMinus[1] != '':
            groups.append((gMflux,gMlabel))
           
        newFlows  = []
        newLabels = []
        for group in groups:
            flux,label = group
            newFlows.append(round(flux.best,2))
            if flux.best >= 0:
                newLabels.append(label+'\n'+'['+str(round(abs(flux.lo),2))+" - "+str(round(abs(flux.hi),2))+']')
            else:
                newLabels.append(label+'\n'+'['+str(round(abs(flux.hi),2))+" - "+str(round(abs(flux.lo),2))+']')

        # pathlengths 
        pathlengths = []        
        for i in range(len(newLabels)):
             pathlengths.append(0.25)

        # orientations
        orientations = []
        
        for i in range(len(newLabels)):
            if i < (len(newLabels))/2:
                 orientations.append(-1)
            else:
                  orientations.append(1)
                  
        maxvalue = newFlows.index(max(newFlows))
        minvalue = newFlows.index(min(newFlows))
        orientations[maxvalue] = 0
        orientations[minvalue] = 0
        patchlabel = metName[0:-2].upper()

        # Sankey plot
        from matplotlib.sankey import Sankey
        matplotlib.rcParams.update({'font.size': 12})

        #fig = figure()    
        fig = figure(figsize=(18,18))
        title = patchlabel+" balances "+titleFig
        ax = fig.add_subplot(1, 1, 1, xticks=[], yticks=[],title=title)   
        sankey = Sankey(ax=ax, scale=0.15, offset=0.12)
        sankey.add(       
                   labels=newLabels,
                   flows=newFlows,
                   orientations = orientations,
                   pathlengths = pathlengths,
                   patchlabel=patchlabel,
                   alpha=0.2, lw=2.0)    # Arguments to matplotlib.patches.PathPatch()
        
        diagrams = sankey.finish()
        diagrams[0].patch.set_facecolor('b')
        diagrams[0].text.set_fontweight('bold')
        diagrams[0].text.set_fontsize(25)
        
        matplotlib.rcParams.update({'font.size': 20})

        # Save fig        
        if save != 'default':
            savefig(save)



class FBAResults(Results):
    def __init__(self,resultsDict,reactionNetwork):
        self.reactionNetwork = copy.deepcopy(reactionNetwork)
        self.resultsDict = resultsDict
        self.rangedFluxDict = resultsDict   # This should be compacted with the previous line
        self.processResults()
        self.titleTag = ' (FBA)'  
        delattr(self,'resultsDict')     # This files are big, better not to keep them
        delattr(self,'rangedFluxDict')

        
    def processResults(self):
        Vout             = self.resultsDict['Vout']
        self.VFBA        = GAMSclasses.fluxGAMSpar(Vout.name,Vout.elements)
        self.successful  = self.resultsDict['Successful']
        
        # Ranged flux dict 
        # TODO(Hector): this should be part of its own outputData() method
        fluxDict = self.VFBA.getReactionList().getFluxDictionary()
        rangedFluxDict = {}
        for name in fluxDict:
            flux = fluxDict[name]
            rangedForw = core.rangedNumber(flux.forward ,flux.forward ,flux.forward) 
            rangedBack = core.rangedNumber(flux.backward,flux.backward,flux.backward) 
            rangedNet  = core.rangedNumber(flux.net     ,flux.net     ,flux.net) 
            rangedExch = core.rangedNumber(flux.exchange,flux.exchange,flux.exchange) 
            
            rangedFlux = core.flux((rangedForw,rangedBack),(rangedNet,rangedExch))
        
            rangedFluxDict[name] = rangedFlux

        self.rangedFluxDict = rangedFluxDict
        self.reactionNetwork.reactionList.addFluxes(fluxDict)


    def printSuccess(self,onlyFail=True):
        string =  ' successful'
        if not self.successful:
            string = ' NOT'+string
        string = 'lp'+string
        
        # Print only if asked
        if onlyFail:
            if not self.successful:
                print string
        else:
            print string



class FVAResults(FBAResults):
    
    def __init__(self,resultsDict,reactionNetwork):
        FBAResults.__init__(self,resultsDict,reactionNetwork)
        self.titleTag = ' (FVA)'   


    def compareFluxes(self, FluxDict, titleFig='',outputFileName="fluxComparison.txt",axes='default',Xlabel='default',Ylabel='default',Xmax='default',Ymax='default'):        
        
        self.compareFluxesBase(self.resultsDict,FluxDict,titleFig=titleFig,outputFileName=outputFileName,axes=axes,Xlabel=Xlabel,Ylabel=Ylabel,Xmax=Xmax,Ymax=Ymax)


    def processResultsOUT(self):
        self.reactionNetwork.reactionList.addFluxes(self.rangedFluxDict)


    def processResults(self):
        
        resultsDict = self.resultsDict

        # Output parsing
        VFBAOut    = resultsDict['Vout']
        VFBAOutMax = resultsDict['Vmax']
        VFBAOutMin = resultsDict['Vmin']
        
        VFBA    = GAMSclasses.fluxGAMSpar(VFBAOut.name,VFBAOut.elements).getReactionList().getFluxDictionary()
        VFBAmax = GAMSclasses.fluxGAMSpar(VFBAOutMax.name,VFBAOutMax.elements).getReactionList().getFluxDictionary()
        VFBAmin = GAMSclasses.fluxGAMSpar(VFBAOutMin.name,VFBAOutMin.elements).getReactionList().getFluxDictionary()
        
        # Ranged fluxes calculation
        attribs = ['forward','backward','net','exchange']
        ranged = {}
        rangedFluxDict  = {}
        for flux in VFBAmax:
            for attrib in attribs:
                ranged[attrib] = core.rangedNumber(getattr(VFBAmin[flux],attrib),getattr(VFBA[flux],attrib),getattr(VFBAmax[flux],attrib))
            rangedFluxDict[flux] = core.flux((ranged['forward'],ranged['backward']),(ranged['net'],ranged['exchange']))
            
        self.reactionNetwork.reactionList.addFluxes(rangedFluxDict)
        
        self.rangedFluxDict = rangedFluxDict



class C13Results(Results):      
    
    def __init__(self,resultsDict,reactionNetwork):
        Results.__init__(self,resultsDict,reactionNetwork)
        self.titleTag = ' (13C)'    


    def processResults(self):
        # Find lowest OFs for each randomization
        bestResults,bestOFs = self.findLowestOF()
        
        # Output data
        rangedFluxDict,EMUlabel,OF  = self.outputData(bestResults)
                
        self.rangedFluxDict = rangedFluxDict
        self.EMUlabel       = EMUlabel
        self.OF             = OF


    def findLowestOF(self):
        # Find lowest objective functions which are feasible
        resultsDict = self.resultsDict
        bestResults = {}
        bestOFs     = {}
        infeasAll   = {}
        for name in resultsDict:
            result = resultsDict[name]
            randomization  = float(name.lstrip('/Run').split('_')[0])
            OF             = result['Info'].OF
            solvestat      = result['Info'].solvestat
            modelstat      = result['Info'].modelstat
            infeasAll[randomization] = result['infeasibilities']
            if solvestat == 1 and (modelstat == 1 or modelstat== 2):
                if randomization not in bestOFs:
                    bestOFs[randomization]     = OF
                    bestResults[randomization] = result
                else:
                    if OF < bestOFs[randomization]:
                        bestOFs[randomization] = OF
                        bestResults[randomization] = result
        
        # Should find a better way to handle this exception
        if not bestOFs:  # No feasible solutions
            print "NO FEASIBLE SOLUTIONS!!!"
            for line in infeasAll[randomization]:
                print line

            raise Exception('NO FEASIBLE SOLUTIONS!!!')            

        return bestResults,bestOFs
        
    def plotMetaboliteFlux(self,metName, minFluxGroup = 0.05, maxFluxGroup=1.5,titleFig='',save='default'):        
        print "Method not available"


    def plotExpvsCompLabelXvsY(self,titleFig='',axes='default',color='b',fmt='.',save='default'):
        # Obtaining experimental data
        try:
            fragDict = self.fragDict
        except AttributeError:
            raise Exception('GCMS Dictionary not available')
    
        # Obtaining fit data
        EMULabelDict = self.EMUlabel
        
        # Obtaining intersection of keys
        allKeys = list(set(fragDict.keys()).intersection(set(EMULabelDict.keys())))        
        
        # Obtaining fragments not included in the fit
        nonFitKeys = []
        fitKeys    = []
        for key in allKeys:
            if fragDict[key].inFit:
                fitKeys.append(key)
            else:
                nonFitKeys.append(key)
        
        # Creating axes if needed
        if axes =='default':
            axes = figure(figsize=(8,8)).add_subplot(111)        
        axes.plot([0,1],[0,1],'k')
        axes.axis([0,1,0,1])
        
        xlabel('Experimental GCMS')
        ylabel('Computational GCMS')
    
        utils.XvsYcomp(fragDict,EMULabelDict,fitKeys,axes,color,fmt)
        utils.XvsYcomp(fragDict,EMULabelDict,nonFitKeys,axes,'g',fmt)

        
        title(titleFig,fontweight='bold')
        # Save fig        
        if save != 'default':
            savefig(save)  

    def plotExpvsCompLabelFragment(self,fragment = 'all',outputFileName='ExpvsCompLabelFragment.txt',titleFig='',save='default'):
        # Obtaining experimental data
        try:
            fragDict = self.fragDict
        except AttributeError:
            raise Exception('GCMS Dictionary not available')
    
        # Obtaining fit data
        EMULabelDict = self.EMUlabel

        # Adapting input
        if fragment == 'all':
            fragments = EMULabelDict.keys()
        elif isinstance(fragment,str):
            fragments = [fragment]

        # Doing plot            
        utils.histcomp(EMULabelDict,fragDict,fragments,outputFileName,titleFig,self.OF,save)
 




class TSResults(C13Results):
    def __init__(self,resultsDict,reactionNetwork):
        self.reactionNetwork = copy.deepcopy(reactionNetwork)
        self.fragDict = self.reactionNetwork.C13ReacNet.fragDict()
        self.resultsDict = resultsDict   
        self.processResults()
        self.titleTag = ' (2S)'    
        delattr(self,'resultsDict')     # This files are big, better not to keep them

    
    def processResults(self):
        bestResults,bestOFs = self.findLowestOF()
        rangedFluxDictAll,rangedFluxDictCore,EMUlabel,OF,fluxDictEnsemble  = self.outputData(bestResults)
        
        # Add fluxes to reaction list        
        
        self.rangedFluxDict     = rangedFluxDictCore
        self.rangedFluxDictAll  = rangedFluxDictAll
        self.EMUlabel           = EMUlabel
        self.OF                 = OF
        self.fluxDictEnsemble   = fluxDictEnsemble
        
            
    def outputData(self,bestResults):
        # 13C part
        rangedFluxDict13C,EMUlabel,OF = C13Results.outputData(self,bestResults)
        # Genome-scale part (repeated code, it should be improved)
        FITFluxes   =  bestResults[0]['VGEMout'].getReactionList().getFluxDictionary()
        #     Std calculation
        attribs = ['forward','backward','net','exchange']
        fluxesList = []
        for randomization in bestResults:          
            if randomization != 0: 
                fluxes   =  bestResults[randomization]['VGEMout'].getReactionList().getFluxDictionary()
                fluxesList.append(fluxes)
        fluxNames = bestResults[0]['VGEMout'].getReactionList().getFluxDictionary().keys()
        fluxDictEnsemble = fluxDictEns(fluxesList,fluxNames,fluxDictRef=bestResults[0]['VGEMout'].getReactionList().getFluxDictionary())
        FITFluxesVar     = fluxDictEnsemble.getFluxDictStd(attribs)
        FITFluxesVarHi,FITFluxesVarLo  = fluxDictEnsemble.getFluxDictStdRef(attribs)

        #     Ranged fluxes calculation
        ranged = {}
        rangedFluxDict  = {}
        for flux in sorted(FITFluxesVar[attribs[0]]):
            for attrib in attribs:
                ranged[attrib] = core.rangedNumber(FITFluxesVarLo[attrib][flux], getattr(FITFluxes[flux],attrib) , FITFluxesVarHi[attrib][flux])
        
            rangedFlux = core.flux((ranged['forward'],ranged['backward']),(ranged['net'],ranged['exchange']))
            rangedFluxDict[flux] = rangedFlux
        
        return rangedFluxDict,rangedFluxDict13C,EMUlabel,OF,fluxDictEnsemble


    def compareFluxesCore(self,FluxDict,titleFig='',outputFileName="fluxComparisonCore.txt",axes='default'):
        # Using function from C13
        C13Results.compareFluxes(self,FluxDict, titleFig = titleFig, outputFileName = outputFileName, axes = axes)


    def compareFluxesAll(self, FluxDict, titleFig='',outputFileName="fluxComparisonAll.txt",axes='default',Xlabel='default',Ylabel='default',Xmax='default',Ymax='default'):        
        
        self.compareFluxesBase(self.rangedFluxDictAll,FluxDict,titleFig=titleFig,outputFileName=outputFileName,axes=axes,Xlabel=Xlabel,Ylabel=Ylabel,Xmax=Xmax,Ymax=Ymax)


    def printFluxesAll(self,fileName='None',names='default'):
        return self.printFluxes(fileName,names,fluxDict = self.rangedFluxDictAll)



# Classes for results in ELVA       
class ELVAResults(Results):
        
    def processResults(self):
        labelMinMax = self.outputData()
        
        self.labelMinMax = labelMinMax


    def outputData(self):
        labelMinMax = {}
        resultsDict = self.resultsDict
        regex       = re.compile('([\w\-]+)(M\d+)')
        
        for dirName in resultsDict:
            result = resultsDict[dirName]
            matches = regex.match(dirName.lstrip('/Run'))
            abbrev = matches.group(1)
            mVal   = matches.group(2).lstrip('M')
            if not abbrev in labelMinMax:
                labelMinMax[abbrev]= {}
            labelMinMax[abbrev][mVal]=[result['bounds'].elements[('min',)],result['bounds'].elements[('max',)]]
        
        return labelMinMax


    def getxDXyDyInfo(self):
        "Produces plotting info for plotExpvsCompLabelxvsy method"
        
        # Obtaining experimental data
        try:
            fragDict = self.fragDict
        except AttributeError:
            raise Exception('GCMS Dictionary not available')
    
        # Obtaining fit data
        labelMinMax = self.labelMinMax
        
        # Obtaining plot data
        x  = []  # Computational data (labelMinMax)
        Dx = []
        y  = []  # Experimental data (fragDict)  
        DyU= []
        DyL= []
    
        # File output lines initialization
        labelComp     = {}
        labelCompDict = {}
    
        # Getting data for comparison
        for abbrev in labelMinMax:
            if abbrev in fragDict:
                frag = fragDict[abbrev]
                mdv    = frag.mdv
                MDVfit = frag.MDVfit
                std  = frag.std
                                
                for m in range(min(len(mdv),len(MDVfit))):
                    x.append(mdv[m])
                    Dx.append(std[m]/2)
                    y.append(MDVfit[m])
                    DyU.append( (labelMinMax[abbrev][str(m)][1]-MDVfit[m]))
                    DyL.append(-(labelMinMax[abbrev][str(m)][0]-MDVfit[m]))
                    #print "y,Dys: "+str(MDVfit[m])+" ; "+str(labelMinMax[abbrev][str(m)][1]-MDVfit[m])+" , "+str(-(labelMinMax[abbrev][str(m)][0]-MDVfit[m]))
                    labelComp[abbrev+'M'+str(m)] = str(MDVfit[m])+' '+str(labelMinMax[abbrev][str(m)][0])+'-'+str(labelMinMax[abbrev][str(m)][1])+'---->'+str(labelMinMax[abbrev][str(m)][1]-labelMinMax[abbrev][str(m)][0])
                    labelCompDict[abbrev+'M'+str(m)] = ( MDVfit[m], labelMinMax[abbrev][str(m)][0], labelMinMax[abbrev][str(m)][1] )
    
        # Plotting
        x = numpy.array(x)
        Dx= numpy.array(Dx)
        y = numpy.array(y)
        DyU= numpy.array(DyU)
        DyL= numpy.array(DyL)
        
        return x,Dx,y,DyU,DyL,labelComp,labelCompDict


    def plotExpvsCompLabelxvsy(self,titleFig='',axes='default',outputFileName="ELVAComparison.txt",save='default'):
        """
        Makes ELVA plot showing computational error
        """
        x,Dx,y,DyU,DyL,labelComp,labelCompDict = self.getxDXyDyInfo()
        
        if axes =='default':
            axes = figure(figsize=(8,8)).add_subplot(111)
        if (len(x) != 0) and (len(y) != 0):
            axes.errorbar(x,y,yerr=[DyL,DyU],xerr=[Dx,Dx],fmt='ko')    
            axes.plot([0,1],[0,1],'k')
            axes.axis([0,1,0,1])
            
            xlabel('Experimental Labeling')
            ylabel('Computational Labeling')
    
            title(titleFig,fontweight='bold')
        else:
            print "No data in plot"
            print x
            print y

        # Output to file
        outputFile = open(outputFileName,'w')
        outputFile.write('Y axis of ELVA plot:value, low value and high value\n')
        for name in sorted(labelComp):
            outputFile.write(str(name)+": "+str(labelComp[name])+"\n")

        # Save fig        
        if save != 'default':
            savefig(save)  


    def printLabeling(self,fileName='None',names='default'):
        labelMinMax = self.labelMinMax
        
        lines = []
        
        # Metabolites to print
        if names == 'default':
            names = labelMinMax.keys()
        
        
        for fragment in names:
            lines.append(fragment+':\n')
            lines.append(labelMinMax[fragment].__str__())
            lines.append('\n')
        
        # Output
        if fileName == 'None':    # Print to screen
            for line in lines:
               print line
        elif fileName == 'toString':
            return ''.join(lines)
        else:                   # Print to file
            file = open(fileName,'w')
            for line in lines:
                file.write(line)
            file.close()        



class C13FVAResults(Results):

    def __init__(self,resultsDict,reactionNetwork):
        self.reactionNetwork = copy.deepcopy(reactionNetwork)
        self.resultsDict = resultsDict
        self.fragDict = reactionNetwork.fragDict()
        self.processResults()
        self.titleTag = ' (13C FVA)'    
        delattr(self,'resultsDict')     # This files are big, better not to keep them

    
    def processResults(self):
        # Barebones, we probably need to output OF and EMUlabel (from SUGGESTED SOLUTION) as well.
        rangedFluxDict = self.outputData(self.resultsDict)
    
        self.rangedFluxDict = rangedFluxDict    


    def outputData(self,resultsDict):
                
        # Ranged fluxes dictionary
        rangedFluxDict = {}
        for name in resultsDict:
            result = resultsDict[name]
            topVal    = result['fluxBounds'].elements[('max',)]
            bottomVal = result['fluxBounds'].elements[('min',)]
            rangedFluxDict[name.replace('/Run_','',1)] = (bottomVal,topVal)
            
        return rangedFluxDict    



class TSFVAResults(C13FVAResults,TSResults):

    def __init__(self,resultsDict,reactionNetwork):
        TSResults.__init__(self,resultsDict,reactionNetwork)
        self.titleTag = ' (TS FVA)'
    
    def outputData(self,resultsDict):
        
        # Ranged fluxes dictionary
        rangedFluxDict = {}
        for name in resultsDict:
            result = resultsDict[name]
            topVal    = result['fluxBounds'].elements[('max',)]
            bottomVal = result['fluxBounds'].elements[('min',)]
            
            fluxName    = name.replace('/Run_','',1)
            fluxNameNew = fluxName            
            
            rangedFluxDict[fluxNameNew] = (bottomVal,topVal)
        
        # Debug
        debug = False
        if debug:      
            print "in debug"
            fluxName = 'G6PDH2r'
            dir ='/scratch/hgmartin_gams_files/tests/batch/bigtestp02/'
            
            # max
            result = resultsDict['/Run_'+fluxName] 
            rangedFluxDictAll = result['VFBAmaxout'].getReactionList().getFluxDictionary(level=1,fluxComp='net',rangeComp='all',norm=self.reactionNetwork.getHighestCarbonInput())
            self.drawFluxesBase(self.reactionNetwork, rangedFluxDictAll, dir+fluxName+'Maxtest.svg', svgInFileName='KEIO3.svg')

            self.fragDict = self.reactionNetwork.C13ReacNet.fragDict()
            self.EMUlabel = result['labelcompMax'].getLabelingDict(self.fragDict)
            self.OF       = result['InfoMax'].OF

            self.plotExpvsCompLabelFragment(titleFig='Max')

            # min
            result = resultsDict['/Run_'+fluxName]
            rangedFluxDictAll = result['VFBAminout'].getReactionList().getFluxDictionary(level=1,fluxComp='net',rangeComp='all',norm=self.reactionNetwork.getHighestCarbonInput())
            self.drawFluxesBase(self.reactionNetwork, rangedFluxDictAll, dir+fluxName+'Mintest.svg', svgInFileName='KEIO3.svg')

            self.fragDict = self.reactionNetwork.C13ReacNet.fragDict()
            self.EMUlabel = result['labelcompMin'].getLabelingDict(self.fragDict)
            self.OF       = result['InfoMin'].OF

            self.plotExpvsCompLabelFragment(titleFig='Min')

        return rangedFluxDict





class fluxDictEns():
    " Class for ensemble of flux dictionaries"
    
    def __init__(self,fluxDictList,fluxNames,fluxDictRef=[]):
        self.fluxDictList = fluxDictList
        self.fluxNames    = fluxNames
        self.fluxDictRef  = fluxDictRef


    # Gets dictionary of ranged fluxes
    # attribs is the attributes to average over (e.g. 'forward','backward', 'net', 'exchange')
    def getFluxDictStd(self,attribs):
        # Std calculation
        fluxDictList = self.fluxDictList 
        
        # Initialization
        FluxesDictVar = {}
        for attrib in attribs:
            FluxesDictVar[attrib] = {}
            for flux in self.fluxNames:
                FluxesDictVar[attrib][flux] = []
        if fluxDictList:        
            # Getting all attribute values into a single list
            for fluxDict in fluxDictList:
                for flux in fluxDict:
                    for attrib in attribs:
                        FluxesDictVar[attrib][flux].append(getattr(fluxDict[flux],attrib))
            # Finding std for whole list
            for attrib in attribs:
                for flux in FluxesDictVar[attrib]:
                    FluxesDictVar[attrib][flux] = numpy.std(FluxesDictVar[attrib][flux])        
        else:
            for attrib in attribs:
                for flux in FluxesDictVar[attrib]:
                    FluxesDictVar[attrib][flux] = 0  

        return FluxesDictVar


    def getFluxDictStdRef(self,attribs):
        "Getting mean deviations from above and below reference"
        fluxDictList = self.fluxDictList 
        fluxDictRef  = self.fluxDictRef

        # Initialization
        FluxesDictVarHi = {}
        FluxesDictVarLo = {}
        FluxesDictRef   = {}
        for attrib in attribs:
            FluxesDictVarHi[attrib]  = {}
            FluxesDictVarLo[attrib]  = {}
            FluxesDictRef[attrib] = {}
            for flux in self.fluxNames:
                FluxesDictVarHi[attrib][flux] = []
                FluxesDictVarLo[attrib][flux] = []
                FluxesDictRef[attrib][flux] = getattr(fluxDictRef[flux],attrib)
        if fluxDictRef:
            if fluxDictList:    
                # Getting all attribute values into a single list
                for fluxDict in fluxDictList:
                    for flux in fluxDict:
                        for attrib in attribs:
                            value    = getattr(fluxDict[flux],attrib)
                            valueRef = getattr(fluxDictRef[flux],attrib)
                            if value > valueRef:
                                FluxesDictVarHi[attrib][flux].append(value)                                
                            elif value < valueRef:
                                FluxesDictVarLo[attrib][flux].append(value)    
                    
                # Finding std for whole list
                for attrib in attribs:
                    for flux in FluxesDictVarHi[attrib]:
                        if not FluxesDictVarHi[attrib][flux]:
                            FluxesDictVarHi[attrib][flux] = [getattr(fluxDictRef[flux],attrib)]
                        if not FluxesDictVarLo[attrib][flux]:
                            FluxesDictVarLo[attrib][flux] = [getattr(fluxDictRef[flux],attrib)]
                        FluxesDictVarHi[attrib][flux] = numpy.mean(FluxesDictVarHi[attrib][flux])   
                        FluxesDictVarLo[attrib][flux] = numpy.mean(FluxesDictVarLo[attrib][flux])                           
            else:
                for attrib in attribs:
                    for flux in FluxesDictVarHi[attrib]:
                        FluxesDictVarHi[attrib][flux] = 0  
                        FluxesDictVarLo[attrib][flux] = 0  
        else:
            raise Exception('Flux reference not available.')

        return FluxesDictVarHi,FluxesDictVarLo



def infoOut(input_):
    "Function that provides info for results from 13C flux fit GAMS problem"
    filename = input_[0]
    try:
        fileInfo = open(filename)
    except IOError as e:
        if e.errno == 2:
            raise Exception('Filename '+filename+' not found in directory '+os.getcwd() )
    
    infoList = []
    for line in fileInfo:
        tmp = line.split('=')
        infoList.append(float(tmp[1]))

    InfoStruct = GAMSclasses.GAMSProblemInfo(infoList[0],infoList[1],infoList[2],infoList[3],infoList[4],infoList[5],infoList[6])   
    return InfoStruct


def getInfeasibilities(input_):
    "Function that gets the infeasibilities from the .lst file in the 13C problem"
    filename = input_[0]
    infesFileName = 'infeasibilites.txt'
    # grep infeasibilities to file    
    call   = "grep INFES "+ filename +" | grep '\[' | grep _ -v | grep 'INFES =' -v > " + infesFileName  
    status = os.system(call)  # Do not remove
    #if status != 0:
    #        raise Exception('Collecting infeasibilities failed')
    # Collect infesibilites
    infeasibilitiesFile = open(infesFileName)
    infeasibilities = []
    infeasibilities.append('Infeasibilities: ')
    for line in infeasibilitiesFile:
        infeasibilities.append(line)    

    return infeasibilities


def getSuccessful(input):
    "Function that gets the result of GAMS problem running"
    filename = input[0]
    successFileName  = 'successfile.txt'
    # grep success to file
    call   = "grep 'Optimal solution found' " + filename + " > " + successFileName
    status = os.system(call)  # Do not remove
    # Check if file has content
    if os.stat(successFileName)[6] == 0:
        successful = False
    else:
        successful = True
    
    return successful
    

    
def fluxParOut(input):
    "Function that provides Flux parameter for results from GAMSProblem"
    fileName = input[0]
    parName  = input[1]
    Parameter = GAMSclasses.fluxGAMSpar(parName,[])
    Parameter.read(fileName)
    return Parameter    


def labelingParOut(input):
    "Function that provides labeling parameter"
    fileName = input[0]
    parName  = input[1]
    Parameter = GAMSclasses.labelingGAMSpar(parName,[])
    Parameter.read(fileName)
    return Parameter


############### Tests ##################


##### UNIT TESTS ########

# Expected result from tests below
expectedResult =  '153.68'  


class FBAtests(unittest.TestCase):  
    """ Testing that FBA works as expected"""

    def setUp(self):
        
        qmodeldir         = os.environ['QUANTMODELPATH']
        basedir           = qmodeldir+'/data/tests/Toya2010/2S/'
        testdir           = basedir + '2SpaperTest/'

        FBAfileName = testdir+'XML/'+'EciJR904wt5hGrowthPartial.xml'    
        FBAmodel    = FBAModel(FBAfileName)        
        
        self.FBAmodel = FBAmodel

    def testFBA(self):
        """ Checks results from FBA"""
        
        FBAresults   = self.FBAmodel.findFluxes()
        fluxDictFBA  = FBAresults.reactionNetwork.reactionList.getFluxDictionary()

        self.FBAmodel.changeObjective('BiomassEcoli',0)
        self.FBAmodel.changeObjective('PDH',1)
        FBAresults2  = self.FBAmodel.findFluxes()
        fluxDictFBA2 = FBAresults2.reactionNetwork.reactionList.getFluxDictionary()
        
        self.assertTrue(abs(fluxDictFBA['PGI'].net -4.8848561646) < 0.1)
        self.assertTrue(abs(fluxDictFBA['PDH'].net -10.2) < 0.1)
        self.assertTrue(abs(fluxDictFBA2['PDH'].net-42.5) < 0.1)

    def testFVA(self):
        """ Checks results from FVA"""
        
        FVAresults  = self.FBAmodel.FVA(['PDH','PGI','MDH'])
        fluxDictFVA = FVAresults.reactionNetwork.reactionList.getReactionDictionary()

        self.assertTrue(abs(fluxDictFVA['PDH'].flux.net.lo - 0.00) < 0.1)
        self.assertTrue(abs(fluxDictFVA['PGI'].flux.net.lo + 6.98) < 0.1)
        self.assertTrue(abs(fluxDictFVA['MDH'].flux.net.lo + 8.28) < 0.1)

        self.assertTrue(abs(fluxDictFVA['PDH'].flux.net.hi - 18.63) < 0.1)
        self.assertTrue(abs(fluxDictFVA['PGI'].flux.net.hi - 11.53) < 0.1)
        self.assertTrue(abs(fluxDictFVA['MDH'].flux.net.hi - 29.21) < 0.1)


    def tearDown(self):
        pass    


class C13MFAtests(unittest.TestCase):  
    """ Testing that 13CMFA works as expected"""

    def setUp(self):
        
        qmodeldir         = os.environ['QUANTMODELPATH']
        dirDATA           = qmodeldir+'/data/tests/TCAtoy/' 


        REACTIONSfilename   = dirDATA+'REACTIONStca.txt' 
        FEEDfilename        = dirDATA+'FEEDtca.txt'
        CEMSfilename        = dirDATA+'GCMStca.txt'
        CEMSSTDfilename     = dirDATA+'GCMSerrtca.txt'
        FLUXESFreefilename  = dirDATA+'FLUXtca.txt'

        atomTransitions = enhancedLists.AtomTransitionList(REACTIONSfilename)
        reactionNetwork = atomTransitions.getReactionNetwork('E. coli wt5h 13C MFA')

        reactionNetwork.addLabeling(CEMSfilename,'LCMSLabelData',CEMSSTDfilename,minSTD=0.001)
        reactionNetwork.addFeed(FEEDfilename)
        reactionNetwork.loadFluxBounds(FLUXESFreefilename)

        self.RNtca = reactionNetwork
        
        
        strain = 'wt5h'
        qmodeldir         = os.environ['QUANTMODELPATH']    
        dirDATA           = qmodeldir+'/data/tests/Toya2010/C13/'+strain+'/' 
    
        # Get sbml file
        REACTIONSfilename   = dirDATA + 'REACTIONS'+strain+'.txt' 
        FEEDfilename        = dirDATA + 'FEED'+strain+'.txt'
        CEMSfilename        = dirDATA + 'GCMS'+strain+'.txt'
        CEMSSTDfilename     = dirDATA + 'GCMSerr'+strain+'.txt'
        FLUXESFreefilename  = dirDATA + 'FLUX'+strain+'.txt'

        atomTransitions = enhancedLists.AtomTransitionList(REACTIONSfilename)
        reactionNetwork = atomTransitions.getReactionNetwork('E. coli wt5h 13C MFA')
        reactionNetwork.addLabeling(CEMSfilename,'LCMSLabelData',CEMSSTDfilename,minSTD=0.001)
        reactionNetwork.addFeed(FEEDfilename)
        reactionNetwork.loadFluxBounds(FLUXESFreefilename)
            
        self.RNtoya = reactionNetwork
         

    def testTCAtoy(self):
        """ Checks results from TCA toy model"""
        
        self.RNtca.write('TCA.sbml')  
        C13modelTCA = C13Model('TCA.sbml')           

        results    = C13modelTCA.findFluxesStds(Nrep=10, procString='proc', erase=False)

        self.assertTrue('Fit:' +str(results.EMUlabel['Glu']) == 'Fit:[ 0.34635  0.26953  0.27083  0.08073  0.02865  0.00391]')
        self.assertTrue('Exp:' +str(results.fragDict['Glu'].mdv) == 'Exp:[ 0.346  0.269  0.27   0.081  0.028  0.004]')
        self.assertTrue(str(results.reactionNetwork.reactionList.printFluxes(fileName='toString',brief="False")) == 
        'co2Out: \t150.0\nr6: \t125.0\nr1: \t100.0\nr2: \t100.0\nr7: \t75.0\nr4: \t50.0\nr5: \t50.0\nr3: \t50.0\nr8: \t50.0')

    def testToya2010(self):
        """ Checks results using Toya 2010 data"""
        
        names = ['2DDA7Pbm','3PHPbm','ACACCT','ALAbm','AcCoabm','CO2bm','F6Pbm','G6Pbm'
         ,'GLCpts','GLUbm','OAAbm','PEPbm','PYRbm','R5Pbm','VALbm']
 
       
        self.RNtoya.write('TOYAnew.sbml')    

        # Get FIT fluxes
        C13modelNew = C13Model('TOYAnew.sbml')   
        C13modelNew.reactionNetwork.fixFluxes(names)
        fluxNames = C13modelNew.reactionNetwork.reactionList.getReactionNameList(level=1)
        resultsNew = C13modelNew.findFluxesRanges(Nrep=10,fluxNames=fluxNames, procString='proc') 

        
        fluxesOut = 'PPC: \t31.9044364374\nPPCK: \t29.3433684566\nPGK: \t-19.6954408125\nGAPD: \t19.6954408125\nPGM: \t-18.3546408125\nENO: \t18.3546408125\nCO2bm: \t16.33\nPDH: \t12.2951639317\nGLCpts: \t11.69\nFBA: \t9.79922701254\nPFK: \t9.79922701254\nTPI: \t9.79922701254\nVALbm: \t9.3156429689\nPGI: \t9.23618253763\nACACCT: \t4.2882\nAcCoabm: \t4.0618\nCS: \t3.38238627254\nACONT: \t3.38238627254\nPEPbm: \t3.043\nICDHyr: \t2.81960861335\nMDH: \t2.77981829174\nFUM: \t2.21704063254\nSUCD1i: \t2.21704063254\nG6PDH2r: \t2.12203246237\nGND: \t2.12203246237\nPGL: \t2.12203246237\nOAAbm: \t1.9585\nSUCOAS: \t-1.65426297335\nAKGDH: \t1.65426297335\nALAbm: \t1.50843590346\nRPI: \t-1.47265998746\n3PHPbm: \t1.3408\nPGCD: \t1.3408\nGLUDy: \t-1.16534564\nGLUbm: \t1.16534564\nR5Pbm: \t0.9202743\nRPE: \t0.649372474911\nTK1: \t0.649372474911\nPYK: \t0.605173931743\nMALS: \t0.562777659198\nICL: \t0.562777659198\nTA1: \t-0.552385687456\nTA2: \t0.552385687456\nTK3: \t-0.552385687456\nDDPA: \t0.4553989\n2DDA7Pbm: \t0.4553989\nG6Pbm: \t0.331785\nTK2: \t-0.0969867874556\nF6Pbm: \t0.086328\nPYRbm: \t1e-05\nME1: \t0.0\nF6PA: \t0.0\nEDA: \t0.0\nEDD: \t0.0\nDHAPT: \t0.0'
        #self.assertTrue(resultsNew.reactionNetwork.reactionList.printFluxes(fileName='toString',brief="False") == fluxesOut)
        self.assertTrue('Fit:' +str(resultsNew.EMUlabel['fdp']) == 'Fit:[ 0.43061  0.24366  0.08718  0.0772   0.05093  0.01888  0.09155]')
        self.assertTrue('Exp:' +str(resultsNew.fragDict['fdp'].mdv) == 'Exp:[ 0.381  0.244  0.081  0.116  0.041  0.017  0.119]')  
        self.assertTrue('Fit:' +str(resultsNew.EMUlabel['pep']) == 'Fit:[ 0.62729  0.16237  0.0599   0.15045]')
        self.assertTrue('Exp:' +str(resultsNew.fragDict['pep'].mdv) == 'Exp:[ 0.624  0.165  0.06   0.151]')


" Additional Tests for 13C MFA and TS-13C MFA can found in ToyaData module, which is devoted fully to these tests "
