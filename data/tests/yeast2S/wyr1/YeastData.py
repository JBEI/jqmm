# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 17:21:03 2014

This file contains the methods to use and test the yeast test data obtained by Amit

@author: Hector Garcia Martin
"""

import SBMLclasses, FluxModels, GAMSclasses
import os
import ToyaData as TD
import shelve

inputTypes = ['Glc','Etoh']

qmodeldir         = os.environ['QUANTMODELPATH']
dirDATA           = qmodeldir+'/data/tests/yeast2S/wyr1/'
testdir           = dirDATA + 'testFiles/'



class testFiles(TD.testFiles):
    """ Testing production of SBML files """

    def setUp(self):
        # Getting reference files in dictionary form"
        files = {}
        for inputType in inputTypes:
            fileName = 'SciMM904wyr1'+inputType+'TS.xml'
            files[fileName] = GAMSclasses.file2tuple(testdir+fileName) 
            
        self.refFiles = files
           
    def testStrain(self,inputType):
        """ Tests a single inputType """
        
        newFile = getSBMLfile(inputType)        
        refFile = self.refFiles['SciMM904wyr1'+inputType+'TS.xml']
        
        self.compareFileLists([newFile],[refFile])
            
    def testQuick(self):
        # Only tests the first of the strains
        self.testStrain(inputTypes[0])
    
    def testFull(self):  
        # Compares all SBML files with previously generated version
        for inputType in inputTypes:
            self.testStrain(inputType)
            

    def refreshSavedFiles(self):
        # Refresh the files used for the comparison        
        for inputType in inputTypes:
            fileName,fileString = getSBMLfile(inputType)
            outFile = file(testdir+fileName,'w')
            outFile.write(fileString)
                    
    def tearDown(self):
        pass    
    

class testTSFluxCalculation(TD.testTSFluxCalculation):
    """ Testing flux profiles obtained from TS 13C MFA """

    def setUp(self):
        # Getting reference files in dictionary form"
        pass
        self.TSresults = {}
        self.fluxDicts = {}
        for inputType in inputTypes:
            db = shelve.open(testdir+'TSresults'+inputType)
            self.TSresults[inputType] = db['TSresult'] 
            self.fluxDicts[inputType] = db['fluxDict']   
            db.close()
        
            
    def testStrain(self,inputType):
        """ Tests a single strain """
        fluxNames = ['PGI','PGL','FBA','GAPD','PYK','CSm','AKGDbm','FUM']
      
        
        fluxNames = sorted(list(set(fluxNames))) # To avoid repetitions
        
        
        # New results
        TSresult, TSmodel = getTSresults(inputType)
        
        # Reference results
        RFresults = self.TSresults
        RFresult = RFresults[inputType]
        
        fluxDictN  = TSresult.reactionNetwork.reactionList.getReactionSubSet(rxnNames=fluxNames).getFluxDictionary(rangeComp='all')
        fluxDictRF = RFresult.reactionNetwork.reactionList.getReactionSubSet(rxnNames=fluxNames).getFluxDictionary(rangeComp='all')
        
        delta = 0.10
        errMsgAdd = '\n in strain '+inputType
        self.compareFluxDict(fluxDictRF,fluxDictN,delta,elements='lohi',errMsg=errMsgAdd)        
        
            
    def testQuick(self):
        # Only tests the first of the strains
        self.testStrain(inputTypes[0]) 
    
    def testFull(self):  
        # test all strains
        for inputType in inputTypes:
            self.testStrain(inputType)
            

    def refreshSavedFiles(self):
        # Refresh the files used for the comparison                
        for inputType in inputTypes:
            results, model = getTSresults(inputType)
            fluxDict = results.reactionNetwork.reactionList.getFluxDictionary(rangeComp='all') 
            TSresult = results  
            TSmodel  = model
            
           #Eliminate resultsDict because it is huge      
            if hasattr(TSresult,'resultsDict'):
                delattr(TSresult,'resultsDict')            
            
            db = shelve.open(testdir+'TSresults'+inputType)
            db['TSresult'] = TSresult
            db['fluxDict'] = fluxDict        
            db.close()
                                          
    def tearDown(self):
        pass        




## Functions to get data    
def getSBMLfile(inputType):
        
    REACTIONSfilename = dirDATA + 'REACTIONSwyr1'+inputType+'.txt' 
    FEEDfilename      = dirDATA + 'FEEDwyr1'+inputType+'.txt'
    CEMSfilename      = dirDATA + 'LCMSwyr1'+inputType+'.txt'    #TODO: proper name for file!!
    FLUXESfilename    = dirDATA + 'FLUXwyr1'+inputType+'.txt'
    SBMLfilename      = qmodeldir+'/data/sbmlFiles/iMM904TKs.xml'        
                
    reactionNetwork = SBMLclasses.TSReactionNetwork(SBMLfilename)
    reactionNetwork.changeModelName('SciMM904wyr1')        
        
    reactionNetwork.loadFluxBounds(FLUXESfilename)
        
    reactionNetwork.addReactions(REACTIONSfilename,translate2SBML=True)        
        
    reactionNetwork.addLabeling(CEMSfilename,'LCMS',minSTD=0.001)        
    
    reactionNetwork.addFeed(FEEDfilename)    

    file_ = ('SciMM904wyr1'+inputType+'TS.xml',reactionNetwork.write('toString'))
    
    return file_
         
def getTSmodel(inputType):
        
    SBMLfile = getSBMLfile(inputType)        
    TSmodel  = FluxModels.TwoSC13Model(SBMLfile)        
               
    return TSmodel
    
def getTSresults(inputType):      
        
    TSmodel   = getTSmodel(inputType)    

    coreFluxes = TSmodel.reactionNetwork.C13ReacNet.reactionList.getReactionNameList(level=1)
    fluxNames = [name for name in coreFluxes if not 'EX' in name]
    fluxNames = list(set(fluxNames)) # Eliminate redundancies

    TSresult = TSmodel.findFluxesRanges(Nrep=30,fluxNames=fluxNames)  
                        
    return TSresult,TSmodel    