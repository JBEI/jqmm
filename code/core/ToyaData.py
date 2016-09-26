# Part of the JBEI Quantitative Metabolic Modeling Library (JQMM)
# Copyright (c) 2016, The Regents of the University of California.
# For licensing details see "license.txt" and "legal.txt".
"""
All functions to obtain the Toya data and a test suite.
This module effectively reproduces all figures in Garcia Martin et al 2015
"""

import unittest, os, copy, shelve, sys, traceback, re
import matplotlib, numpy, difflib
import matplotlib.pyplot as plt
import core, FluxModels, predictions, ReactionNetworks, enhancedLists
from pylab import figure, title, savefig, plot, xlabel, ylabel
from IPython.display import SVG
import utilities as utils

### Needed paths
qmodeldir         = os.environ['QUANTMODELPATH']
basedir           = qmodeldir+'/data/tests/Toya2010/2S/'
testdir           = basedir + '2SpaperTest/'

### 
# Strains for which data is available in 
strains    = ['wt5h','wt6h','wt7h','pyk5h','pyk6h','pyk7h','pgi16h','pgi21h','pgi23h']
# Dictionary that maps the base strain to make predictions from (i.e. pyk5h MOMA predictions are based of wt5h)
strainBase = {'pyk5h':'wt5h','pyk6h':'wt6h','pyk7h':'wt7h','pgi16h':'wt5h','pgi21h':'wt6h','pgi23h':'wt7h'} 

class generalTest(unittest.TestCase):
    """ Class with some general methods of use to all other tests"""
    def compareFileLists(self,files1,files2):
        """
        Compares lists of files to check whether they are identical or not
        files1 and files2 are expected to be lists containing files to be compared:
        
        files1 = [file1, file2, file3, .... etc]
        
        Files in the list can be given in the form of filenames or tuples for the form (filename, filestring)
        Alternatively, the inputs can be dicts:
        
        files1 = {'filename1':file1, 'filename2':file2, 'filename3':file3... etc}
        
        """
        # Put file lists in form of dictionaries for easier manipulation if needed
        dicts = []
        for thing in [files1,files2]: 
            # if list convert to dict
            if type(thing) is list:   
                fileList = thing
                dict_= {}
                for file_ in fileList:
                    (name,string) = utils.file2tuple(file_)
                    dict_[name]=string.splitlines()
            # if dict keep as it is
            elif type(thing) is dict: 
                dict_ = thing
            else:              # unknown input error  
                print files1
                print files2
                raise Exception('Inputs must be in list or dict type')
            dicts.append(dict_)

        dict1 = dicts[0]
        dict2 = dicts[1]                
        
        # See if we have the same file names on each list
        names1 = set(dict1.keys())
        names2 = set(dict2.keys())
        errmsg1 = "File list is different: \n"+str(names1)+" vs \n"+str(names2)
        self.assertEqual(names1,names2,msg=errmsg1)
        
        # Make comparison of file content
        for fileName in dict1:
            diff = difflib.unified_diff(dict1[fileName],dict2[fileName],n=10)  # TODO: figure out how to get a better output out of the diff
            diffOut = '\n'.join(list(diff))
            errmsg2 = "New version of files "+fileName+" is different\n"+"diff:\n"+diffOut
            self.assertEqual(diffOut,'',msg=errmsg2)

   
    def compareFluxDict(self,dict1,dict2,delta,elements='all',errMsg=''):
        """
        Compares two FLUX dictionaries of the type (e.g.):
        
        fluxDict['PGI'] = flux class instance
        
        -)delta is the relative change permitted
        -)Elements decides which elements to compare (lo, best or hi) from the ranged number in core 
        -)errMsg is the additional string to add to the error messages
        """

        # See if we have the same fluxes on each dict
        errmsg1 = "Fluxes are different:\n"+str(dict1.keys())+"\n"+str(dict2.keys())
        self.assertEqual(set(dict1.keys()),set(dict2.keys()),msg=errmsg1)

        noChangeTotal = True
        errmsg        = '\n'
        for name in dict1.keys():
            fluxRF = dict1[name]
            fluxN  = dict2[name]
            
            # Test all values
            bestChanged = changed(fluxN.net.best,fluxRF.net.best,delta)
            loChanged   = changed(fluxN.net.lo,fluxRF.net.lo,delta)
            hiChanged   = changed(fluxN.net.hi,fluxRF.net.hi,delta)
            
            if elements == 'all':
                noChange = (not bestChanged) and (not loChanged) and (not hiChanged)
            elif elements == 'lohi':
                noChange = (not loChanged) and (not hiChanged)
            elif elements == 'best':
                noChange = (not bestChanged)                 

            noChangeTotal = noChangeTotal and noChange
            if not noChange:
                errmsg = errmsg + 'Flux '+name+' changed from: \n'+str(fluxRF.net)+' to: \n'+ str(fluxN.net)+'\n'

        self.assertTrue(noChangeTotal,msg=errmsg) 
    
    def testFull(self):
        "If not defined, defined as testQuick by default"
        self.testQuick()



class testInputFiles(generalTest):  # TODO: change this to testInputFiles
    """ Testing production of SBML files used as input for everything else """

    def setUp(self):
        # Getting reference files in dictionary form"
        try:
            files = {}
            for strain in strains:
                fileName = 'EciJR904TKs'+strain+'TS.xml'   
                files[fileName] = utils.file2tuple(testdir+'XML/'+fileName) 
            self.refFiles = files
        except:
            e = sys.exc_info()[0]
            print 'Problems loading stored data in testFiles: \n'+str(e)   
           
           
    def testStrain(self,strain):
        """ Tests a single strain by comparing with reference stored files"""
        
        newFile = getSBMLfile(strain)        
        refFile = self.refFiles['EciJR904TKs'+strain+'TS.xml']
        
        self.compareFileLists([refFile],[newFile])


    def testQuick(self):
        """ Only tests the first of the strains """
        self.testStrain(strains[0])


    def testFull(self):  
        """ Tests all of the strains """
        for strain in strains:
            self.testStrain(strain)
            

    def refreshSavedFiles(self):
        """ Refresh the files used as reference for the comparison """        
        for strain in strains:
            fileName,fileString = getSBMLfile(strain)
            outFile = file(testdir+'XML/'+fileName,'w')
            outFile.write(fileString)


    def tearDown(self):
        pass    



class testOptimizationFiles(generalTest):
    """ Testing input files for TS optimization problem """

    def setUp(self):
        # Getting reference files in dictionary form"
        try:
            files = {}
            for strain in strains:
                files[strain] = readFiles(testdir+strain+'/') 
            self.refFiles = files
        except:
            e = sys.exc_info()[0]
            print 'Problems loading stored data in testFiles: \n'+str(e)   


    def testStrain(self,strain):
        """ Tests a single strain by comparing with reference stored files"""
        
        TSmodel  = getTSmodel(strain)
        newFiles = [filterFiles(x) for x in TSmodel.reactionNetwork.getTwoS13CMFAfiles()] 
        refFiles = [filterFiles(y) for y in sorted(self.refFiles[strain]) if filterFiles(y)!=None]
        
        self.compareFileLists(sorted(refFiles),sorted(newFiles))


    def testQuick(self):
        """ Only tests the first of the strains """
        self.testStrain(strains[0])


    def testFull(self):  
        """ Tests all of the strains """
        for strain in strains:
            self.testStrain(strain)
            

    def refreshSavedFiles(self):
        """ Refresh the files used as reference for the comparison """        
        for strain in strains:
            TSmodel = getTSmodel(strain)
            files   = TSmodel.reactionNetwork.getTwoS13CMFAfiles() 
            for file_ in files:
                fileName,fileString = file_
                outFile = file(testdir+strain+'/'+fileName,'w')
                outFile.write(fileString)
                    
    def tearDown(self):
        pass        



class testTSFluxCalculation(generalTest):
    """ Testing flux profiles obtained from TS 13C MFA """

    def setUp(self):
        """ Getting reference files in dictionary form """
        self.TSresults = {}
        self.fluxDicts = {}
        for strain in strains:
            try:
                db = shelve.open(testdir+strain+'/'+'TSresults'+strain)
                self.TSresults[strain] = db['TSresult'] 
                self.fluxDicts[strain] = db['fluxDict']   
                db.close()
            except:
                print 'Problems loading stored data in testTSFluxCalculation for strain %s:\n%s' % (strain, traceback.format_exc())


    def testStrain(self,strain,saved=False):
        """ Tests a single strain 
        -) saved decides whether to used saved solutions (if True) or calculate them anew (if False) 
        """
        print '\nTesting strain: '+str(strain)+'\n'
        # Names for fluxes to be tested
        fluxNames      = ['G6PDH2r','GAPD','CS','FUM','ICDHyr','ACODA','ACKr','PDH']
        fluxNames.extend(['THD2','ICDHyr','G6PDH2r','GND','MTHFD','GLUDy','KARA1i','C160SN','C181SN','ASAD'])
        fluxNames.extend(['MDH','PGCD','AKGDH','NADH6'])
        fluxNames.extend(['EX_h2o_e_','EX_o2_e_','EX_co2_e_','EX_nh4_e_','EX_h_e_','EX_glc_e_','EX_ac_e_','EX_pi_e_','EX_urea_e_'])
        fluxNames.extend(['EX_so4_e_','EX_glyclt_e_','EX_acald_e_','EX_fum_e_'])
        # Avoid repetitions
        fluxNames = sorted(list(set(fluxNames))) 
           
        # Obtain results
        TSresult, TSmodel = getTSresults(strain,saved=saved)
        
        # Reference results
        RFresults = self.TSresults
        RFresult = RFresults[strain]
        
        fluxDictN  = TSresult.reactionNetwork.reactionList.getReactionSubSet(rxnNames=fluxNames).getFluxDictionary(rangeComp='all')  # new results
        fluxDictRF = RFresult.reactionNetwork.reactionList.getReactionSubSet(rxnNames=fluxNames).getFluxDictionary(rangeComp='all')  # reference results
        
        # Flux comparison
        delta = 0.10
        errMsgAdd = '\n in strain '+strain        
        self.compareFluxDict(fluxDictRF,fluxDictN,delta,elements='lohi',errMsg=errMsgAdd) # Comparing only lowest and highest flux values       


    def testQuick(self):
        """ Only tests the first of the strains """
        self.testStrain(strains[0]) 


    def testFull(self):  
        """ Test all strains """
        AllCasesWork = True
        errmsg       = ''
        
        # Run tests and compile all failures by catching the exceptions
        failures = []
        for strain in strains:
            try:
                self.testStrain(strain)
            except AssertionError as e:
                failures.append("\n"+"Strain:"+str(strain)+"\n"+str(e)+"\n")
        
        # Raise failure if any fail but give report on all of them
        if failures:
            errmsg = "\n".join([str(x) for x in failures])
            AllCasesWork = False
            
        self.assertTrue(AllCasesWork,msg=errmsg)


    def refreshSavedFiles(self):
        """ Refresh the files used for the comparison """               
        for strain in strains:
            results, model = getTSresults(strain)
            fluxDict = results.reactionNetwork.reactionList.getFluxDictionary(rangeComp='all') 
            TSresult = results  

            #Eliminate resultsDict because it is huge      
            if hasattr(TSresult,'resultsDict'):
                delattr(TSresult,'resultsDict')
            
            db = shelve.open(testdir+strain+'/'+'TSresults'+strain)
            db['TSresult'] = TSresult
            db['fluxDict'] = fluxDict        
            db.close()
                                          
    def tearDown(self):
        pass        



class testELVA(generalTest):
    """ Tests ELVA (External Labeling Variability Analysis) results """

    def setUp(self):
        """ Getting reference files in dictionary form """
        try:
            db = shelve.open(testdir+'ELVArefs')
            self.labelCompDict1 = db['labelCompDict1'] 
            self.labelCompDict2 = db['labelCompDict2']   
            db.close()
        except:
            e = sys.exc_info()[0]
            print 'Problems loading stored data in testELVA: \n'+str(e)  
            
            
    def getResults(self,saved=False):
        """ Results to use for the comparison """
        strain = 'wt5h'        

        ### Obtain ELVA results used for the comparison                
        # ELVA with new reaction set
        TSresult, TSmodel = getTSresults(strain,saved=saved)
        resultELVA  = TSmodel.ELVA(TSresult,procString='proc')  
         
        # ELVA with old reaction set
        TSresultOld, TSmodelOld = getTSresults(strain,oldVersion=True)
        resultELVAOld  = TSmodelOld.ELVA(TSresultOld,procString='proc')     

        ### Obtain labeling component dictionaries to be compared
        x1,Dx1,y1,DyU1,DyL1,labelComp1,labelCompDict1 = resultELVA.getxDXyDyInfo()
        x2,Dx2,y2,DyU2,DyL2,labelComp2,labelCompDict2 = resultELVAOld.getxDXyDyInfo()
        
        return labelCompDict1,labelCompDict2


    def compLabelCompDict(self,dict1,dict2,errAdd):
        """
        Compares labelCompDicts. Assuming dict1 to be reference 
        -) errAdd is the extra string added at the beginning of the error string      
        """
        # Acceptable relative variation fraction
        delta = 0.1 
        
        # See if we have the same files on each list
        errmsg1 = "Fragment lists keys to be compared are different: \n" +str(dict1.keys())+'\n'+'vs.'+str(dict2.keys())
        self.assertEqual(set(dict1.keys()),set(dict2.keys()),msg=errmsg1)
        
        for mvalue in dict1:
            best1,lo1,hi1 = dict1[mvalue]
            best2,lo2,hi2 = dict2[mvalue]

            bestChanged = changed(best1,best2,delta)
            loChanged   = changed(lo1,lo2,delta)
            hiChanged   = changed(hi1,hi2,delta)
            
            noChange = (not bestChanged) and (not loChanged) and (not hiChanged)
            errmsg = errAdd+' M value '+mvalue+' changed from '+','.join([str(best1),str(lo1),str(hi1)]) +' to '+ ','.join([str(best2),str(lo2),str(hi2)])
            self.assertTrue(noChange,msg=errmsg)        


    def testQuick(self):
        """ Quick test for wt5h only """
        # Obtain reference labeling dictionaries
        labelCompDictRef    = self.labelCompDict1
        labelCompDictOldRef = self.labelCompDict2

        # Obtain new labeling dictionaries
        labelCompDict,labelCompDictOld = self.getResults()

        # Error checks
        errAdd = "In new ELVA:"
        self.compLabelCompDict(labelCompDictRef,labelCompDict,errAdd)    
        errAdd = "In old ELVA:"
        self.compLabelCompDict(labelCompDictOldRef,labelCompDictOld,errAdd)
        
    
    def refreshSavedFiles(self):
        """ Refreshes saved files used for comparison """
        
        labelCompDict1,labelCompDict2 = self.getResults(saved=True)

        db = shelve.open(testdir+'ELVArefs')
        db['labelCompDict1'] = labelCompDict1
        db['labelCompDict2'] = labelCompDict2
        db.close()



class testConstPower(generalTest):
    """ Test relative constraining power of different approaches (2S, FVA, FVA with 2S approximation) """
    
    def setUp(self):
        """ Getting reference saved reaction lists """
        try:
            db = shelve.open(testdir+'ConstPowerRefs')
            self.TSreactionList       = db['TSreactionList'] 
            self.constFVAreactionList = db['constFVAreactionList']   
            self.FVAreactionList      = db['FVAreactionList']   
            db.close()
        except:
            print 'Problems loading stored data in testConstPower:\nTest dir: %s:\n%s' % (testdir, traceback.format_exc())
 
 
    def getResults(self,strain='wt5h',saved=True):
        """ New results to be compared with saved references """ 
        ## TS calculations          
        TSresult, TSmodel = getTSresults(strain,saved=saved)
        fluxNames = getFluxNames(TSmodel)

        ## FVA calculations
        SBMLfile = getSBMLfile(strain)
        FVAmodel   = FluxModels.FBAModel(SBMLfile)
        FVAresults = FVAmodel.FVA(reactions=fluxNames)        

        ## FVA with 2S constraints
        SBMLfile = TSresult.reactionNetwork.write('toString')
        constFVAModel   = FluxModels.FBAModel(SBMLfile)
        constFVAresults = constFVAModel.FVA(reactions=fluxNames)         

        # Getting reactions lists from results
        TSreactionList       = TSresult.reactionNetwork.reactionList
        constFVAreactionList = constFVAresults.reactionNetwork.reactionList
        FVAreactionList      = FVAresults.reactionNetwork.reactionList

        # Providing a name for each reaction list (for plotting purposes)
        TSreactionList.name       = '2S-$^{13}$C MFA'
        constFVAreactionList.name = '2S FVA'
        FVAreactionList.name      = 'FVA'

        return TSreactionList, constFVAreactionList, FVAreactionList


    def testQuick(self):
        """ Quick test for wt5h only """
        
        # Reference reaction lists
        TSreactionListRF       = self.TSreactionList
        constFVAreactionListRF = self.constFVAreactionList
        FVAreactionListRF      = self.FVAreactionList

        # Newly calculated reaction lists
        TSreactionListN, constFVAreactionListN, FVAreactionListN = self.getResults(saved=False)

        # Names of fluxes to be compared
        fluxNames = ['G6PDH2r','GND','PGL','RPE','TK1','TA2','EDA','EDD','TK2','TA1','TK3','RPI']

        # Limiting reactions lists to desired fluxes as given per fluxNames above and obtaining flux dictionaries
        TSfluxDictRF = TSreactionListRF.getReactionSubSet(rxnNames=fluxNames).getFluxDictionary(rangeComp='all')
        TSfluxDictN  = TSreactionListN.getReactionSubSet(rxnNames=fluxNames).getFluxDictionary(rangeComp='all')
        
        constFVAfluxDictRF = constFVAreactionListRF.getReactionSubSet(rxnNames=fluxNames).getFluxDictionary(rangeComp='all')
        constFVAfluxDictN  = constFVAreactionListRF.getReactionSubSet(rxnNames=fluxNames).getFluxDictionary(rangeComp='all')
        
        FVAfluxDictRF = FVAreactionListRF.getReactionSubSet(rxnNames=fluxNames).getFluxDictionary(rangeComp='all')
        FVAfluxDictN  = FVAreactionListRF.getReactionSubSet(rxnNames=fluxNames).getFluxDictionary(rangeComp='all')
        
        # Flux value comparisons
        delta = 0.05
        self.compareFluxDict(TSfluxDictRF,TSfluxDictN,delta,elements='lohi')
        self.compareFluxDict(constFVAfluxDictRF,constFVAfluxDictN,delta,elements='lohi')
        self.compareFluxDict(FVAfluxDictRF,FVAfluxDictN,delta,elements='lohi')
        
        
    def refreshSavedFiles(self):   
        """ Refresh saved files used as reference """
        TSreactionList, constFVAreactionList, FVAreactionList  = self.getResults(saved=False)

        db = shelve.open(testdir+'ConstPowerRefs')
        db['TSreactionList']       = TSreactionList
        db['constFVAreactionList'] = constFVAreactionList
        db['FVAreactionList']      = FVAreactionList
        db.close()    



class testPredictions(generalTest):
    """ Test flux and labeling predictions through COBRA methods """    

    def setUp(self):
        """ Getting saved reference files """

        try:
            db = shelve.open(testdir+'PredictionRefs')
            #self.resultsV1 = db['resultsV1']          
            #self.resultsV2 = db['resultsV2']           
            self.TSMOMA    = db['TSMOMA']
            self.TS13CMOMA = db['TS13CMOMA']
            self.results13CMOMAFixed = db['results13CMOMAFixed']
            self.resultsMOMAFixed    = db['resultsMOMAFixed']
            db.close()
        except:
            print 'Problems loading stored data in testPredictions:\n%s' % (traceback.format_exc())


    def getFVAresults(self,strain):
        """
        Produces FVA results for predictions        
        """
        ### FBA results
        TSmodel   = getTSmodel(strain)    
        fluxNames = getFluxNames(TSmodel)
        # TODO: produce the initial files anew here
        #dirDATAsave     = basedir+strain+'/'+'testFiles/'
        dirDATAsave     = testdir+'XML/'
        
        fileNameGP = dirDATAsave+'EciJR904'+strain+'Growth'+'Partial'+'.xml'
        fileNameGF = dirDATAsave+'EciJR904'+strain+'Growth'+'Full'+'.xml'
        fileNameAP = dirDATAsave+'EciJR904'+strain+'Atp'+'Partial'+'.xml'
        fileNameAF = dirDATAsave+'EciJR904'+strain+'Atp'+'Full'+'.xml'
        
        FBAmodelGP = FluxModels.FBAModel(fileNameGP)
        FBAmodelGF = FluxModels.FBAModel(fileNameGF)
        FBAmodelAP = FluxModels.FBAModel(fileNameAP)
        FBAmodelAF = FluxModels.FBAModel(fileNameAF)

        resultsVGP = FBAmodelGP.FVA(reactions=fluxNames)   # Growth maximization
        resultsVGF = FBAmodelGF.FVA(reactions=fluxNames)   # Growth maximization
        resultsVAP = FBAmodelAP.FVA(reactions=fluxNames)   # ATP maximization
        resultsVAF = FBAmodelAF.FVA(reactions=fluxNames)   # ATP maximization
        
        return resultsVGP,resultsVGF,resultsVAP,resultsVAF


    def getFluxPredResults(self,strain):
        """
        Produces MOMA and ROOM flux prediction results
        """
        ### Calculating base flux profiles
        # Standard MOMA
        FBAfileName = testdir+'XML/'+'EciJR904'+strainBase[strain]+'Growth'+'Partial'+'.xml'
        FBAmodel = FluxModels.FBAModel(FBAfileName)

        FBAresults = FBAmodel.findFluxes()

        TSmodel = getTSmodel(strainBase[strain])
        reactionNetwork = copy.deepcopy(TSmodel.reactionNetwork)

        # TS 13C MOMA      
        TSresultsRand = TSmodel.findFluxesStds(Nrep=30,Nrand=10)   # this takes a while

        ### Predictions
        # Change flux bounds to avoid biasing prediction
        reactionNetwork.changeFluxBounds('GLCpts'      ,core.fluxBounds( 0, 25  ,False)[1])      
        reactionNetwork.changeFluxBounds('EX_glc_e_'   ,core.fluxBounds(-15,-0  ,True,True)[1])
        reactionNetwork.changeFluxBounds('BiomassEcoli' ,core.fluxBounds( 0, 25  ,False)[1])      
        reactionNetwork.changeFluxBounds('EX_ac_e_'    ,core.fluxBounds( 0, 25  ,True,True)[1])
            
        # Find reactions to knock out
        KO = re.sub('\d+h','',strain).upper()   # 'pyk5h' --> 'PYK'
            
        TSMOMA    = predictions.predict(FBAresults   , KO, 'MOMA', reactionNetwork.getSBMLString())   # TODO: find better way to store these predictions (in reactionList?)
        TS13CMOMA = predictions.predict(TSresultsRand, KO, 'MOMA', reactionNetwork.getSBMLString())  
        TSROOM    = predictions.predict(FBAresults   , KO, 'ROOM', reactionNetwork.getSBMLString())   
        TS13CROOM = predictions.predict(TSresultsRand, KO, 'ROOM', reactionNetwork.getSBMLString())  
        
        return TSMOMA,TS13CMOMA,TSROOM,TS13CROOM


    def getLabelPredResults(self,TS13CMOMA,TSMOMA):   
        """
        Produces labeling predictions
        """
        ## Labeling predictions
        # 13CMOMA
        TS13CMOMAFix = copy.deepcopy(TS13CMOMA)
        TS13CMOMAFix.reactionNetwork.reactionList.fixFluxes(level=1,epsilon=0.01,maximum=900)        
        TS13CMOMAFix.reactionNetwork.C13ReacNet.reactionList.fixFluxes(level=2,epsilon=0.01) 
        # 13CROOM   # leave in here for possible future use   
        #TS13CROOMFix = copy.deepcopy(TS13CROOM)
        #TS13CROOMFix.reactionNetwork.reactionList.fixFluxes(level=1,epsilon=0.01,maximum=900)        
        #TS13CROOMFix.reactionNetwork.C13ReacNet.reactionList.fixFluxes(level=2,epsilon=0.01) 
        
        # MOMA
        TSMOMAFix = copy.deepcopy(TSMOMA)
        TSMOMAFix.reactionNetwork.C13ReacNet.reactionList.zeroFluxExchange()
        TSMOMAFix.reactionNetwork.reactionList.fixFluxes(level=1,epsilon=0.01,maximum=900) 
        TSMOMAFix.reactionNetwork.C13ReacNet.reactionList.fixFluxes(level=2,epsilon=0.01) 
        # ROOM
        #TSROOMFix = copy.deepcopy(TSROOM)
        #TSROOMFix.reactionNetwork.C13ReacNet.reactionList.zeroFluxExchange()
        #TSROOMFix.reactionNetwork.reactionList.fixFluxes(level=1,epsilon=0.01,maximum=900) 
        #TSROOMFix.reactionNetwork.C13ReacNet.reactionList.fixFluxes(level=2,epsilon=0.01) 
    
        # Calculating labeling
        results13CMOMAFixed = TS13CMOMAFix.findFluxesStds(Nrep=1,Nrand=0,limitFlux2Core=False,maxFlux13C=300) 
        resultsMOMAFixed    = TSMOMAFix.findFluxesStds(Nrep=1,Nrand=0,limitFlux2Core=False,maxFlux13C=300)         
        #results13CROOMFixed = TS13CROOMFix.findFluxesStds(Nrep=1,Nrand=0,limitFlux2Core=False,maxFlux13C=300) 
        #resultsROOMFixed    = TSROOMFix.findFluxesStds(Nrep=1,Nrand=0,limitFlux2Core=False,maxFlux13C=300)       
        
#        # Alternative labeling calculation
#        holderModel1   = getTSmodel('wt5h') 
#        holderModel1.reactionNetwork = TS13CMOMAFix.reactionNetwork
#        results13CMOMAFixed = holderModel1.findFluxesStds(Nrep=1,Nrand=0,limitFlux2Core=False,maxFlux13C=300) 
#
#        holderModel2   = getTSmodel('wt5h') 
#        holderModel2.reactionNetwork = TSMOMAFix.reactionNetwork
#        resultsMOMAFixed = holderModel2.findFluxesStds(Nrep=1,Nrand=0,limitFlux2Core=False,maxFlux13C=300) 
        
        return results13CMOMAFixed,resultsMOMAFixed


    def getFluxNames(self):
        """
        Obtain flux names for plots in supp. figures 17 and 18
        """
        
        # Finding order for fluxes (fluxes for first strain from highest to lowest)
        subsystems  = ['S_Pentose_Phosphate_Pathway','S_GlycolysisGluconeogenesis','S_Citric_Acid_Cycle','S_OUT']
        TSresult, TSmodel = getTSresults('wt5h',saved=True)    
        reactionListRef   = TSresult.reactionNetwork.reactionList
        first   = enhancedLists.ReactionList([])
        for subsys in subsystems:
            subReactionList = reactionListRef.getReactionSubSet(rxnNames=subsys)
            subReactionList.sort('by net flux')
            first.extend(subReactionList)
        rxnNames = first.getReactionNameList(level=1,sort=False)      
        
        return rxnNames


    def getFullCompFig(self, strain, saved=True):
        """
        Obtain figure of full flux prediction comparisons (Fig S17)
        """
        subsystems  = ['S_Pentose_Phosphate_Pathway','S_GlycolysisGluconeogenesis','S_Citric_Acid_Cycle','S_OUT']
        rxnNames3   = self.getFluxNames()

        # Get results
        TSresult, TSmodel = getTSresults(strain,saved=saved)    
        resultsVGP,resultsVGF,resultsVAP,resultsVAF = self.getFVAresults(strain)
        TSMOMA,TS13CMOMA,TSROOM,TS13CROOM           = self.getFluxPredResults(strain)
            
        # Get reaction lists
        TSreactionList     = TSresult.reactionNetwork.reactionList.getReactionSubSet(rxnNames=rxnNames3)
        FVAreactionList    = resultsVGF.reactionNetwork.reactionList.getReactionSubSet(rxnNames=rxnNames3)
        ATPreactionList    = resultsVAF.reactionNetwork.reactionList.getReactionSubSet(rxnNames=rxnNames3)

        TSreactionList.name     = '2S-$^{13}$C MFA'    
        FVAreactionList.name    = 'FBA max. growth'
        ATPreactionList.name    = 'FBA max. ATP'
        
        C13MOMAreactionList  = TS13CMOMA.reactionNetwork.reactionList.getReactionSubSet(rxnNames=rxnNames3)
        C13ROOMreactionList  = TS13CROOM.reactionNetwork.reactionList.getReactionSubSet(rxnNames=rxnNames3)
        C13MOMAreactionList.name = '$^{13}$C MOMA'
        C13ROOMreactionList.name = '$^{13}$C ROOM'

        MOMAreactionList  = TSMOMA.reactionNetwork.reactionList.getReactionSubSet(rxnNames=rxnNames3)
        ROOMreactionList  = TSROOM.reactionNetwork.reactionList.getReactionSubSet(rxnNames=rxnNames3)
        MOMAreactionList.name = 'MOMA'
        ROOMreactionList.name = 'ROOM'

        # Plot fluxes
        allLists = [TSreactionList, FVAreactionList, ATPreactionList, C13MOMAreactionList, C13ROOMreactionList, MOMAreactionList, ROOMreactionList]
        utils.fluxPlotGraph(allLists,subsystems,titleFig=strain,save='compFull'+strain+'A.png',plotType=1,legendFontSize=17)  


    def getPartCompFig(self, strain, saved=True):
        """
        Obtain figure of partial flux prediction comparisons  (Fig S18)
        """
        subsystems  = ['S_Pentose_Phosphate_Pathway','S_GlycolysisGluconeogenesis','S_Citric_Acid_Cycle','S_OUT']
        rxnNames3   = self.getFluxNames()

        # Get results
        TSresult, TSmodel = getTSresults(strain,saved=saved)    
        resultsVGP,resultsVGF,resultsVAP,resultsVAF = self.getFVAresults(strain)
              
        # Get reaction lists
        TSreactionList     = TSresult.reactionNetwork.reactionList.getReactionSubSet(rxnNames=rxnNames3)
        FVAreactionList    = resultsVGP.reactionNetwork.reactionList.getReactionSubSet(rxnNames=rxnNames3)
        ATPreactionList    = resultsVAP.reactionNetwork.reactionList.getReactionSubSet(rxnNames=rxnNames3)
        
        TSreactionList.name     = '2S-$^{13}$C MFA'    
        FVAreactionList.name    = 'FBA max. growth'
        ATPreactionList.name    = 'FBA max. ATP'
                              
        # Plot fluxes
        utils.fluxPlotGraph([TSreactionList, FVAreactionList, ATPreactionList], subsystems, titleFig=strain, save='compPart'+strain+'A.png', plotType=1)        


    def testQuick(self):
        """ Quick test for pyk5h only """
        strain ='pyk5h'

        # Reference data
        TS13CMOMARF           = self.TS13CMOMA      
        results13CMOMAFixedRF = self.results13CMOMAFixed 

        # New data
        #resultsV1,resultsV2,TSMOMA,TS13CMOMA,results13CMOMAFixed,resultsMOMAFixed = self.getResults()  
        TSMOMA,TS13CMOMA,TSROOM,TS13CROOM           = self.getFluxPredResults(strain)           
        results13CMOMAFixed,resultsMOMAFixed        = self.getLabelPredResults(TS13CMOMA,TSMOMA)

        # Flux sets for desired flux names
        fluxNames = ['G6PDH2r','GND','PGL','RPE','TK1','TA2','EDA','EDD','TK2','TA1','TK3','RPI']
        
        fluxDictRF = TS13CMOMARF.reactionNetwork.reactionList.getReactionSubSet(rxnNames=fluxNames).getFluxDictionary(rangeComp='all')
        fluxDictN  = TS13CMOMA.reactionNetwork.reactionList.getReactionSubSet(rxnNames=fluxNames).getFluxDictionary(rangeComp='all')
        
        # Flux comparison
        delta = 0.15
        self.compareFluxDict(fluxDictRF,fluxDictN,delta,elements='best')      
        
        # Labeling comparison
        EMUlabelRF = results13CMOMAFixedRF.EMUlabel
        EMUlabel   = results13CMOMAFixed.EMUlabel
        
        delta = 0.05
        for name in EMUlabel:
            noChange = True
            MDVRF = EMUlabelRF[name]
            mdv   = EMUlabel[name]
            for n in range(len(mdv)):
                noChange = noChange and (not changed(MDVRF[n],mdv[n],delta,relative=False))
                
            errmsg = 'Fragment '+name+' changed from '+str(MDVRF)+' to '+ str(mdv)
            self.assertTrue(noChange,msg=errmsg)         


    def refreshSavedFiles(self):    
        """ Refreshes saved files used for comparison """
        strain ='pyk5h'
        #resultsV1,resultsV2,TSMOMA,TS13CMOMA,results13CMOMAFixed,resultsMOMAFixed = self.getResults()  
        TSMOMA,TS13CMOMA,TSROOM,TS13CROOM           = self.getFluxPredResults(strain)           
        results13CMOMAFixed,resultsMOMAFixed        = self.getLabelPredResults(TS13CMOMA,TSMOMA)

        db = shelve.open(testdir+'PredictionRefs')
        db['TSMOMA']              = TSMOMA
        db['TS13CMOMA']           = TS13CMOMA
        db['results13CMOMAFixed'] = results13CMOMAFixed
        db['resultsMOMAFixed']    = resultsMOMAFixed
        db.close()        



############# Auxliary functions to standarized and ease obtaining Toya data and objects #################   

def getSBMLfile(strain,oldVersion=False):
    """
    Provides SBML files for given strain
    -) oldVersion indicates whether to use the original version of transition reactions (if True) or the version obtained after the ELVA (if False)
    """
        
    datadir = basedir+strain+'/'
            
    BASEfilename      = datadir + 'EciJR904TKs.xml'
    FLUXESfilename    = datadir + 'FLUX'+strain+'.txt'
    REACTIONSfilename = datadir + 'REACTIONS'+strain+'.txt'      
    MSfilename        = datadir + 'GCMS'+strain+'.txt'
    FEEDfilename      = datadir + 'FEED'+strain+'.txt'
    MSSTDfilename     = datadir + 'GCMSerr'+strain+'.txt'
            
    if oldVersion:
        REACTIONSfilename = datadir + 'REACTIONS'+strain+'Old.txt'      
                            
    reactionNetwork = getTSReacNet(BASEfilename,FLUXESfilename,REACTIONSfilename,MSfilename,MSSTDfilename,FEEDfilename)

    file_ = ('EciJR904TKs'+strain+'TS.xml',reactionNetwork.write('toString'))

    return file_


def getTSmodel(strain,oldVersion=False):
    """ Provides two-scale model for given strain
    -) oldVersion indicates whether to use the original version of transition reactions (if True) or the version obtained after the ELVA (if False)    
    """
        
    SBMLfile = getSBMLfile(strain,oldVersion=oldVersion)        
    TSmodel  = FluxModels.TwoSC13Model(SBMLfile)        
               
    return TSmodel


def getTSresults(strain,saved=False,oldVersion=False): 
    """
    Provides two-scale results and model for specified strain
    -) saved decides whether to provide saved results (if True) or to calculate them anew (if False)
    -) oldVersion decides wehther to use the initial set of reactions (if True) or the reactions determined after the ELVA (if False)    
    """    
        
    if not saved:  
        TSmodel   = getTSmodel(strain,oldVersion=oldVersion)    
        fluxNames = getFluxNames(TSmodel)
        TSresult = TSmodel.findFluxesRanges(Nrep=30,fluxNames=fluxNames,procString='proc')  
    else:   # get saved version (quicker for tests)
        if oldVersion:
            raise Exception('Old version of models not saved')
        else:
            db = shelve.open(testdir+strain+'/'+'TSresults'+strain)
            TSresult = db['TSresult'] 
            TSmodel   = getTSmodel(strain,oldVersion=oldVersion)    
            db.close()               
        
    return TSresult,TSmodel    


def get13CMFAresults(strain): 
    """
    Provides 13CMFA results and model for specified strain.   
    """    
        
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

    names = ['2DDA7Pbm','3PHPbm','ACACCT','ALAbm','AcCoabm','CO2bm','F6Pbm','G6Pbm','GLCpts','GLUbm','OAAbm','PEPbm','PYRbm','R5Pbm','VALbm']
    reactionNetwork.write('TOYA.sbml')    

    # Get FIT fluxes
    C13model = FluxModels.C13Model('TOYA.sbml')   
    C13model.reactionNetwork.fixFluxes(names)
    results    = C13model.findFluxesRanges(Nrep=10,fluxNames=C13model.reactionNetwork.reactionList.getReactionNameList(level=1), procString='proc') 

    return results,C13model


def getFluxNames(TSmodel):
    """ 
    Gives names for a set of interesting fluxes to look at 
    """
    coreFluxes = TSmodel.reactionNetwork.C13ReacNet.reactionList.getReactionNameList(level=1)
    fluxNames = [name for name in coreFluxes if 'EX' not in name]    
 
    fluxNames.append('THD2')
    fluxNames.append('MTHFD')
    fluxNames.append('ACt2r')
    fluxNames.append('ACCOACr')

    fluxNames.extend(TSmodel.reactionNetwork.reactionList.getReactionSubSet('exchange').getReactionNameList(level=1))
    fluxNames.extend(TSmodel.reactionNetwork.reactionList.getReactionSubSet(cofactor='nadph_c').getReactionNameList(level=1))
    fluxNames.extend(TSmodel.reactionNetwork.reactionList.getReactionSubSet(cofactor='nadh_c').getReactionNameList(level=1))
      
    fluxNames = list(set(fluxNames)) # Eliminate redundancies 
    
    return fluxNames



def getTSReacNet(BASEfilename,FLUXESfilename,REACTIONSfilename,MSfilename,MSSTDfilename,FEEDfilename):
    """
    Produces the TS reaction Network given the necessary input files
    """

    # Load initial SBML file
    reactionNetwork = ReactionNetworks.TSReactionNetwork(BASEfilename)
    
    # Add Measured fluxes
    reactionNetwork.loadFluxBounds(FLUXESfilename)
    # Add carbon transitions
    reactionNetwork.addTransitions(REACTIONSfilename)
    # Add measured labeling information
    reactionNetwork.addLabeling(MSfilename,'LCMSLabelData',MSSTDfilename,minSTD=0.001)
    # Add feed labeling information
    reactionNetwork.addFeed(FEEDfilename)
 

    # Limit fluxes to 500       
    reactionNetwork.capFluxBounds(500)
    
    return reactionNetwork  


def changed(old,new,delta,relative=True):
    """ Tests if a number changed significantly
    -) delta is the maximum change allowed    
    -) relative decides if the delta given indicates relative changes (if True) or absolute change (if False)
    """
    delta = abs(delta)
    epsilon = 1.0
    
    if old > epsilon:
        if relative:
            notChanged = (new <= (1+delta)*old) and (new >= (1-delta)*old) 
        else:
            notChanged = (new <= old+delta) and (new >= old-delta)             
    elif old < -epsilon:
        if relative:
            notChanged = (new >= (1+delta)*old) and (new <= (1-delta)*old) 
        else:
            notChanged = (new >= old-delta) and (new <= old+delta)             
    else:
        notChanged = (new >= old-epsilon) and (new <= epsilon+old)         
    return not notChanged


# Auxiliary function used to change cofactors, e.g.: nadh --> nad   
def changeMetabolite(reactionNetworkOld,rxnName,changes):
    """
    Produces the changes given in the changes dictionary for reaction rxnName in the reactionNetwork given.
    For example:
    
    changes = {'nadp_c':'nad_c','nadph_c':'nadh_c'}    

    substitutes NADP by NAD, and NADPH by NADH.    
    """    
    #import core, copy
    
    # Deep copy of original network, which will be untouched
    reactionNetwork = copy.deepcopy(reactionNetworkOld)
    # Getting needed dictionaries
    reactionDict = reactionNetwork.reactionList.getReactionDictionary()
    metDict      = reactionNetwork.metList.metDict
    
    # Getting reaction to be changed
    rxn = reactionDict[rxnName]
    
    # Getting products and reactants for desired reaction
    products  = rxn.products
    reactants = rxn.reactants
    
    # Creating new products and reactants list to be added
    products2  = []
    reactants2 =[]

    for product in products:
        if product.name in changes:
            newmet = metDict[changes[product.name]]
            newProduct = core.Product(newmet, product.stoichiometry)
            products2.append(newProduct)
        else:
            products2.append(product)
    for reactant in reactants:
        if reactant.name in changes:
            newmet = metDict[changes[reactant.name]]
            newReactant = core.Reactant(newmet, reactant.stoichiometry)
            reactants2.append(newReactant)
        else:
            reactants2.append(reactant)
            
    # Replacing old products and reactants lists
    rxn.products  = products2
    rxn.reactants = reactants2
    
    return reactionNetwork


def compareGraph(reactionList1, reactionList2, reactionList3, maxPlot, titleFig='', save=''):
    """ 
    Auxiliary function to create paper figure  
    """
    font = {'size'   : 20}
    matplotlib.rc('font', **font)
    
    # Obtain data
    reactionDict  = reactionList2.getReactionDictionary()
    reactionDict2 = reactionList3.getReactionDictionary()

    names  = []
    lower  = []
    upper  = []
    total  = []
    lower2 = []
    upper2 = []
    total2 = []
    lower3 = []
    upper3 = []
    total3 = []
    for reaction in reactionList1:           
        if len(names) < maxPlot:
            names.append(reaction.name.replace('EX_','',1).replace('_e_','',1).replace('(e)','',1).replace('_','-'))
            
            lower.append(reaction.flux.net.lo)
            upper.append(reaction.flux.net.hi-reaction.flux.net.lo)
            total.append(reaction.flux.net.hi)
 
            reaction2 = reactionDict[reaction.name]
            lower2.append(reaction2.flux.net.lo)
            upper2.append(reaction2.flux.net.hi-reaction2.flux.net.lo)
            total2.append(reaction2.flux.net.hi)
            
            reaction3 = reactionDict2[reaction.name]
            lower3.append(reaction3.flux.net.lo)
            upper3.append(reaction3.flux.net.hi-reaction3.flux.net.lo)
            total3.append(reaction3.flux.net.hi)
            
    # bar data
    N = len(names)
    ind = numpy.arange(N)+0.5    # the x locations for the groups
    width = 0.25              # the width of the bars
        
    fig = figure(figsize=(16,8))
    axes = fig.add_subplot(111) 

    axes.axis([0,N+width,min(1.1*min(min(lower3),0),min(total3)),max(total3)*1.1])
         
    p1 = plt.bar(ind-0.29, lower, width, color='b', alpha=0.1)
    p2 = plt.bar(ind-0.29, upper, width, color='b', bottom=lower,alpha=0.8)
    p3 = plt.bar(ind+0.0, lower2, width, color='k', alpha=0.1)
    p4 = plt.bar(ind+0.0, upper2, width, color='k', bottom=lower2,alpha=0.8)
    p5 = plt.bar(ind+0.29, lower3, width, color='r', alpha=0.1)
    p6 = plt.bar(ind+0.29, upper3, width, color='r', bottom=lower3,alpha=0.8)
    
    plt.ylabel('Flux (mmol/gdw/h)')
    plt.title('Minimum and maximum values for exchange fluxes')
    plt.xticks(ind+width/2.,tuple(names),rotation='vertical')
    plt.legend( (p2[0], p1[0], p4[0], p3[0], p6[0], p5[0]),
               ('Max. stoich. + LFTC + $^{13}$C data', 'Min. stoich. + LFTC + $^{13}$C data','Max. stoich. + lim. flux to core (LFTC)', 'Min. stoich. + lim. flux to core (LFTC)','Max. only stoich.', 'Min. only stoich.'),
               prop={'size':18})

    title(titleFig)
    # Change x tick label size
    for label in axes.get_xticklabels():
        label.set_fontsize(20)
    # Inclination for xlabels
    fig.autofmt_xdate()   
    
    # Save Figures
    if save:
        savefig(save)    


def getFigure(number,saved=True):
    """
    Gets figure given by number
    -) saved  argument determines wether to used saved results (faster) or not    
    """
    
    if number == 4:
        "ELVA figure"
        
        strain ='wt5h'
        
        ### ELVA with new reaction set
        TSresult, TSmodel = getTSresults(strain,saved=saved)
        resultELVA  = TSmodel.ELVA(TSresult)  
         
        ### ELVA with old reaction set
        TSresultOld, TSmodelOld = getTSresults(strain,oldVersion=True)
        resultELVAOld  = TSmodelOld.ELVA(TSresultOld)            
    
        # Plots
        resultELVA.plotExpvsCompLabelxvsy(titleFig=strain,outputFileName="ELVAComparison"+strain+".txt",
                                               save="ELVA"+strain+".eps")            
                                               
        resultELVAOld.plotExpvsCompLabelxvsy(titleFig=strain,outputFileName="ELVAComparison"+strain+"Old.txt",
                                      save="ELVA"+strain+"Old.eps")     
        
    elif number == 5:
        "Relative flux constraint effect figure"        
        
        subsystems   = ['S_Pentose_Phosphate_Pathway']  

        test = testConstPower('refreshSavedFiles')
        TSreactionList, constFVAreactionList, FVAreactionList = test.getResults(saved=False)
        
        fluxNamesPPP  = ['G6PDH2r','GND','PGL','RPE','TK1','TA2','EDA','EDD','TK2','TA1','TK3','RPI'] 
        fluxNamesGLYC = ['GLCpts','PGI','PFK','FBA','TPI','F6PA','DHAPT','GAPD','PGK','PGM','ENO','PYK','PDH']        
        
        ##### PPP plot
        TSreactionListPPP       = TSreactionList.getReactionSubSet(rxnNames=fluxNamesPPP)
        constFVAreactionListPPP = constFVAreactionList.getReactionSubSet(rxnNames=fluxNamesPPP)
        FVAreactionListPPP      = FVAreactionList.getReactionSubSet(rxnNames=fluxNamesPPP)     
        
        TSreactionListPPP.name       = TSreactionList.name 
        constFVAreactionListPPP.name = constFVAreactionList.name
        FVAreactionListPPP.name      = FVAreactionList.name

        utils.fluxPlotGraph([TSreactionListPPP, constFVAreactionListPPP, FVAreactionListPPP], subsystems,
                        titleFig='', save='Comparison2.png', plotType=2, axesRanges=[-20,40], legendFontSize=14)
                        
        # Quantify degree of constraining
        fluxDictTSPPP  = TSreactionListPPP.getFluxDictionary(rangeComp='all')
        fluxDictFVAPPP = FVAreactionListPPP.getFluxDictionary(rangeComp='all')      
        QfluxPPP=[]
        for name in fluxDictTSPPP:
            DTS = abs(fluxDictTSPPP[name].net.hi  - fluxDictTSPPP[name].net.lo)
            DFVA= abs(fluxDictFVAPPP[name].net.hi - fluxDictFVAPPP[name].net.lo)
            if DFVA !=0 :
                QfluxPPP.append(DTS/DFVA)
        QPPP = float(sum(QfluxPPP))/len(QfluxPPP)
        print 'TS is '+str(QPPP)+' % of FVA'
      
        ##### Glycolysis plot                                              # Duplicated code, but ok for the moment being
        TSreactionListGLYC       = TSreactionList.getReactionSubSet(rxnNames=fluxNamesGLYC)
        constFVAreactionListGLYC = constFVAreactionList.getReactionSubSet(rxnNames=fluxNamesGLYC)
        FVAreactionListGLYC      = FVAreactionList.getReactionSubSet(rxnNames=fluxNamesGLYC)     
        
        TSreactionListGLYC.name       = TSreactionList.name 
        constFVAreactionListGLYC.name = constFVAreactionList.name
        FVAreactionListGLYC.name      = FVAreactionList.name

        utils.fluxPlotGraph([TSreactionListGLYC, constFVAreactionListGLYC, FVAreactionListGLYC], subsystems,
                        titleFig='', save='Comparison2.png', plotType=2, axesRanges=[-30,50], legendFontSize=14)
                       
        # Quantify degree of constraining
        fluxDictTSGLYC  = TSreactionListGLYC.getFluxDictionary(rangeComp='all')
        fluxDictFVAGLYC = FVAreactionListGLYC.getFluxDictionary(rangeComp='all')      
        QfluxGLYC=[]
        for name in fluxDictTSGLYC:
            DTS = abs(fluxDictTSGLYC[name].net.hi  - fluxDictTSGLYC[name].net.lo)
            DFVA= abs(fluxDictFVAGLYC[name].net.hi - fluxDictFVAGLYC[name].net.lo)
            if DFVA !=0 :
                QfluxGLYC.append(DTS/DFVA)
        QGLYC = float(sum(QfluxGLYC))/len(QfluxGLYC)
        print 'TS is '+str(QGLYC)+' % of FVA'
    
    elif number == 6:
        "Cofactor balances figure"
        # Get results
        TSresults = {}
        TSmodels  = {}
        TSresults['wt5h']  , TSmodels['wt5h']   = getTSresults('wt5h',saved=saved)
        TSresults['pgi21h'], TSmodels['pgi21h'] = getTSresults('pgi21h',saved=saved)

        # Plot fluxes
        TSresults['wt5h'].plotMetaboliteFlux('nadph_c',titleFig='wt5h',minFluxGroup = 0.3, maxFluxGroup=1.5,save='cofNADPHwt5h.png')
        TSresults['pgi21h'].plotMetaboliteFlux('nadph_c',titleFig='pgi21h',minFluxGroup = 0.3, maxFluxGroup=1.5,save='cofNADPHpgi21h.png')           
        TSresults['wt5h'].plotMetaboliteFlux('nadh_c', titleFig='wt5h',minFluxGroup = 0.3, maxFluxGroup=0.6,save='cofNADHwt5h.png')
        TSresults['pgi21h'].plotMetaboliteFlux('nadh_c',titleFig='pgi21h',minFluxGroup = 0.3, maxFluxGroup=0.6,save='cofNADHpgi21h.png')

    elif number == 7:    
        "Robustness with respect to measurement error in labeling profile"

        # Get model and results
        TSmodelRob = getTSmodel('wt5h')        
        resultsRob = TSmodelRob.findFluxesStds(Nrep=30,Nrand=30)        
        
        # Plotting results
        subsystems    = ['S_Pentose_Phosphate_Pathway']
        fluxNamesPlot = ['G6PDH2r', 'GND', 'PGL', 'RPE', 'TK1', 'TA2', 'EDA', 'EDD', 'TK2', 'TA1', 'TK3', 'RPI']

        TSreactionListRob       = resultsRob.reactionNetwork.reactionList.getReactionSubSet(rxnNames=fluxNamesPlot)
        TSreactionListRob.name = 'Stds'

        utils.fluxPlotGraph([TSreactionListRob], subsystems,
                        titleFig='', save='LabelingRobust.png', plotType=2, axesRanges=[-2,4])

    elif number == 8:    
        "Robustness with respect to genome annotations errors"

        fluxNamesTCA = ['CS', 'ACONT', 'ICDHyr', 'AKGDH', 'SUCOAS', 'SUCD1i', 'FUM', 'MDH']
        subsystems  = ['S_Pentose_Phosphate_Pathway']
        changes = {'nadp_c':'nad_c','nadph_c':'nadh_c'}
        
        #### FBA figure

        # Get models
        #cofModelsFBA1= FluxModels.FBAModel('EciJR904TKswt5hTS.xml')
        #cofModelsFBA2= FluxModels.FBAModel('EciJR904TKswt5hTS_G6PDH2r.xml')
        strain        = 'wt5h'
        datadir       = basedir+strain+'/'
        cofModelsFBA1 = FluxModels.FBAModel(datadir + 'EciJR904TKs.xml')
        cofModelsFBA1.reactionNetwork.capFluxBounds(500)        
        
        reactionNetworkNew = changeMetabolite(cofModelsFBA1.reactionNetwork,'G6PDH2r',changes)
        reactionNetworkNew.write('EciJR904TKswt5hTS_G6PDH2r.sbml') 
        cofModelsFBA2= FluxModels.FBAModel('EciJR904TKswt5hTS_G6PDH2r.sbml')

        # Get results
        cofResultsFBA1 = cofModelsFBA1.findFluxes()
        cofResultsFBA2 = cofModelsFBA2.findFluxes()

        # Get reactionLists
        reactionListFBA1    = cofResultsFBA1.reactionNetwork.reactionList.getReactionSubSet(rxnNames=fluxNamesTCA)
        reactionListFBA2    = cofResultsFBA2.reactionNetwork.reactionList.getReactionSubSet(rxnNames=fluxNamesTCA)
        reactionListFBA1.name = 'Original'
        reactionListFBA2.name = 'Cofactor changed'

        # Plot fluxes in reactionLists
        utils.fluxPlotGraph([reactionListFBA1, reactionListFBA2], subsystems,
                        titleFig='FBA', plotType=1, axesRanges=[-8,8])    

        #### 2S-13C MFA figure
        # Get models
        TSmodel = getTSmodel('wt5h')
        TSmodel.reactionNetwork.write('Base2.sbml')
        reactionNetworkNew = changeMetabolite(TSmodel.reactionNetwork,'G6PDH2r',changes)
        reactionNetworkNew.write('Base2_G6PDH2r.sbml')  
        
        cofModelsTS1= FluxModels.TwoSC13Model('Base2.sbml')
        cofModelsTS2= FluxModels.TwoSC13Model('Base2_G6PDH2r.sbml')

        # Get results
        cofResultsTS1 = cofModelsTS1.findFluxesRanges(Nrep=30,fluxNames=fluxNamesTCA, procString='proc') 
        cofResultsTS2 = cofModelsTS2.findFluxesRanges(Nrep=30,fluxNames=fluxNamesTCA, procString='proc') 

        # Get reactionLists
        reactionListTS1    = cofResultsTS1.reactionNetwork.reactionList.getReactionSubSet(rxnNames=fluxNamesTCA)
        reactionListTS2    = cofResultsTS2.reactionNetwork.reactionList.getReactionSubSet(rxnNames=fluxNamesTCA)
        reactionListTS1.name = 'Original'
        reactionListTS2.name = 'Cofactor changed'

        # Plot fluxes in reactionLists
        utils.fluxPlotGraph([reactionListTS1, reactionListTS2], subsystems,
                        titleFig='2S 13C MFA', plotType=1, axesRanges=[-10,14])

    elif number == 9:    
        "Extracellular metabolite prediction figure"
        
        # Get results
        test = testConstPower('refreshSavedFiles')
        TSreactionList, constFVAreactionList, FVAreactionList = test.getResults(saved=saved)        
        
        # Compare fluxes
        allReacNames  = TSreactionList.getReactionSubSet('exchange').getReactionNameList(1)
        
        excReactionList = TSreactionList.getReactionSubSet(rxnNames=allReacNames)     
        excReactionList.sort(function='by absolute low flux')

        compareGraph(excReactionList, constFVAreactionList, FVAreactionList, 12, save='Figure9')       

    elif number == 10:
        "Flux comparison with COBRA predictions figure"
        strain = 'pyk5h'  
        
        # Get results
        TSresult, TSmodel = getTSresults(strain,saved=saved)
        test = testPredictions('refreshSavedFiles')
        #resultsV1,resultsV2,TSMOMA,TS13CMOMA,results13CMOMAFixed,resultsMOMAFixed = test.getResults(saved=saved)  
        resultsVGP,resultsVGF,resultsVAP,resultsVAF = test.getFVAresults(strain)
        TSMOMA,TS13CMOMA,TSROOM,TS13CROOM           = test.getFluxPredResults(strain)        
        
        # Compare fluxes
        subsystems  = ['S_Pentose_Phosphate_Pathway']
        rxnNames = ['G6PDH2r','GND','PGL','RPE','TK1','TA2','EDA','EDD','TK2','TA1','TK3','RPI'] 

        TSreactionList      = TSresult.reactionNetwork.reactionList.getReactionSubSet(rxnNames=rxnNames)
        FVAreactionList     = resultsVGF.reactionNetwork.reactionList.getReactionSubSet(rxnNames=rxnNames)
        ATPreactionList     = resultsVAF.reactionNetwork.reactionList.getReactionSubSet(rxnNames=rxnNames)
        C13MOMAreactionList = TS13CMOMA.reactionNetwork.reactionList.getReactionSubSet(rxnNames=rxnNames)

        TSreactionList.name      = '2S-$^{13}$C MFA'    
        FVAreactionList.name     = 'FBA max. growth'
        ATPreactionList.name     = 'FBA max. ATP'
        C13MOMAreactionList.name = '$^{13}$C MOMA'

        utils.fluxPlotGraph([TSreactionList, FVAreactionList, ATPreactionList, C13MOMAreactionList], subsystems,
                        titleFig=strain, save='compFull'+strain+'FinalA.png', plotType=1, axesRanges=[-8,8])

    elif number == 11:
        "Metabolite labeling prediction figure"
        strain = 'pyk5h'          
        
        # Get results
        TSresult, TSmodel = getTSresults(strain,saved=saved)
        test = testPredictions('refreshSavedFiles')
        TSMOMA,TS13CMOMA,TSROOM,TS13CROOM           = test.getFluxPredResults(strain)           
        results13CMOMAFixed,resultsMOMAFixed        = test.getLabelPredResults(TS13CMOMA,TSMOMA)

        # Compare lagbeling predictions
        expData = [results13CMOMAFixed.fragDict,resultsMOMAFixed.fragDict]
        fitData = [results13CMOMAFixed.EMUlabel,resultsMOMAFixed.EMUlabel]
        labels  = ['$^{13}$C MOMA, OF='+'{0:4.1f}'.format(results13CMOMAFixed.OF),'MOMA, OF='+'{0:4.1f}'.format(resultsMOMAFixed.OF)]
        utils.labelingGraph(expData,fitData,titleFig='',labels=labels,save='labelFinal.png'
                        ,colors=['g','y'],formats=['*','v'],markersizes=[13,7])

    else:
        raise Exception("Figure "+str(number)+" not available")


def getSuppFigure(number,saved=True):
    """Gets supplementary figure given by number
        saved  argument determines wether to used saved results (faster) or not    
    """
    
    if number in [1,2,3]:
        "Labeling fits"
        
        for strain in strains:   
            TSresult, TSmodel = getTSresults(strain,saved=saved)            
            TSresult.plotExpvsCompLabelFragment(titleFig=strain,save="fit"+strain+".eps")
    elif number in range(4,13):
        "Fluxmaps"
        
        for strain in strains:
            TSresult, TSmodel = getTSresults(strain,saved=saved)            
            TSresult.reactionNetwork.name = 'E. coli '+strain
            TSresult.drawFluxes('EcoliTS'+strain+'.svg',svgInFileName='TOYAexpOLD2.svg',oldVersion=True)
    elif number == 13:
        "Extracellular fluxes prediction for pyk 6 hrs"
        # Get results
        test = testConstPower('refreshSavedFiles')
        TSreactionList, constFVAreactionList, FVAreactionList = test.getResults(strain='pyk6h',saved=saved)  
        
        # Compare fluxes
        allReacNames  = TSreactionList.getReactionSubSet('exchange').getReactionNameList(1)
        
        excReactionList = TSreactionList.getReactionSubSet(rxnNames=allReacNames)     
        excReactionList.sort(function='by absolute low flux')

        compareGraph(excReactionList, constFVAreactionList, FVAreactionList, 12, save='FigureS13') 
        
    elif number == 14:
        "Extracellular fluxes prediction for pgi 21 hrs"
        # Get results
        test = testConstPower('refreshSavedFiles')
        TSreactionList, constFVAreactionList, FVAreactionList = test.getResults(strain='pgi21h',saved=saved)  
        
        # Compare fluxes
        allReacNames  = TSreactionList.getReactionSubSet('exchange').getReactionNameList(1)
        
        excReactionList = TSreactionList.getReactionSubSet(rxnNames=allReacNames)     
        excReactionList.sort(function='by absolute low flux')

        compareGraph(excReactionList, constFVAreactionList, FVAreactionList, 12, save='FigureS14') 
        
    elif number == 15:
        "Robustness with respect to stoichiometric errors for GLUDy and ICDHyr"

        fluxNamesTCA = ['CS', 'ACONT', 'ICDHyr', 'AKGDH', 'SUCOAS', 'SUCD1i', 'FUM', 'MDH']
        subsystems  = ['S_Pentose_Phosphate_Pathway']
        changes = {'nadp_c':'nad_c','nadph_c':'nadh_c'}

        #### FBA figure
        # Create models
        strain        = 'wt5h'
        datadir       = basedir+strain+'/'
        cofModelsFBA1 = FluxModels.FBAModel(datadir + 'EciJR904TKs.xml')
        cofModelsFBA1.reactionNetwork.capFluxBounds(500)    
        
        reactionNetworkGLUDy = changeMetabolite(cofModelsFBA1.reactionNetwork,'GLUDy',changes)
        reactionNetworkGLUDy.write('Base_GLUDy.sbml')  
        cofModelsFBA3= FluxModels.FBAModel('Base_GLUDy.sbml')
        reactionNetworkICDHyr = changeMetabolite(cofModelsFBA1.reactionNetwork,'ICDHyr',changes)
        reactionNetworkICDHyr.write('Base_ICDHyr.sbml')         
        cofModelsFBA4= FluxModels.FBAModel('Base_ICDHyr.sbml')    

        # Get results
        cofResultsFBA1 = cofModelsFBA1.findFluxes()
        cofResultsFBA3 = cofModelsFBA3.findFluxes()
        cofResultsFBA4 = cofModelsFBA4.findFluxes()

        # Get reaction lists
        reactionListFBA1    = cofResultsFBA1.reactionNetwork.reactionList.getReactionSubSet(rxnNames=fluxNamesTCA)
        reactionListFBA3    = cofResultsFBA3.reactionNetwork.reactionList.getReactionSubSet(rxnNames=fluxNamesTCA)
        reactionListFBA4    = cofResultsFBA4.reactionNetwork.reactionList.getReactionSubSet(rxnNames=fluxNamesTCA)
        reactionListFBA1.name = 'Original'
        reactionListFBA3.name = 'Cofactor changed'
        reactionListFBA4.name = 'Cofactor changed'

        # Make flux plots
        utils.fluxPlotGraph([reactionListFBA1, reactionListFBA3], subsystems,
                        titleFig='FBA', plotType=1, axesRanges=[-8,8])

        utils.fluxPlotGraph([reactionListFBA1, reactionListFBA4], subsystems,
                        titleFig='FBA', plotType=1, axesRanges=[-8,8])
                        

        #### 2S-13C MFA figure
        TSmodel = getTSmodel('wt5h')
        fluxNames = getFluxNames(TSmodel)
        TSmodel.reactionNetwork.write('Base.sbml')
        reactionNetworkGLUDy = changeMetabolite(TSmodel.reactionNetwork,'GLUDy',changes)
        reactionNetworkGLUDy.write('Base_GLUDy.sbml')
        reactionNetworkICDHyr = changeMetabolite(TSmodel.reactionNetwork,'ICDHyr',changes)
        reactionNetworkICDHyr.write('Base_ICDHyr.sbml')

        # Create models
        cofModelsTS1= FluxModels.TwoSC13Model('Base.sbml')
        cofModelsTS3= FluxModels.TwoSC13Model('Base_GLUDy.sbml')
        cofModelsTS4= FluxModels.TwoSC13Model('Base_ICDHyr.sbml')

        # Get results
        cofResultsTS1 = cofModelsTS1.findFluxesRanges(Nrep=30, fluxNames=fluxNamesTCA, procString='proc') 
        cofResultsTS3 = cofModelsTS3.findFluxesRanges(Nrep=30, fluxNames=fluxNamesTCA, procString='proc') 
        cofResultsTS4 = cofModelsTS4.findFluxesRanges(Nrep=30, fluxNames=fluxNamesTCA, procString='proc') 

        # Get reaction lists
        reactionListTS1    = cofResultsTS1.reactionNetwork.reactionList.getReactionSubSet(rxnNames=fluxNamesTCA)
        reactionListTS3    = cofResultsTS3.reactionNetwork.reactionList.getReactionSubSet(rxnNames=fluxNamesTCA)
        reactionListTS4    = cofResultsTS4.reactionNetwork.reactionList.getReactionSubSet(rxnNames=fluxNamesTCA)

        reactionListTS1.name = 'Original'
        reactionListTS3.name = 'Cofactor changed'
        reactionListTS4.name = 'Cofactor changed'

        # Make flux plots
        utils.fluxPlotGraph([reactionListTS1, reactionListTS3], subsystems,
                        titleFig='2S 13C MFA', plotType=1, axesRanges=[-10,16])                        

        utils.fluxPlotGraph([reactionListTS1, reactionListTS4], subsystems,
                        titleFig='2S 13C MFA', plotType=1, axesRanges=[-10,16])        

    elif number == 16:
        """
        13C MFA fluxes for wt5h
        From 'Final Figures' ipython notebook        
        """              

        strain = 'wt5h'
        results, model = get13CMFAresults(strain)
        
        # Draw map
        results.drawFluxes('ecoli'+strain+'13C.svg',svgInFileName='TOYA13CMFA.svg',norm='GLCpts',oldVersion=True)
        
        # Show map
        SVG(filename='ecoli'+strain+'13C.svg') 

    elif number == 17:
        "Comparison of full predictions"
        
        strainsFull = ['pyk5h','pgi16h']        
        test = testPredictions('refreshSavedFiles')
        
        for strain in strainsFull:
            test.getFullCompFig(strain,saved=saved)            

    elif number == 18:
        "Comparison of partial predictions"
        
        strainsPart = ['wt6h','pyk6h','pgi21h']       
        test = testPredictions('refreshSavedFiles')
        
        for strain in strainsPart:
            test.getPartCompFig(strain,saved=saved)            

    elif number == 19:
        """
        OF vs N Figure
        This one takes forever!!!
        """  
        
        # Get fluxNames
        strain = 'wt5h'
        TSmodel   = getTSmodel(strain)    
        fluxNames = getFluxNames(TSmodel)
        
        TSresults = {}
        
        Ns    = [1,2,5,10,20,30,40,50,100]
        Nrand = 10
        OFsAv= {}
        
        # Get results for different number of replicates
        for Nreplicates in Ns:
            OFs = []
            for index in range(1,Nrand+1):
                TSmodel   = getTSmodel(strain)    
                TSresults = TSmodel.findFluxesRanges(Nrep=Nreplicates,fluxNames=fluxNames, procString='proc') 
                OFs.append(TSresults.OF)
            OFsAv[Nreplicates] = sum(OFs)/float(len(OFs))

        # Assemble graph
        x = []
        y = []
        for reps in Ns:
            x.append(reps)
            y.append(OFsAv[reps])       

        plot(x,y,'o-')
        xlabel('Number of runs')
        ylabel('Average objective Function (OF)')        

        print x
        print y

    else:
        raise Exception("Figure "+str(number)+" not available")


# Auxiliary functions used in the optimization files test to eliminate expected variability such as the random seed        
def filterFiles(fileTuple):
    output       = fileTuple
    name, string = fileTuple
    
    if name == 'condensation_rxns_eqns.txt' or name == 'condensation_rxns_names.txt':
        string = '\n'.join(sorted(string.split('\n')))
        output = (name, string)
    elif name == 'randseed.txt':  # eliminate randseeds
        output =(name, 'randseed\n')
    elif 'TSresults' in name: # eliminate xml file which is not really part of input, just stored in the same place
        output = None
    elif 'nfs' in name: # eliminate damn nfs handles
        output = None
        
    return output


def readFiles(dir_):  # TODO: innline this one
    """Reads all files in a directory"""
    output = []
    names = os.listdir(dir_)
    for name in names:
        file_  = open(dir_+'/'+name)
        string = ''.join(file_.readlines())
        output.append((name,string))
     
    #output = [( name ,''.join(open(dir_+'/'+name).readlines()) ) for name in os.listdir(dir_)]   --> inlined version 
    
    return output
