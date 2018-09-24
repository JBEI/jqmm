# Part of the JBEI Quantitative Metabolic Modeling Library (JQMM)
# Copyright (c) 2016, The Regents of the University of California.
# For licensing details see "license.txt" and "legal.txt".
"""
This module contains all the methods used to treat labeling data (i.e. Mass Distribution Vectors)
"""
from __future__ import print_function
from __future__ import division

from builtins import str
from builtins import range
from builtins import object
from past.utils import old_div
import re
import numpy, random
import DB
import utilities as utils


##################################
# ALL TYPES OF LABELING DATA
##################################


# Labeling data classes
class labelingData(object):
    "General class for storing labeling data data"

    def __init__(self,lines = 'default',dataTup='default',minSTD=0):  # TODO: This init can be improved by using class methods      
        if lines !='default':
            # Unpacking lines input = tuple of MS line and MS stds line
            if lines.__class__.__name__ == 'tuple':
                line,STDline = lines
            else:
                line    = lines
                STDline = ''
            # Eliminate header from lines
            line      = line.split(':')[-1]
            
            if STDline:
                STDline = STDline.split(':')[-1]   
            # Components from line
            self.fragment = self.getFragmentName(line)    
            self.inFit = True   # This decides wether fragment info is included in the fit            
            self.abbrev   = self.getAbbrev(self.fragment)

            comps    = line.split('\t')
            
            # MDV (and STD if present) values
            values = '\t '.join(comps[2:])
            MDVSTDs = re.findall('([\.\d\-eE]+)\(([\.\d\-eE]+)\)',values)      # input can be in the form mdv(std): 0.234(0.030) or the form mdv: 0.234
            MDVs    = re.findall('([\.\d\-eE]+)',values)

            
            MDVlist = []
            STDlist = []
            if MDVSTDs:
                for tup in MDVSTDs:
                    mdv, STD = tup
                    #mdv = 0 if mdv == '-' else mdv
                    if not mdv == '-':
                        MDVlist.append(float(mdv))
                        STDlist.append(float(STD))
            else:
                if MDVs:
                    for mdv in MDVs:
                    #    mdv = 0 if mdv == '-' else mdv
                        if not mdv == '-':
                            MDVlist.append(float(mdv))
                else:
                    raise Exception('Line input cannot be parsed: '+line)
                    
                    
            # STD values
            if STDline:
                STDfragment = self.getAbbrev(STDline.split('\t')[0].strip())
                if STDfragment != self.fragment:
                    raise Exception('STD line fragment name does not correspond to MDV line fragment: '+STDfragment+' vs '+self.fragment)
                else:
                    STDs    = re.findall('([\.\d]+)','\t'.join(STDline.split('\t')[2:]))
                    STDlist = []
                    for STD in STDs:
                        STDlist.append(float(STD))
            else:
                if not STDlist:                    # Fill up the standard deviation with default result if nothing else is available
                    STDlist = list(0.03*numpy.ones(shape=(len(MDVlist))))
            # Capping minimum value of STDs
            STDlistMod = []
            for value in STDlist:
                value = max(value,minSTD)
                STDlistMod.append(value)
                    

            self.mdv = numpy.array(MDVlist)
            
            STDlist = STDlist[0:len(MDVlist)]
            self.std = numpy.array(STDlistMod)

        elif dataTup !='default':                    # This part of labeling data initialization should be updated
            fragment, mdv, STD = dataTup
            
            self.fragment = fragment
            self.mdv = mdv
            if STD:
                self.std = STD
            else:
                self.std = numpy.array(list(0.03*numpy.ones(shape=(len(self.mdv)))))
            self.abbrev = self.getAbbrev(fragment)
                    
        else:
            raise Exception('Incorrect input: must input lines or MDV info')
    
        self.line = self.getLine()


    def getMval(self,fragment):
        [aa,mVal] = fragment.split('M-')
        return mVal


    def getLine(self,maxVals=13):
        "Converts MDVs and STDs into a single string line"
        base     = '-\t'
        fragment = self.fragment
        mdv      = self.mdv
        try:
            std    = self.std
            stdStr = []
            for i in range(len(std)):
                stdStr.append('{0:2.3f}'.format(std[i]))
            for i in range(len(std),len(mdv)):
                stdStr.append('0')                
        except AttributeError:
            # No std provided
            for i in range(len(mdv)):
                stdStr.append('--')
        
        MDVlist = []
        
        for i in range(len(mdv)):
            MDVlist.append('{0:2.3f}'.format(mdv[i])+'('+stdStr[i]+')')
        mVal = self.getMval(fragment)
        name = self.getName()               # This should probably be done through the fragment class
        line = name+'\t'+'M-'+mVal+'\t'+'\t'.join(MDVlist)+'\t'+base*(maxVals-len(MDVlist))
    
        return line       
                
    
    def randomize(self):
        mdv    = self.mdv
        STD    = self.std
        newMDV = []
        for i in range(len(mdv)):
            val = mdv[i]
            std = STD[i]
            newVal = numpy.maximum(val+(2*random.random()-1)*std,0)
            newMDV.append(newVal)
        self.mdv = numpy.array(newMDV)
        
        # Update line        
        self.line = self.getLine()

        
    def normalize(self):
        total = sum(self.mdv)
        if total != 0:
            self.mdv   = old_div(self.mdv, total)
         
        # Update line        
        self.line = self.getLine()
    


class LCMSLabelData(labelingData):
    "Class for storing LCMS data"

    def getName(self):
        return self.abbrev


    def getMval(self,fragment):
        return '0'


    def getAbbrev(self,fragment,standard = False):
        "Creates abbreviation for fragment"   
        abbrev = fragment.strip()      # e.g. 3pg 
        
        if standard:     # Generate standarized alias
            abbrev, warning = getStdAlias(abbrev)  
            if not abbrev:
                print(warning)
        
        return abbrev    


    def getFragmentName(self,line):
        "Creates abbreviation for fragment"   
        comps  = line.split('\t')
        
        part1 = comps[0].strip()
        part2 = comps[1].strip().replace('M-0','',1).replace('M-','',1).strip()

        fragmentName = part1+part2
        
        return fragmentName



class CEMSLabelData(LCMSLabelData):
    "Class for storing CEMS data"
    def nothing(self):
        "I can't have an inherited class without adding something, it seems."
        pass
        


class GCMSLabelData(labelingData):
    "Class for storing GCMS data"

    def getAbbrev(self,fragment,standard = False):
        "Creates abbreviation for fragment"  
        abbrev = fragment.replace('M-','',1)      # e.g. Ala159 
        number = re.search('(\d+)',abbrev).group(1)
        abbrev = abbrev.replace(number,'',1)        
        
        if standard:     # Generate standarized alias
            abbrev, warning = getStdAlias(abbrev)  
            if not abbrev:
                print(warning)
                
        return abbrev+number


    def getName(self):
        name = re.sub('(\d+)','',self.abbrev)                # This should probably be done through the fragment class
        return name


    def getFragmentName(self,line):
        "Creates abbreviation for fragment"   
        comps  = line.split('\t')
        fragmentName = comps[0].strip()+comps[1]        
        return fragmentName


    def getMval(self,fragment):
        [aa,mVal] = fragment.split('M-')
        return mVal


# Classes for fragments (GCMS,LCMS....)     

class fragment(object):
    """ A base class for all fragments"""

    def __init__(self,emu,name,abbrev,ncarbons,composition):
        """
        emu: emu name (e.g. alaL_2_3)
        name: fragment name (e.g. Alanine,159)
        abbrev: fragment abbreviation (e.g. Ala159)
        ncarbons: number of carbons (e.g. 2) 
        composition: fragment composition minus carbon backbone (e.g. C6H20NSi)
        """    
        
        self.emu      = emu
        self.name     = name
        self.abbrev   = abbrev
        self.ncarbons = ncarbons
        self.comp     = composition

    
    def __str__(self):
        string = self.name
        if hasattr(self, 'mdv'):
            string = string + ': ' + str(self.mdv)
        if hasattr(self, 'std'):
            string = string + ': ' + str(self.std)
        return string



class GCMSfragment(fragment):
    """ A class for GCMS fragments"""

    def __init__(self,emu,name,abbrev,ncarbons,composition,antonApproved,derWeight,derAtoms):
        """
        emu: emu name (e.g. alaL_2_3)
        name: fragment name (e.g. Alanine,159)
        abbrev: fragment abbreviation (e.g. Ala159)
        ncarbons: number of carbons (e.g. 2) 
        composition: fragment composition minus carbon backbone (e.g. C6H20NSi)
        antonApproved: approved by antoniewicz
        derWeight: weight of derivatized fragment
        derAtoms: composition of derivatized fragment
        """    
        
        fragment.__init__(self,emu,name,abbrev,ncarbons,composition)
        
        self.antonAp  = antonApproved
        self.derWeigth= derWeight
        self.derAtoms = derAtoms



class LCMSfragment(fragment):
    """ A class for LCMS fragments"""
    
    def nothing(self):
        "I can't have an inherited class without adding something, it seems."
        pass
    

class CEMSfragment(fragment):
    """ A class for CEMS fragments"""
    
    def nothing(self):
        "I can't have an inherited class without adding something, it seems."
        pass


        
class LabelinglString(object):
    """ An class for strings containing labeling info
        Example input:
                 lstring: U
                 ncarbons: 6
    """
    # Adding data input as part of the init                        
    def __init__(self, lstring, ncarbons):
        if lstring == 'U':
            # Fully labeled
            self.labelCarbons = ''.zfill(ncarbons).replace('0','1') 
        elif lstring == '1-C':
            # First Carbon labeled
            self.labelCarbons = '1'+''.zfill(ncarbons-1)
        elif lstring == '2-C':
            # Second Carbon labeled
            self.labelCarbons = '01'+''.zfill(ncarbons-2)
        elif lstring == '1,2-C':
            # First 2 carbons labeled
            self.labelCarbons = '11'+''.zfill(ncarbons-2)
        elif lstring == '2,3-C':
            # Second and third carbons labeled
            self.labelCarbons = '011'+''.zfill(ncarbons-3)
        elif lstring == 'UN':
            # Unlabeled
            self.labelCarbons = ''.zfill(ncarbons)
        else:
            raise Exception('Feed labeling code not recognized: ' + str(lstring))
        if self.labelCarbons == '':
            raise Exception('Labeling string %s with ncarbons %s resulted in empty self.labelCarbons' % (lstring, str(ncarbons)))


    def mdv(self, emu):    
        # Find list of atom positions involved in emu.

        # This turns a string like 'glc_D_e_3_4_5_6' into an array ['3', '4', '5', '6'] 
        elList = re.search('((_\d)+)', emu).group(0).split('_')[1:]
        mval = 0
        for el in elList:
            try:
                mval = mval + int(self.labelCarbons[int(el)-1])
            except IndexError:
                raise Exception('Trying to address index %s in self.labelCarbons "%s", based on emu "%s"' % (int(el)-1, self.labelCarbons, emu))                
        mdv = []
        for i in range(len(elList)+1):
            mdv.append(0)
        mdv[mval] = 1
        return numpy.array(mdv)



########
# Auxiliary functions
########

def getStdAlias(name):
    """
    Generates standard alias for a metabolite name, as from fragdictbasic    
    
    For example:
    
    getStdAlias('Asp-L') = ('Asp,'')
    """
        
    aliases        = genAliases(name)
    fragDictBasic_ = DB.fragDictBasic()
    presenceVector = [i in fragDictBasic_ for i in aliases]  # Vector with presence or absence for each alias
    if sum(presenceVector)>0:                                # checks if any of the aliases is in fragDictBasic keys
        newName = aliases[presenceVector.index(True)]        # Abbreviation name present in fragDictBasic
        warning = ''
    else:
        newName = ''
        warning = 'Warning: fragment '+str(name)+' not found in fragDictBasic:\n'+str(list(fragDictBasic_.keys()))

    return newName, warning
    

# Function which generates aliases to metabolite names
def genAliases(name):
    """ 
    Generates aliases for metabolite names, e.g.:
        
        val --> set(['Val-L', 'Val', 'val', 'val-L'])
    
    """
    name = name.replace('-L','').replace('-l','')
    output = []
    output.append(name)
    output.append(name.lower())
    output.append(name.lower()+'-L')
    output.append(name.lower()+'_L')
    output.append(name.capitalize())
    output.append(name.capitalize()+'-L')
    output.append(name.capitalize()+'_L')
    
    return output
    

def createGammaMatrix(distleng,string): # TODO: find a better place for this 
    """
    Function to create gamma matrix for the atoms found in the input string as per Wahl et al Biotechnol Bioeng 85: 259-268 (2004)
        -) distleng is the length of the distribution the matrix has to be multiplied by.
        Example: gammarray(10,'C3H2O')
    """
    # Natural isotope abundances
    abundDict = {
        'C' :numpy.array([[12,0.98918],[13,0.01082]]),
        'O' :numpy.array([[16,0.99756],[17,0.00038],[18,0.00206]]),
        'H' :numpy.array([[1,0.999844],[2,0.000156]]),
        'N' :numpy.array([[14,0.99634],[15,0.00366]]),
        'Si':numpy.array([[28,0.9222],[29,0.0469],[30,0.0309]]), 
        'S' :numpy.array([[32,0.95039],[33,0.00749],[34,0.04197],[36,0.00015]]), 
        'P' :numpy.array([[31,1]]) 
        }

    # Getting atom species and number
    atoms    = utils.getAtoms(string)

    # Initialization
    gammatot = numpy.eye(distleng)

    # Main loop
    for atom in atoms:
        symbol, number = atom
        if symbol in abundDict:
            for j in range(number):
                # Make gamma matrix
                distribution = abundDict[symbol][:,1]
                gamma = numpy.zeros((distleng+len(distribution)-1,distleng))
                for j in range(distleng):
                    gamma[j:j+len(distribution),j] = distribution.transpose()
                # Add to total gamma
                gammatot = numpy.matrix(gamma)*gammatot
                distleng = gamma.shape[0]
        else:
            raise Exception ('Atom '+symbol+' not in table, please add')

    return gammatot


