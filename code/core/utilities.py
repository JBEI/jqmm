# Part of the JBEI Quantitative Metabolic Modeling Library (JQMM)
# Copyright (c) 2016, The Regents of the University of California.
# For licensing details see "license.txt" and "legal.txt".
"""
Set of utilities used in the rest of the JQMM library
"""
from __future__ import print_function
from __future__ import division

from builtins import str
from builtins import range
from past.utils import old_div
import re, os
import core, DB
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from pylab import figure, plot, xlabel, ylabel, title, text, savefig
import numpy


def is_float(s):
    """
    Determines if a string is a float.

    For example:
    is_float('3.3')    = True
    is_float('potato') = False
    
    """
    try:
        float(s)
        return True
    except TypeError:
        return False
    except ValueError:
        return False
    except AttributeError:
        return False


def symsplit(item):
    """
    Splits symmetric molecule atom mapping expression into its component parts.
    Returns a list of atom mappings in the symmetric molecule.

    For example:

    symsplit('(abc;cba)') = ['abc','cba']
    """

    # Removing white space
    tmp=item.replace(' ','')
    # checking format      
    if (item[0] != '(') or (item[-1] != ')'):
        print("error in", item,"!\n")
    else:
        tmp=tmp[1:-1].split(';')
        # checking sizes to make sure they are all the same
        sizet=len(tmp[0])
        for t in tmp[1:]:
            if len(t) != sizet:
                print("Sizes of the symmetric molecules are not the same!\n")
                print(tmp)
    return tmp


def indexMap(inds,prodLabel,reactLabel):
    """
    Transforms a labeling mapping of type:
    
    abc --> cba
    
    to an indexed labeling of the type
    
    [1,2,3]  --> [3,2,1]
    
    Inputs:

    -) inds: numerical indices, e.g.: [1,2,3]
    -) prodLabel: product labeling, e.g. : cba
    -) reactLabel: reactant labeling, e.g.: abc
    
    """

    carbonsAll = prodLabel  
    carbons = ''.join(carbonsAll[i-1] for i in inds)        

    if not carbons:
        raise Exception('Indices: '+str(inds)+' not in product label: '+prodLabel)

    indP = []
    indR = []
    for carb in carbons:
        index = reactLabel.find(carb)+1
        if index > 0:
            indP.append(carb)
            indR.append(index)   

    success = True if indR else False

    return success,(indP,indR)


def getAtoms(string):
    """
    Extracts symbols and number from a string representing elemental composition
    
    For example:
        getAtoms('H6NO2Si') = [('H', 6), ('N', 1), ('O', 2), ('Si', 1)]
    
    """
    comps = re.findall('[A-Z][a-z]*\d*',string)
    atoms = []
    for comp in comps:
        symbol = re.match('([a-zA-Z]+)',comp).group(1)
        numberst = comp.lstrip(symbol)
        if numberst == '':
            number = 1
        else:
            number = int(numberst)
        AtTuple= symbol,number 
        atoms.append(AtTuple)
    return atoms    


def eraseFilesIn(dirname):     # TODO: check if this is needed anymore, substitute by os function?
    """
    Erases all files in a dir 
    """
    os.chdir(dirname)
    allfiles = os.listdir('.')
    for file in allfiles:
        os.remove(file)
    dir = os.getcwd()
    os.chdir('..')
    os.removedirs(dir)

            
def substituteInFile(fileIn,changes):    # Only used in one place so far
    """
    Makes substitutions in file given by the list of list changes (e.g.):

    changes = [['maxflow13CU = (\d+)','maxflow13CU = '+str(125)],['BiomassEcoli','biomassSC4']]

    """
    # Do some input check for sanity's sake
    if changes:
        assert changes.__class__.__name__ == 'list'   ,  "Changes input must be a list of lists"
        assert changes[0].__class__.__name__ == 'list',  "Changes input must be a list of lists"
    else:
        assert changes.__class__.__name__ == 'list'   ,  "Changes input must be a list of lists"        


    # Get file lines
    fileName,string = file2tuple(fileIn)
    lines = string.split('\n')
               
    newlines = lines
    # Do substitutions
    for change in changes:
        lines = newlines
        # Regular expression     
        p = re.compile(change[0])
        value = str(change[1])   # Convert into string

        # Substitution
        newlines=[]
        for line in lines:        
            newlines.append(p.sub(value, line))
    
    newString = '\n'.join(newlines)

    return (fileName,newString)               
            

def file2tuple(thing):
    """
    Function that converts a filename to a tuple of the type: (filename,string containing file)
    Input can be filename or tuple with filename and content in a string.
    """
    
    if isinstance(thing,tuple):  # it is a tuple with filename and content in a string
        filename,string = thing
        output = (filename.split('/')[-1],string)
    else:    # it is a filename
        filename = thing
        file = open(filename)
        output = (filename.split('/')[-1],file.read())
        file.close()   
        
    return output            


# Parsing functions for reactions
def getNameComp(nameComp):
    "Takes a name of the type ac[c] and separates the name(ac) from the compartment(c =cytosol)"
    nameComp.strip()
    if len(nameComp)==1:
        stName = nameComp
        comp = ''
    else:
        if nameComp[-1] == ']' and nameComp[-3] == '[': 
            stName  = nameComp[0:-3]
            comp  = nameComp[-2:-1]
        elif nameComp[-2] == '_': 
            stName  = nameComp[0:-2]
            comp  = nameComp[-1] 
        else:
            stName = nameComp
            comp = '' 
    stName = stName.strip()
    
    if ' ' in stName:    #look stochiometry info
        st,name=stName.strip().split(' ')
        try:
            stoichiometry = int(st)
        except:
            print("invalid input: "+nameComp)
    else:
        name = stName
        stoichiometry = 1
   
    return stoichiometry,name.strip(),comp


def parseLine(line):    # TODO: avoid duplication with method or integrate with cameo/cobrapy!!!!!
    """
        Parses input lines
        This function takes two possible inputs, the 13C MFA input:

        GLCpts  glc-D[e] + pep --> g6p + pyr  abcdef + ABC : abcdef + ABC             
            
        and the COBRA input:
            
        TES1[c]: aacoa --> acac + coa
            
        Eventually we should create a single type of input
    """
    
    # First detect type of input
    #       Finds wether  ':' shows up before of after the arrow symbol
    COBRA = line.rfind(':') < max(line.find('-->'),line.find('<==>'))   # TODO: this will do for the time being, but we need a unified input type 
    
    if COBRA:
        # COBRA input parsing
        # Line of type:    TES1[c]: aacoa --> acac + coa
        reversible = False
        exchange   = False
        if 'EX_' in line:   # Exchange reactions
            exchange   = False
            reversible = True
            namepcomp, rest=line.split(':')
            reacName = namepcomp.strip()
            comp = 'e'
            reacsign = '-->'
            reacLine,prodLine=rest.split(reacsign)
            reactsAll = reacLine.split('+')
            newReacts = []
            for react in reactsAll:
                stoich = 1
                name = react.strip(' ').replace('[e]','')
                newReacts.append((stoich,name,comp))
            reacts = newReacts
            prods  = []
        else:
            compAll = ''
            nameComp, rest=line.split(':')
            # Name
            [stoich,reacName,compAll] = getNameComp(nameComp.strip())   # compAll is the compartment provided for all metabolites
            # Producs and reactants
            if '<==>' in line:
                reacsign = '<==>'
                reversible = True
            else:
                reacsign = '-->'
            reacLine,prodLine=rest.split(reacsign)
            reactsAll = reacLine.split('+')
            prodsAll  = prodLine.split('+')
        
            reacts = []
            for nameComp in reactsAll:
                [stoich,name,comp] = getNameComp(nameComp.strip())
                comp = comp if comp else compAll
                reacts.append((stoich,name,comp))
            prods = []
            for nameComp in prodsAll:
                [stoich,name,comp] = getNameComp(nameComp.strip())
                comp = comp if comp else compAll
                prods.append((stoich,name,comp))
    
        return reacName, reversible, exchange, reacts, prods
    else:
        # 13C MFA input parsing
        # Line of type:    GLCpts  glc-D[e] + pep --> g6p + pyr  abcdef + ABC : abcdef + ABC  
        transLine = core.AtomTransition(line)
        exchange = False  # For the time being....
        
        # TODO: this is a temporary fix, need to redo reaction input
        reacts = []
        for reactant in transLine.reactants:
            reacts.append((reactant.stoichiometry, reactant.originalMet.name, 'c'))

        prods = []
        for product in transLine.products:
            prods.append((product.stoichiometry, product.originalMet.name, 'c'))

        return transLine.name, transLine.reversible, False, reacts, prods


def getCarbonDict(reactionList):  # TODO: Eliminate?
    """
        Provides carbon dictionary out of the formula dictionary kept in the base directory
    """
    formulaDict = DB.getFormulaDict()

    carbonDict={}
    reacDict = reactionList.getReactionDictionary()
    mets = []
    for name in reacDict:
        reaction = reacDict[name]
        if reaction.exchange:
            flux = reaction.flux
            name = reaction.reactants[0].name
            if flux != 0:
                mets.append(name)

    for name in mets:
        abbr = name
        reg = []
        if abbr in formulaDict:
            reg = re.search('C(\d*)[A-Z]',formulaDict[abbr])
        if reg:
            if reg.group(1):
                final =  int(reg.group(1))
            else:
                final = 1
        else:
            final = 0
        carbonDict[abbr] = final

    return carbonDict            
            
          
   
def fluxPlotGraph(inputList,subsystemList,titleFig='',save='default',plotType=1,axesRanges=[-30,30],legendFontSize=22):
    " Function that creates graph for toya paper"    
    # Indexing data in first method by subsystem for future reference
    inp = inputList[0]
    subsysDict = {}
    for subsys in subsystemList:
        subsysDict[subsys] = inp.getReactionSubSet(rxnNames=subsys)

    # Plotting data
    axes = figure(figsize=(8,8)).add_subplot(111)
    if plotType == 1:
        colors = ['b','k','r','g','m','c','y']
        fmts   = ['o','s','s','*','*','v','v']
        markersizes = [8,5,5,13,13,7,7]
        plt.hold(True)
        n=0
        inp = inputList[n]
        inp.plotFluxes(color=colors[n],axes=axes,fmt=fmts[n],markersize=markersizes[n],label=inp.name,axesRanges=axesRanges)
        inp.plotFluxes(color=colors[n],axes=axes,fmt=fmts[n],markersize=markersizes[n],label=inp.name,plotType='transparent std',axesRanges=axesRanges)
        for n in range(len(inputList))[1:]:
            inp = inputList[n]
            inp.plotFluxes(color=colors[n],axes=axes,fmt=fmts[n],markersize=markersizes[n],label=inp.name,plotType='transparent',axesRanges=axesRanges)
    elif plotType == 2:
        colors = ['b','k','r','g','m','c','y']
        fmts   = ['o','s','s','*','*','v','v']
        markersizes = [8,5,5,13,13,7,7]
        plt.hold(True)
        hands =[]
        labs  =[]
        for inp in inputList:
            n = inputList.index(inp)
            hands.append(inp.plotFluxes(color=colors[n],axes=axes,fmt=fmts[n],markersize=markersizes[n],label=inp.name,plotType='transparent std',axesRanges=axesRanges))
            labs.append(inp.name)
        axes.legend(hands, labs,prop={'size':legendFontSize})


    elif plotType == 3:
        colors = ['b','k','r','g','m','c','y']
        fmts   = ['o','s','s','*','*','v','v']
        markersizes = [8,5,5,13,13,7,7]
        plt.hold(True)
        for inp in inputList:
            n = inputList.index(inp)
            inp.plotFluxes(color=colors[n],axes=axes,fmt=fmts[n],markersize=markersizes[n],label=inp.name,axesRanges=axesRanges)
    # Doing axes
    font0 = FontProperties()
    font1 = font0.copy()
    font1.set_size(12)
    font1.set_weight('bold')
    xlim = 0
    if len(subsystemList)>1:
        plot([xlim,xlim],axesRanges,'k--')   # initial black separator line
        for subsystem in subsystemList:        
            # Separator plot   
            xmid = xlim + old_div(len(subsysDict[subsystem].reactions),2)     # Coordinates for subsystem acronym        
            xlim = xlim + len(subsysDict[subsystem].reactions)       # Coordinates for separator 
            plot([xlim+old_div(1,2.),xlim+old_div(1,2.)],axesRanges,'k--')                 # Black separator lines
            # Subsystem acronym plot
            name = subsystem.replace('S_','').replace('_','')
            name = ''.join(w for w in name if w.isupper())         # Keep only capital letters
            text(xmid+0.9, old_div(axesRanges[1],2), name, fontproperties=font1, horizontalalignment='center')
    else:
        xlim = xlim + len(subsysDict[subsystemList[0]].reactions)     
    
    plot([0,xlim+1],[0,0],'k:')    # Black zero level line        
    
    # Axis labels
    xlabel('Reaction number')
    ylabel('Fluxes(mmol/gdw/h)')
    
    # legend
    handles, labels = axes.get_legend_handles_labels()
                
    axes.legend(loc=3,prop={'size':legendFontSize},numpoints=1)

    # title
    title(titleFig)
    
    # Save fig        
    if save != 'default':
        savefig(save)  
    
    
def labelingGraph(inputListX,inputListY,titleFig='',labels=[],save='default',colors='default',formats='default',markersizes='default'):
    
    # Plotting
    axes = figure(figsize=(8,8)).add_subplot(111)
    axes.plot([0,1],[0,1],'k')     # Black line for comparison
    axes.axis([0,1,0,1])    
    
    if colors == 'default':
        colors = ['k','r','g','m','c','y']
    if formats == 'default':
        formats   = ['s','s','*','*','v','v']
    if markersizes=='default':    
        markersizes = [5,5,13,13,7,7]

    plt.hold(True)
    for i in range(len(inputListX)):
        Xdata = inputListX[i]
        Ydata = inputListY[i]
        allKeys = list(set(Xdata.keys()).intersection(set(Ydata.keys())))        # intersection of both keys
        XvsYcomp(Xdata,Ydata,allKeys,axes,colors[i],formats[i],label=labels[i],alpha=0.6,markersize=markersizes[i])

    title(titleFig)
    
    # legend
    axes.legend(loc=4,prop={'size':16},numpoints=1)
    
    xlabel('Experimental Labeling')
    ylabel('Predicted Labeling')    
    
    # Save fig        
    if save != 'default':
        savefig(save) 
              
def XvsYcomp(Dict1,Dict2,keys,axes,color='b',fmt='.',label='',alpha=1,markersize=10): #),titleFig=''):
    "Function which does a x vs y comparison of labeling dictionaries for given keys"

    #Limiting to given keys    
    DictX = dict((k, Dict1[k]) for k in keys if k in Dict1)    
    DictY = dict((k, Dict2[k]) for k in keys if k in Dict2)
    
    #Obtaining plot data
    x=[]  # Computational data (EMUlabel)
    y=[]  # Experimental data (fragDict)  

    for emu in DictX:
        frag = DictX[emu]
        mdv  = frag.mdv
        #for m in range(len(mdv)):
        for m in range(frag.nmeas):
            x.append(mdv[m])
            y.append(DictY[emu][m])

    #Plotting
    x = numpy.array(x)
    y = numpy.array(y)
    axes.plot(x,y,color+fmt,label=label,alpha=alpha,markersize=markersize)    
    
#    title(titleFig)
    



def histcomp(Dict1,Dict2,fragments,outputFileName,titleFig,OF,save):
    "Function which does a histogram comparison of labeling dictionaries for given keys"

 # Plotting each fragment (and output to file)
    outputFile = open(outputFileName,'w')
    width = 0.35
    fig = figure(figsize=(20,10))
    fig.suptitle(titleFig+" (OF= "+str(OF)+"), Exp(red)/Comp(blue)",fontsize=20,fontweight='bold')

    ymax  = old_div(int(len(fragments)),5)+1
    index = 1
    for fragment in sorted(fragments):
        if fragment in Dict1 and fragment in Dict2:
            if Dict2[fragment].inFit:   # Background difference depending on whether the fragment is in the fit
                BGcolor = [1,1,1]
            else:
                BGcolor = [0.92,1,0.92]

            ax = fig.add_subplot(ymax,5,index,axisbg=BGcolor)
            title(fragment,fontsize=16)
            mdv = Dict2[fragment].mdv
            h1 = mdv
            h2 = Dict1[fragment][0:len(h1)]

            ax.bar(numpy.arange(len(h1))-0.2,h1,width,color='r')
            ax.bar(numpy.arange(len(h2))+0.2,h2,width,color='b')
            # Change label ticks
            ind = list(range(len(mdv)))
            ax.set_xticks(ind)
            ax.set_xticklabels(ind)
            # Labels            
            xlabel('m',fontsize=18)
            ylabel('MDV(m)',fontsize=18)
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontsize(16)
            index = index+1

            # output to file
            outputFile.write("fragment:\n")
            outputFile.write(fragment+"\n")
            for m in range(min(len(h1),len(h2))):
                outputFile.write(str(m)+": "+str(h1[m])+","+str(h2[m])+"\n")

    fig.tight_layout()
    fig.subplots_adjust(top=0.9)

    # Save fig        
    if save != 'default':
        savefig(save)





    
        
          
############### Tests ##################


if __name__ == "__main__":
    
    dir1 = '/scratch/hgmartin_gams_files/tests/batch/20120327T104852GAMS/Run0_1/'
    dir2 = '/scratch/hgmartin_gams_files/tests/amittest/'
    
