# Part of the JBEI Quantitative Metabolic Modeling Library (JQMM)
# Copyright (c) 2016, The Regents of the University of California.
# For licensing details see "license.txt" and "legal.txt".
"""
The GAMSclasses module is the connection to the GAMS system (General Algebraic Modeling System, https://www.gams.com/) used to solve optimizatin problems.
This module should soon be adapted to the GAMS python API: https://www.gams.com/help/index.jsp?topic=%2Fgams.doc%2Fapis%2Findex.html
"""
from __future__ import print_function

from builtins import str
from builtins import range
from builtins import object
import os,time,re,numpy,subprocess, unittest
from utilities import file2tuple
import core, enhancedLists

# GAMS sets, parameters and tables form the basis of problems defined in GAMS. 
# Class for GAMS sets
class GAMSSet(object):
    def __init__(self,name,inset):
        # Checking types
        if not isinstance(name,str):
            raise Exception('First input must be a string')
        if not isinstance(inset,set):
            raise Exception('Second input must be a set')
        # Assigning
        self.set  = inset
        self.name = name


    # Read input file    
    def read(self,filename):
        fileName, string = file2tuple(filename)    
        lines = string.split('\n')        
        newSet = []
        for line in lines:
            name = line.rstrip()
            name = name.rstrip('/')
            name = name.strip('\'')
            if not name=="":
                newSet.append(name)            
        self.set = set(newSet)

    # Write set to a file
    def write(self,filename):

        # Sort function definition        
        def compare_strings(s):
            try:
                float(s)
                return float(s)
            except ValueError:
                return s

        # Get all lines
        lines = ['/ \n']
        for element in sorted(self.set,key=compare_strings):
            stringout = "'"+element+"'"+"\n"
            lines.extend(stringout)
        lines.append('/ \n')
        
         # Output
        if filename == 'toString':        # To string
            return ''.join(lines)
        else:                             # To file
            outfile = open(filename,'w')
            for line in lines:
                outfile.write(line)
            outfile.close()


# Class for GAMS parameters
class GAMSParameter(object):
    def __init__(self, name, elements):
        # Checking types
        if not isinstance(name,str):
            raise Exception('First input must be a string')
        self.name     = name
        # Change format to standard, in terms of tuples (only used in practice for single element parameters)
        newElements = {}
        for element in elements:
            if isinstance(element, tuple):
                newElements[element] = elements[element]
            else:
                newElements[tuple([element])] = elements[element]
        self.elements = newElements


    # Read parameter from a file    
    def read(self,filename):
        # Regular expression
        p = re.compile('(\w+)\(([\'\,\w\s\-\(\)\[\]]+)\)\s*=\s*([\-\.\d\w\+]+);*')
        
        fileName, string = file2tuple(filename)    
        lines = string.split('\n')          
        
        # Opening file and parsing content
        Elements = {}
        for line in lines:
            if line.rstrip():
                a = p.match(line)
                # Name of parameter
                parname = a.group(1)
                if not self.name == parname:
                    self.name = parname
                # Elements of parameter 
                #   Keys    
                namesA = a.group(2).split(',')     
                namesB =[]
                for name in namesA:
                    namesB.append(name.strip('\'').strip()) 
                #   Value
                value = float(a.group(3))
                Elements[tuple(namesB)]= value

        self.elements = Elements


    # Write parameter to a file    
    def write(self, filename):
        # Getting lines
        parname  = self.name
        Elements = self.elements
        lines    = []
        
#        for element in sorted(Elements,key=compare_elements):
        for element in sorted(Elements, key= lambda elem: '*'.join(list(elem))):
            middle = "\'"+'\',\''.join(list(element))+"\'"
            string = parname + "("+middle+")="+str(Elements[element])+";\n"
            lines.append(string)

        # Output
        if filename == 'toString':        # To string
            return ''.join(lines)
        else:                             # To file
            outfile = open(filename,'w')
            for line in lines:
                outfile.write(line)
            outfile.close()


def GAMSParameterFromResult(input):
    """
    Function that provides a GAMSParameter based on results output from a GAMSProblem.
    Example usage:
    outputFuncs = {'Vmax':(GAMSclasses.fromResult,['Vmaxout.txt','Vmax']), 
                   'Vmin':(GAMSclasses.fromResult,['Vminout.txt','Vmin']), 
                  }
    Note: To support pickling without a 'TypeError: can't pickle instancemethod objects',
    this method cannot be a classmethod of GAMSParameter    
    """
    
    fileName = input[0]
    parName  = input[1]
    p = GAMSParameter(parName,[])
    p.read(fileName)
    return p



class GAMSTable(object):
    "Class for GAMS tables"
    def __init__(self, name, xkeys, ykeys, array):
        # Checking types 
        if not isinstance(name,str):
            raise Exception('First input (name) must be a string')
        if not isinstance(xkeys, GAMSTableKeys):
            raise Exception('Second input must be a GAMSTableKeys object')
        if not isinstance(ykeys, GAMSTableKeys):
            raise Exception('Third input must be a GAMSTableKeys object')
        if not isinstance(array, numpy.ndarray):
            raise Exception('Fourth input must be an array')
        else:
            for row in array:
                for value in row:
                    if not isinstance(value,int) and not isinstance(value,float):
                        raise Exception('Array must be composed of floats or integers')
        # Assigning         
        self.name  = name
        self.xkeys = xkeys
        self.ykeys = ykeys
        self.array = array


    def write(self, filename):    
        # Unpacking
        xkeys = self.xkeys
        ykeys = self.ykeys
        oldarray = self.array

        # Sorting y keys
        ykeyDict = {}
        for i in range(len(ykeys.values)):
            ykeyDict[ykeys.values[i]] = oldarray[i]
            
        array = []  
        ykeys.values = sorted(ykeys.values)
        for ykey in ykeys.values:
            array.append(ykeyDict[ykey])
            
        array = numpy.array(array)    

        # storing everything in lines
        lines = []
        # headers
        head1    = "Table "+self.name+"("+ykeys.name+","+xkeys.name+")\n"
        head2 = '\t'+'\t'.join(xkeys.values)+"\n"
        lines.append(head1)
        lines.append(head2)
        # Array
        for i in range(len(ykeys.values)):
            valuelist=[]
            for value in array[i]:
                if value>=10:           # trying to avoid here destruction of table format. There is quite probably a more elegant way to do it
                    valuelist.append('{0:4.4f}'.format(value))                    
                else:
                    valuelist.append('{0:5.5f}'.format(value))
            stringout = ykeys.values[i]+"\t"+"\t".join(valuelist)+"\n"
            lines.append(stringout)
        lines.append(";")

        # Writing file
        if filename == 'toString':
            output = ''.join(lines)
        else:
            fileout  = open(filename,'w')
            for line in lines:
                fileout.write(line)
            fileout.close()
            output = ''
        
        return output



class GAMSTableKeys(object):
    "Auxiliary class for table keys for GAMS table"
    def __init__(self, name, values):
        # Checking types 
        if not isinstance(name, str):
            raise Exception('First input (name) must be a string')
        if not isinstance(values, list):
            raise Exception('Second input (values) must be a list')
        # Assigning
        self.name   = name
        self.values = values



# Class for GAMS problems
globalGAMSRunCount = 0   # Global variable so as to get unique directory names

class GAMSProblem(object):
    def __init__(self,name,gamsfile,includeFiles,outputFuncs,execType='parallel',directory='',erase=False):
        "gamsfile or includeFiles can be either filenames or tuples of filenames and strings with file content: (filename, string)"
        # Storing
        self.name         = name          # Name of GAMS problem
        self.outputFuncs  = outputFuncs   # Output functions for processing results
        self.type         = execType      # Type of execution process: in parallel or series 
        self.dir          = directory     # Directory for problem execution
        self.erase        = erase         # Flag signaling for erasure of GAMS output files after solver completes

        # Storing files as tuples
        # gamsfile
        self.gamsFile = file2tuple(gamsfile)
        # Include files
        newIncludeFiles = [file2tuple(thing) for thing in includeFiles]
        
                        
        self.includeFiles = newIncludeFiles


    def __enter__(self):
        "Used with 'with' statement to set problem up automatically"
        self.writeAll()
        
        return self


    def __exit__(self, type, value, traceback):
        "Used with 'with' statement to clean problem up automatically"     
        # Not supressing any errors (i.e. doing nothing with type, value or traceback)
              
        # Return to original directory
        os.chdir(self.origDir)
        
        # Erase files and base dir
        if self.erase:    
            self.eraseFiles()


    # Write gams file to outfilename        
    def writeGAMSFile(self,outfilename):  
        output = []
        if outfilename == 'toString':
            output = self.gamsFile
        else:
            outfile  = open(outfilename,'w') 
            outfile.write(self.gamsFile[1])
            outfile.close()
            
        return output


    # Write include files needed for GAMS problem
    def writeIncludeFiles(self,toString=False):
        output = []
        if toString:
            output = self.includeFiles
        else:
            for filename,string in self.includeFiles:
                outfile = open(filename,'w')   
                outfile.write(string)
                outfile.close()
        return output


    # Temporary folder name for file holding and execution.
    @staticmethod
    def getGAMSTempFolderName():
        global globalGAMSRunCount
        globalGAMSRunCount = globalGAMSRunCount + 1
        return 'gams-' + time.strftime("d%Y%m%d-t%H%M%S")+'-p'+str(globalGAMSRunCount)


    # Write all files
    def writeAll(self,dirName = 'default'):
        output = []
        toString = False
        self.origDir = os.getcwd()
        if dirName == 'toString':     # Output as strings
            toString = True
            output = [self.writeGAMSFile('toString')]
            output.extend(self.writeIncludeFiles(toString))
        else:                         # Write to dir
            # Decide dir name to store files
            if dirName == 'default':
                # Name of dir is just date,time and microseconds plus GAMS string
                tempFolder = GAMSProblem.getGAMSTempFolderName()
                dirName  = os.getcwd()+'/'+tempFolder
                self.dir = dirName
            else:
                self.dir = dirName
            # Create directory and move in there
#            try:
            os.mkdir(dirName)
#            except OSError as e:
#                if e.errno == 17:                    
#                    pass  # if  directory already exists. NOTE: This is potentially a dangerous situation and should be an error!
                    
            # Change to dir to write files        
            os.chdir(dirName)
            # Write  all files
            self.writeGAMSFile(self.gamsFile[0])
            self.writeIncludeFiles(toString)
        
            # Going back to original directory
            os.chdir(self.origDir)     
            
        return output


    # Run problem
    def run(self):
        # Determining gams command depending on problem type
        if self.type == 'serial':
            command     = 'gams'
            output      = self.gamsFile[0].replace('gms','log')
        elif self.type == 'parallel':
            command = 'gamsbatch'    
            output  = "/dev/null 2>&1"  # throw output out (details in .log file)
        else:
            raise Exception('GAMS problem type: ' + self.type + ' not recognized') 

        # Dir to store files
        try:
            dirName = self.dir
        except AttributeError:
            raise Exception('Files not written yet')
        os.chdir(dirName)
                
        # Run gams process
        currdir = os.getcwd()+'/'
        options = ' Workdir=' + currdir + ' ps=9999 pw=120 MaxProcDir=1000 lo=2 lf='+currdir+'logfile';
        call    = command+' '+ self.gamsFile[0] + options +" > " + output  # throw output out
 
        status  = os.system(call)
        if status != 0:
            print("call:")
            print(call)
            raise Exception('GAMS process launch failed')
        
        # Flagging as started
        self.started = True
        # Going back to original directory
        os.chdir(self.origDir)


    # Check if GAMS solver is done
    def check(self):
        try:
            started = self.started
        except AttributeError:
            raise Exception('GAMS process not started yet')

        if self.type == 'serial':
            # Serial problems are not run in batch
            finished = True
        elif self.type == 'parallel':
            # Note: Inspecting a logfile on a timer is a poor way to track progress of subprocesses.
            # An example of a better way would be to use os.fork and check the status of PIDs.

            # Change to right dir
            dirName = self.dir
            os.chdir(dirName)

            # Check if STATUS is in .log file
            gamsFile = self.gamsFile[0]
            gamsLog = gamsFile.split('.')[0]+'.log'
            command = ["grep","Status",dirName+'/'+gamsLog]
            grepLine = subprocess.Popen(command,stdout=subprocess.PIPE).communicate()[0]
            if (grepLine == ''):
                finished = False
            else:
                finished = True
        else:
            raise Exception('GAMS problem type: ' + self.type + ' not recognized') 
            
        # Going back to original directory
        os.chdir(self.origDir)
        
        return finished


    # Wait until GAMS solver is done
    def waitTilDone(self,pauseTime=0.001):
        # Check problem is finished
        ready     = False
        while ready == False:
            ready = self.check()
            time.sleep(pauseTime)  


    # Collect results
    def collect(self):
        # Change to right dir
        try:
            dirName = self.dir
        except AttributeError:
            raise Exception('Process not started yet')
        os.chdir(dirName)

        # Output functions
        outputFuncs = self.outputFuncs
        output = {}
        for name in outputFuncs:
            function  = outputFuncs[name][0]
            arguments = outputFuncs[name][1]
            output[name] = function(arguments)

        # Going back to original directory
        os.chdir(self.origDir)
        
        self.Results = output


    # Erase files 
    def eraseFiles(self):
        # Change to right dir
        try:
            dirName = self.dir
        except AttributeError:
            raise Exception('Process not started yet')  

        os.chdir(dirName)
        
        # Erase all files and dir
        allfiles = os.listdir('.')
        
        for fileName in allfiles:
            try:
                os.remove(fileName)
            except OSError as e:
                if e.errno == 2:  # 'file not found' error
                    pass
                elif e.errno == 21: # 'file is a directory' error
                    time.sleep(0.1)
                else:
                    raise 
                    
        dirOut = os.getcwd()
        
        os.chdir('..')
        os.rmdir(dirOut)
        
        # Going back to original directory          
        os.chdir(self.origDir)
    
    

# Class for a batch of GAMS problems
class GAMSProblemBatch(object):
    def __init__(self, gamsProblems, erase=True):
        self.problems = gamsProblems
        self.erase = erase

    def __enter__(self):
        "Used with 'with' statement to set problem up automatically"
        self.writeFiles()
        
        return self
        
    
    def __exit__(self, type, value, traceback):
        "Used with 'with' statement to clean problem up automatically"
        # Note: Using __exit__ (to automatically run code when an object drops out of scope)
        # may not be the best way to delete files from the filesystem, especially when it
        # requires changing the working directory, which can affect other file operations.

        # Return to original directory              
        os.chdir(self.origDir)
        
        # Erase files and base dir
        if self.erase:    
            self.eraseFiles()    


    # Run job batch    
    def run(self):
        problems = self.problems
        for problem in problems:
            problem.run() 


    # Write all needed files
    def writeFiles(self):
        self.origDir = os.getcwd()       # This is the original directory for the BATCH, do not confuse with that of each problem
        for problem in self.problems:
            if problem.dir:
                problem.writeAll(problem.dir)
            else:
                problem.writeAll()


    # Check jobs are done
    def check(self):
        for problem in self.problems:
            if not problem.check():
                return False
        return True


    # Wait until all jobs are done
    def waitTilDone(self):
        pauseTime = 1
        ready     = False
        while ready == False:
            ready = self.check()
            time.sleep(pauseTime)


    # Collect data        
    def collect(self):
        problems = self.problems
        for problem in problems:
            problem.collect()


    # Return results in a list
    def getAllResults(self):
        problems    = self.problems
        resultsDict = {}
        for problem in problems:
            resultsDict[problem.name] = problem.Results
        
        return resultsDict    


    # Erase all data
    def eraseFiles(self):
        for problem in self.problems:
            problem.eraseFiles()    



##### UNIT TESTS ########

# Expected result from tests below
expectedResult =  '153.68'  


class testGAMSInstallation(unittest.TestCase):  
    """ Testing that GAMS is installed """

    def testGAMS(self):
        """ Checks if gams can be called at the command line successfully"""
        
        status = os.system('gams')
        errmsg = 'GAMS can not be called at the command line succesfully. \n'

        self.assertTrue(status==0,msg=errmsg)
        

class testGAMSProblem(unittest.TestCase):  
    """ Testing that GAMS problem class provides recorded result"""

    def setUp(self):
         pass  
           
           
    def testProblem(self):
        """ Tests the classic transport problem in GAMS"""
        
        # Get GAMS problem        
        testProblem = getTransportProblem()
        
        # Solve testProblem 
        with testProblem as prob:
            # Run problem
            prob.run()

            # Check if it is done
            prob.waitTilDone()

            # Collect data
            prob.collect()
        
            # Output data 
            resultsDict = prob.Results  
       
        # Extract expected results for Z and compare
        match =  re.search('Z    =       (\d*\.\d*)\n',resultsDict['Zout'])
        self.assertTrue(match.group(1)==expectedResult)       
                    
    def tearDown(self):
        pass    


class testGAMSBatch(unittest.TestCase):  
    """ Testing that GAMS problem batch class provides recorded results"""

    def setUp(self):
        pass  
           
    def testBatch(self):
        """ Tests the classic transport problem in GAMS in batch mode"""
        
        # Get GAMS problems        
        testProblem1 = getTransportProblem('transport1','dir1')
        testProblem2 = getTransportProblem('transport2','dir2')

        testBatch = GAMSProblemBatch([testProblem1,testProblem2],erase=True)
       
       
        # Solve testBatch    
        with testBatch as batch:
            # Run problem(s)
            batch.run()

            # Check if they are done
            batch.waitTilDone()
        
            # Collect data
            batch.collect()
        
            # Put results in local format 
            resultsDict = batch.getAllResults()
       
        # Extract expected results for Z and compare
        match1 =  re.search('Z    =       (\d*\.\d*)\n',resultsDict['transport1']['Zout'])
        match2 =  re.search('Z    =       (\d*\.\d*)\n',resultsDict['transport2']['Zout'])
        
        self.assertTrue(match1.group(1)==expectedResult)  
        self.assertTrue(match2.group(1)==expectedResult)  
                    
    def tearDown(self):
        pass    


### Auxiliary functions for tests

def getTransportLines():
    """
    This function provides the GAMS definition of the classic Dantzig transport problem
        https://www.gams.com/docs/example.htm    
    """
    
    transportLines=[
    'SETS',
    'I   canning plants   / SEATTLE, SAN-DIEGO /',
    'J   markets          / NEW-YORK, CHICAGO, TOPEKA / ;',
    'PARAMETERS',
    'A(I)  capacity of plant i in cases',
    '         /    SEATTLE     350',
    '              SAN-DIEGO   600  /',
    'B(J)  demand at market j in cases',
    '         /    NEW-YORK    325',
    '              CHICAGO     300',
    '              TOPEKA      275  / ;',
    'TABLE D(I,J)  distance in thousands of miles',
    '                    NEW-YORK       CHICAGO      TOPEKA',
    '      SEATTLE          2.5           1.7          1.8   ',
    '      SAN-DIEGO        2.5           1.8          1.4  ;',
    'SCALAR F  freight in dollars per case per thousand miles  /90/ ;',
    'PARAMETER C(I,J)  transport cost in thousands of dollars per case ;',
    '       C(I,J) = F * D(I,J) / 1000 ;',
    'VARIABLES',
    '       X(I,J)  shipment quantities in cases',
    '       Z       total transportation costs in thousands of dollars ;',
    'POSITIVE VARIABLE X ;',
    'EQUATIONS',
    '       COST        define objective function',
    '       SUPPLY(I)   observe supply limit at plant i',
    '       DEMAND(J)   satisfy demand at market j ;',
    'COST ..        Z  =E=  SUM((I,J), C(I,J)*X(I,J)) ;',
    'SUPPLY(I) ..   SUM(J, X(I,J))  =L=  A(I) ;',
    'DEMAND(J) ..   SUM(I, X(I,J))  =G=  B(J) ;',
    'MODEL TRANSPORT /ALL/ ;',
    'SOLVE TRANSPORT USING LP MINIMIZING Z ;',
    '',
    'file foutfile/"Zout.txt"/;',
    'put foutfile',
    'put "Z    = ",Z.l," ";    put /;',    
    'putclose foutfile;'
    ]
    
    return transportLines
    
    
def getTransportProblem(name='transport',directory='dir1'):
    """Obtain GAMS transport problem"""
    
    # Get GAMS description of problem
    transportLines = getTransportLines()
    
    # Functions for output data
    outputFuncs = {'Zout'           :(getZout           ,['Zout.txt'])}
        
    # Start GAMS problem
    testProblem = GAMSProblem(name,('transport.gms','\n'.join(transportLines)),[], outputFuncs,execType='serial',directory=directory,erase=True)
    
    return testProblem


def getZout(input_):
    "Auxiliary function for test"
    
    fileName = input_[0]
    fileInfo = open(fileName)
    
    return ''.join(fileInfo)




# Flux GAMS parameter class (only adds an output method)
class fluxGAMSpar(GAMSParameter):
    
    def __init__(self,name,elements,reactionList='default',fluxDict='default'):
        
        if reactionList != 'default':
            # if reactionList present, use that to initalize parameter elements
            elements = {}
            for reaction in reactionList:
                if reaction.reversible and not reaction.exchange:
                    elements[reaction.name+'_f'] = reaction.flux.forward.best
                    elements[reaction.name+'_b'] = reaction.flux.backward.best                        
                else:
                    elements[reaction.name] = reaction.flux.net.best
        elif fluxDict != 'default':
            elements = {}
            for fluxName in fluxDict:
                try:
                    elements[fluxName] = fluxDict[fluxName].net.best      # This should be generalized
                except AttributeError:
                    try:
                        elements[fluxName] = fluxDict[fluxName].net      
                    except AttributeError:
                            elements[fluxName] = fluxDict[fluxName]                                 
        
        GAMSParameter.__init__(self,name,elements)


    def getReactionList(self):
        reactions = []
        for rxnTup in self.elements:
            name    = rxnTup[0]
            forward = self.elements[rxnTup]
            # Find if it is reversible
            if name[-2:] != '_b':
                if name[-2:] == '_f':
                    # Find backward reaction
                    name = name[0:-2]
                    backward = self.elements[tuple([name+'_b'])]
                    rev = True
                else:
                    backward = 0
                    rev = False
                flux = core.flux(for_back_tup=(forward,backward))
                reaction = core.Reaction(name,reversible=rev,flux=flux)
                reactions.append(reaction)
                
        return enhancedLists.ReactionList(reactions)


    def compare(self,fluxDict):
        pass
        


class labelingGAMSpar(GAMSParameter):
    "GAMS parameter for labeling"

    def getLabelingDict(self,fragDict):
        EMUlabel = self.elements
        EMULabelDict = {}
        maxLength    = {}
        # Gathering data
        for tup in EMUlabel:
            emu = tup[0]
            m   = int(tup[1])
            frag= fragDict[emu]
            nmeas = frag.nmeas
            if m <= nmeas-1:
                if emu not in EMULabelDict:
                    EMULabelDict[emu] = []
                    maxLength[emu] = 0
                EMULabelDict[emu].append((m,EMUlabel[tup]))
                maxLength[emu] = max(maxLength[emu],m)
       
        # Reordering data
        for emu in EMULabelDict:
            old = EMULabelDict[emu]
            EMULabelDict[emu] = numpy.zeros(maxLength[emu]+1)
            for tuple in old:
                m   = tuple[0]
                val = tuple[1]
                EMULabelDict[emu][m] = val
        
        return EMULabelDict



class GAMSProblemInfo(object):
    " Class for info in GAMSProblem "

    def __init__(self,OF,solvestat,modelstat,numequ,numnz,numvar,resusd):
        self.OF        = OF
        self.solvestat = solvestat
        self.modelstat = modelstat
        self.numequ    = numequ
        self.numnz     = numnz
        self.numvar    = numvar
        self.resusd    = resusd


