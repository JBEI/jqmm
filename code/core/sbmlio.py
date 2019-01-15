# Part of the JBEI Quantitative Metabolic Modeling Library (JQMM)
# Copyright (c) 2016, The Regents of the University of California.
# For licensing details see "license.txt" and "legal.txt".
"""
The sbmlio module provides support for input/output of ReactionNetwork structures in SBML format.
"""

from builtins import str
from builtins import range
from past.builtins import basestring
from builtins import object
import re
import core, DB, enhancedLists, labeling
import Genes, Proteins
import utilities as utils
import libsbml


class SBMLImporter(object):
    """
    Helper class to parse the given string as SBML and extract various structures important to a ReactionNetwork.
    """
    def __init__(self, sbmlFileContents=None, exceptionOnError=True):

        self.model_name = None
        self.model_notes = {}

        self.metList = None
        self.reactionList = None
        self.masterGeneSet = None
        self.masterProteinSet = None

        self.initialSetFluxesCount = None
        self.initialBoundedFluxesCount = None
        self.initialCloselyBoundedFluxesCount = None

        self.errorCount = 0
        self.logStrings = []

        if not isinstance(sbmlFileContents, basestring):
            self.errorCount = 1
            self.logStrings = ['No SBML string provided (object is not a string).']
        elif sbmlFileContents.strip() == '':
            self.errorCount = 1
            self.logStrings = ['No SBML string provided (string is empty).']
        else:
            # Not saving a local pointer to sbmlDoc, in order to circumvent this error when attempting to deepcopy:
            # TypeError: object.__new__(SwigPyObject) is not safe, use SwigPyObject.__new__()
            sbmlDoc = libsbml.readSBMLFromString(sbmlFileContents)
            self.errorCount = sbmlDoc.getNumErrors()
            if self.errorCount:
                errs = sbmlDoc.getErrorLog()
                for x in range(self.errorCount):
                    err = errs.getError(x)
                    self.logStrings.append(err.getMessage())
            self.parseSBMLDoc(sbmlDoc)

        if self.errorCount > 0 and exceptionOnError:
            raise Exception('Errors parsing SBML:\n%s' % ('\n'.join(self.logStrings)))


    def parseSBMLDoc(self, sbmlDoc):
        model = sbmlDoc.getModel()
        self.model_name = model.getName()
        self.logStrings.append('Parsed model name ' + self.model_name)

        # ###### Get notes (= everything stored in SBML as notes: GCMS labeling and so on...)
        self.model_notes = convertSBMLNotes(notesSBML=model.getNotes())

        # Change labeling abbreviations to those found in fragDictBasic. This should probably go into its own method of a Notes class
        types = ['GCMSLabelData','GCMSLabelData(fit)','LCMSLabelData','LCMSLabelData(fit)','CEMSLabelData','CEMSLabelData(fit)']  #TODO: create a more general way to do this
        for t in types:
            if t in self.model_notes:
                MDVDict = self.model_notes[t]
                for abbrev in MDVDict:
                    labelData = MDVDict[abbrev]
                    abbrevNew = labelData.getAbbrev(abbrev,standard = True)
                    if abbrevNew != abbrev:
                        MDVDict[abbrevNew] = MDVDict[abbrev]
                        MDVDict.pop(abbrev)

        # ####### Species
        listOfSpecies = model.getListOfSpecies()

        mets = []
        anyFeed   = False
        anySource = False

        for sp in listOfSpecies:
            # Changing met name to appropiate form
            species_id  = re.sub('^M_','', sp.getId()) # ID must exist in valid SBML   
            full_name   = sp.getName()
            if not isinstance(full_name, str):
                full_name = species_id
            elif full_name.strip() == '':
                full_name = species_id
            else:
                full_name = re.sub('^M_','', full_name)

            ncarbons    = 0
            excluded    = False
            source      = False
            feed        = False
            destination = False
            carbonsOut  = None
            formula     = ''

            # Note: Converts all names to SBML
            spName     = core.MetaboliteName(species_id)

            spNotes = {}
            # Getting metabolite names for each category
            if sp.isSetNotes():
                notes    = sp.getNotes()
                spNotes  = convertSBMLNotes(notesSBML=notes)

                destination   = True                       if 'LABELING_DEST'   in spNotes else destination
                carbonsOut    = spNotes['LABELING_DEST']   if 'LABELING_DEST'   in spNotes else carbonsOut                        
                source        = spNotes['LABELING_SOURCE'] if 'LABELING_SOURCE' in spNotes else source
                excluded      = spNotes['EXCLUDED']        if 'EXCLUDED'        in spNotes else excluded
                feed          = spNotes['FEED']            if 'FEED'            in spNotes else feed
                ncarbons      = spNotes['CARBON NUMBER']   if 'CARBON NUMBER'   in spNotes else ncarbons
                formula       = spNotes['FORMULA']         if 'FORMULA'         in spNotes else formula

                anySource = True if source else anySource
                anyFeed   = True if feed else anyFeed
          
            met = core.Metabolite(spName, ncarbons, excluded, source, feed, destination, carbonsOut, full_name=full_name, formula=formula, extraNotes=spNotes)
            mets.append(met)

        self.logStrings.append('Detected ' + str(len(listOfSpecies)) + ' species, parsed ' + str(len(mets)) + ' metabolites')
            
        # Creating met list    
        metList = enhancedLists.MetaboliteList(mets)

        # Warning if no feed
        if not anyFeed and anySource:
            self.logStrings.append('Feed not present!')

        # ####### Reactions
        listOfReactions = model.getListOfReactions()
        Cvec            = {}
        metabolite_dict = metList.metDict

        set_fluxes_count = 0
        bounded_fluxes_count = 0
        closely_bounded_fluxes_count = 0

        reactions = []
        masterGeneSet    = Genes.GeneSet()
        masterProteinSet = Proteins.ProteinSet()

        # For each reaction in the list of reactions, extract the associated data
        for rxn in listOfReactions:
            rxnNameOrig     = rxn.getId().replace('R_','',1).replace('_e_','(e)',1)
            rxnName         = core.ReactionName(rxnNameOrig)
            exchangeReaction = True if rxnName[0:3] == 'EX_' else False 

            # Stoichiometry (reactants and products)
            sbml_reactants = rxn.getListOfReactants()
            reacts = []
            for sbml_reactant in sbml_reactants:
                # Convert the name to the standard format we expect, by instantiating it as a MetaboliteName
                reactant_name = core.MetaboliteName(sbml_reactant.getSpecies().replace('M_','',1))
                metabolite = metabolite_dict[reactant_name]
                react = core.Reactant(metabolite, sbml_reactant.getStoichiometry())
                reacts.append(react)

            sbml_products = rxn.getListOfProducts()
            prods =[]
            for sbml_product in sbml_products:
                product_name = core.MetaboliteName(sbml_product.getSpecies().replace('M_','',1))
                metabolite = metabolite_dict[product_name]
                prod = core.Product(metabolite, sbml_product.getStoichiometry())
                prods.append(prod)
                    
            #  Carbon Transitions and other notes including genes and protein sets
            transitionLine = 'None'
            subsystem = 'None'
            flux = None
            geneString = None
            geneValueString = None
            proteinString = None
            proteinValueString = None
            rxnNotes = {}
            if rxn.isSetNotes():
                rxnNotesSBML = rxn.getNotes()
                rxnNotes     = convertSBMLNotes(notesSBML=rxnNotesSBML)
                if 'CARBON_TRANSITIONS' in rxnNotes:
                    tl = rxnNotes['CARBON_TRANSITIONS'].replace('==','<==>').replace('--','-->')
                    transitionsLine = core.AtomTransition(rxnName + "\t" + tl)
                    transitionLine = transitionsLine.fullLine
                if 'FLUX VALUE' in rxnNotes:
                    flux = rxnNotes['FLUX VALUE']
                    set_fluxes_count += 1            
                if 'SUBSYSTEM' in rxnNotes:
                    subsystem = rxnNotes['SUBSYSTEM']
                if 'GENE_ASSOCIATION' in rxnNotes:
                    geneString = rxnNotes['GENE_ASSOCIATION']
                if 'GENE_TRANSCRIPTION_VALUES' in rxnNotes:
                    geneValueString = rxnNotes['GENE_TRANSCRIPTION_VALUES']
                if 'PROTEIN_ASSOCIATION' in rxnNotes:
                    proteinString = rxnNotes['PROTEIN_ASSOCIATION']
                if 'PROTEIN_COPY_VALUES' in rxnNotes:
                    proteinValueString = rxnNotes['PROTEIN_COPY_VALUES']

            geneSets = []
            if geneString:
                geneString = re.sub('^\s*:','', geneString) # Chop off any bogus colon
                geneString = geneString.strip()
                if geneString.lower() in ["none", "n.a.", "n/a"]:
                    geneString = ""
                if geneString != "":
                    # is a list of GeneSets
                    geneSets = Genes.GeneSet.createSetsFromStrings(geneString, geneValueString)
                    for geneSet in geneSets:
                        masterGeneSet.recastSet(geneSet)

            proteinSets = []
            if proteinString:
                proteinString = re.sub('^\s*:','', proteinString) # Chop off any bogus colon
                proteinString = proteinString.strip()
                if proteinString.lower() in ["none", "n.a.", "n/a"]:
                    proteinString = ""
                if proteinString != "":
                    # is a list of ProteinSets
                    proteinSets = Proteins.ProteinSet.createSetsFromStrings(proteinString, geneString, proteinValueString)
                    for proteinSet in proteinSets:
                        masterProteinSet.recastSet(proteinSet)

            #  Measured fluxes (MeasFluxes for 13CMFA, UB and LB for 2S-13CMFA or FBA)
            cval = 0
            LB = None
            UB = None

            kLaw = rxn.getKineticLaw()
            if kLaw:
                listOfParameters = kLaw.getListOfParameters()
                for par in listOfParameters:
                    if par.getId() == 'LOWER_BOUND':
                        LB = par.getValue()
                    elif par.getId() == 'UPPER_BOUND':
                        UB = par.getValue()
                    elif par.getId() == 'OBJECTIVE_COEFFICIENT':
                        cval = par.getValue()
                        Cvec[rxnName] = par.getValue()
                if (LB != None) and (UB != None):
                    bounded_fluxes_count += 1
                    if (LB == UB):
                        closely_bounded_fluxes_count += 1
                
            measured,fluxBounds = core.fluxBounds(LB, UB, rxn.getReversible(), exchangeReaction)

            # If we don't have a determined flux, use the flux bounds.
            if not flux:
                flux = fluxBounds


            # append gene and protein sets to the master gene and protein sets
            # !!!!!!!ALERT!!!! masterGeneSet is a list of list of geneSets now. Should be a geneSet with all genes in it
            # masterGeneSet.append(geneSets)
            # masterProteinSet.append(proteinSets)

            # Instantiating reaction object
            new_reaction = core.Reaction(rxnName, rxn.getReversible(), False, measured, flux=flux, fluxBounds=fluxBounds,
                                     transitionLine=transitionLine, exchange=exchangeReaction, reactants=reacts, products=prods,
                                     cval=cval, subsystem=subsystem, gene=geneSets, protein=proteinSets, extraNotes=rxnNotes)
            reactions.append(new_reaction)

        self.metList      = metList
        self.reactionList = enhancedLists.ReactionList(reactions)
        self.masterGeneSet = masterGeneSet
        self.masterProteinSet = masterProteinSet

        self.initialSetFluxesCount     = set_fluxes_count
        self.initialBoundedFluxesCount = bounded_fluxes_count
        self.initialCloselyBoundedFluxesCount = closely_bounded_fluxes_count

        self.logStrings.append('Detected ' + str(len(listOfReactions)) + ' reactions, parsed ' + str(len(reactions)) + '.')
        self.logStrings.append('Found ' + str(bounded_fluxes_count) + ' bounded fluxes (with ' + str(closely_bounded_fluxes_count) + ' closely bounded) and ' + str(set_fluxes_count) + ' set fluxes.')



class SBMLExporter(object):
    """
    Helper class to render the given ReactionNetwork as an SBML string on-demand.
    """
    def __init__(self, sourceReactionNetwork=None, exceptionOnError=True):
        self.sourceReactionNetwork=sourceReactionNetwork
        self.exceptionOnError=exceptionOnError
        self.errorCount = 0
        self.logStrings = []


    def getSBMLString(self):
        """
        Generate SBML string from this ReactionNetwork.
        """
        sbmlDoc = self.getSBMLDocument()
        return libsbml.writeSBMLToString(sbmlDoc)


    def writeSBMLToFile(self, fileName):
        """
        Generate SBML from this ReactionNetwork and write it to the given file.
        """
        doc    =  self.getSBMLDocument()
        status = libsbml.writeSBMLToFile(doc,fileName)
        return status
            

    def getSBMLDocument(self):
        """
        Generate SBML document from this ReactionNetwork.
        """
        name = self.sourceReactionNetwork.name
        notes = self.sourceReactionNetwork.notes
        metList = self.sourceReactionNetwork.metList
        reactionList = self.sourceReactionNetwork.reactionList

        # #### Creating model
        nspace = libsbml.SBMLNamespaces(2,1)
        nspace.addNamespace("http://www.w3.org/1999/xhtml","html")
        newDoc = libsbml.SBMLDocument(nspace)
        model = newDoc.createModel()
        model.setId(name)
        model.setName(name)

        # #### Add notes
        notesSBML = convertSBMLNotes(notesDict=notes)
        model.appendNotes(notesSBML)

        # #### Create compartments
        compartmentDict = DB.getCompartmentDict()
        compartmentObjDict = {}
        for c in metList.getCompartments():
            cObj = model.createCompartment()
            cObj.setId(c)
            if c in compartmentDict:
                cObj.setName(compartmentDict[c])
            else:
                cObj.setName(c)
            compartmentObjDict[c] = cObj

        # #### Create Units
        fluxunit = model.createUnitDefinition()
        fluxunit.setId('mmol_per_gDW_per_hr')
        mole = model.createUnit()   
        mole.setKind(libsbml.UNIT_KIND_MOLE)
        mole.setScale(-3)
        gram = model.createUnit()
        gram.setKind(libsbml.UNIT_KIND_GRAM)
        gram.setExponent(-1)
        second = model.createUnit()
        second.setKind(libsbml.UNIT_KIND_SECOND)
        second.setExponent(-1)
        second.setMultiplier(0.00027777)

        # #### Add Species
        for met in metList.mets:
            sp  = model.createSpecies()

            metNameSBML = met.name
            full_name = met.name
            if hasattr(met,'full_name'):
                if met.full_name is not None:
                    full_name = met.full_name

            sp.setId("M_"+metNameSBML)
            sp.setName("M_"+full_name)
            # We might be dealing with old shelved objects that lack a compartment
            if hasattr(met, 'compartment'):
                sp.setCompartment(met.compartment)
            elif isinstance(met.name, core.MetaboliteName):
                sp.setCompartment(met.name.compartment)
            else:
                # There is a slim chance here that the resulting compartment won't be in the compartments list
                # at the top of the document.  We leave it to the SBML library to decide the appropriate response.
                temp_name = core.MetaboliteName(met.name)
                sp.setCompartment(temp_name.compartment)

            spNotes = met.extraNotes
            # Metabolite notes
            if met.destination:
                spNotes['LABELING_DEST'] = met.carbonsOut
            if met.source:
                spNotes['LABELING_SOURCE'] = 'True'
            if met.excluded:
                spNotes['EXCLUDED'] = ''
            if met.source:
                spNotes['FEED'] = met.feed
            spNotes['CARBON NUMBER'] = met.ncarbons
            
            spSBMLNotes = convertSBMLNotes(notesDict=spNotes)
            sp.appendNotes(spSBMLNotes)
            
            # Adding species to model    
            model.addSpecies(sp)
        
        # ### Add Reactions

        for reaction in reactionList.reactions:
            sbmlRxn = model.createReaction()
            SBMLsyntax = "JR904"      #change SBML syntax to the name of the model being used
            if SBMLsyntax == "JR904":
                rxnname = "R_"+reaction.name.replace('(e)','_e_',1)  
            elif SBMLsyntax == "JO1366":
                rxnname = "R_"+reaction.name.replace('(e)','_e',1)    
            elif SBMLsyntax == "iAF1260":
                rxnname = "R_"+reaction.name.replace('(e)','tex',1) 
            sbmlRxn.setId(rxnname)
            sbmlRxn.setName(rxnname)
            sbmlRxn.setReversible(False)
            if reaction.reversible:
                sbmlRxn.setReversible(True)

            #  Reactants
            for react in reaction.reactants:
                stoich = react.stoichiometry
                spref  = sbmlRxn.createReactant()
                spref.setSpecies("M_"+react.name)   # M_ added in front because libsbml screws up if metabolites start with a number
                spref.setStoichiometry(stoich)
            #  Products
            for prod in reaction.products:
                stoich = prod.stoichiometry
                spref  = sbmlRxn.createProduct()
                spref.setSpecies("M_"+prod.name)         
                spref.setStoichiometry(stoich)

            #  Kinetic law (upper and lower bounds)
            kinlaw = sbmlRxn.createKineticLaw()
            parLB = kinlaw.createParameter()
            parLB.setId('LOWER_BOUND')
            parLB.setUnits('mmol_per_gdW_per_hr')
        
            parUB = kinlaw.createParameter()
            parUB.setId('UPPER_BOUND')
            parUB.setUnits('mmol_per_gdW_per_hr')
        
            try:
                UB = reaction.fluxBounds.net.hi    
            except AttributeError:
                UB = core.DEFAULT_UPPER_BOUND

            try:
                LB = reaction.fluxBounds.net.lo
            except AttributeError:
                LB = core.DEFAULT_LOWER_BOUND if reaction.reversible else 0
        
            parLB.setValue(LB)
            parUB.setValue(UB)

            # Kinetic law (cvec = objective vector for FBA)
            parC = kinlaw.createParameter()
            parC.setId('OBJECTIVE_COEFFICIENT')
            objective = 0
            if hasattr(reaction,'cval'):
                objective = reaction.cval
            parC.setValue(objective)

            rNotes = reaction.extraNotes
            # Add subsystem identity, gene and protein association, to the pre-existing notes dictionary
            rNotes['SUBSYSTEM'] = reaction.subsystem
            # NOTE: We can't depend on these attributes being set, because old pickled data may not contain them.
            if hasattr(reaction,'geneSets'):
                if reaction.geneSets:   # Non-empty
                    rNotes['GENE_ASSOCIATION'] = Genes.GeneSet.createStringFromSets(reaction.geneSets)
            if hasattr(reaction,'proteinSets'):
                if reaction.proteinSets:
                    rNotes['PROTEIN_ASSOCIATION'] = Proteins.ProteinSet.createStringFromSets(reaction.proteinSets)
            # Add flux value as note
            if hasattr(reaction,'flux'):
                rNotes['FLUX VALUE'] = reaction.flux
            try:
                transLine  = reaction.transitionLine.fullLine
                # Eliminate reaction name
                newLine = transLine.split('\t')
                newLine = '\t'.join(newLine[1:])   
                # Eliminate > and < because SBML does not like them
                rNotes['CARBON_TRANSITIONS'] = newLine.replace('<==>','==').replace('-->','--')                  
            except AttributeError:    # if no transitions attached to reaction
                pass

            rxnNotesSBML = convertSBMLNotes(notesDict=rNotes)
            sbmlRxn.appendNotes(rxnNotesSBML)

        # Returning sbml file as a string
        errorCount = newDoc.getNumErrors()
        self.errorCount = errorCount
        if errorCount:
            logStrings = []        
            errs = newDoc.getErrorLog()
            for x in range(errorCount):
                err = errs.getError(x)
                logStrings.append(err.getMessage())
            self.logStrings = logStrings
        return newDoc



def convertSBMLNotes(notesSBML='unset', notesDict='unset'):
    """
    Function that translates notes sbml node to Dictionary and vice-versa
    """
    # We use 'unset' instead of None as the deault here because sometimes None is given as a valid input value for notesSBML,
    # and we want to return an empty dictionary for that.

    # The following functions are used to convert to and from the types of records that
    # can be conveyed in an SBML notes object.

    rPassthrough    = lambda line: line
    rToInt          = lambda line: int(line)
    rToStrippedStr  = lambda line: line.strip()
    rToTrue         = lambda line: True

    rFromNumber     = lambda numb: [str(numb)]

    rGetAbbrev      = lambda GCMS: GCMS.abbrev 


    def rFromLinesInDict(fragDict):
        lines = []
        for abbrev in fragDict:
            lines.append(fragDict[abbrev].line)
        return lines


    # looks like e.g. [[0,1],[2,3]]. It means first fragment has carbon 0 and 1 eliminated and second fragment carbons 2 and 3
    def rToDestInfo(line):
        indsnew = []
        inds = line.split(';')
        for ind in inds:
            indsnew.append(ind.strip().split(','))
        return indsnew


    def rFromDestInfo(carbonsOut):
        if carbonsOut is None:
            return []
        newgroups = []
        for group in carbonsOut:
            newgroups.append(','.join(group))
        return [';'.join(newgroups)]


    def rToFlux(line):
        # Examples of acceptable inputs:
        #   <flux forward="240.34" backward="5.0" net="[230:235.34:237.004]" exchange="5.0" />
        # Or:
        #   "240.34, 5.0, [230:235.34:237.004], 5.0"
        # Or:
        #   "240.34, 5.0, 235.34, 5.0"

        if isinstance(line, libsbml.XMLNode):
            forward  = None
            backward = None
            net      = None
            exchange = None
            if line.hasAttr("forward"):
                forward  = line.getAttrValue("forward")
            if line.hasAttr("backward"):
                backward = line.getAttrValue("backward")
            if line.hasAttr("net"):
                net      = line.getAttrValue("net")
            if line.hasAttr("exchange"):
                exchange = line.getAttrValue("exchange")
        else:
            v = line.split(',')
            if len(v) != 4:
                raise Exception('Wrong input for flux value from SBML')
            forward  = v[0]
            backward = v[1]
            net      = v[2]
            exchange = v[3]

        processed = []
        for item in [forward, backward, net, exchange]:
            # If they look like rangedNumbers, turn them into rangedNumber objects
            # before passing them to the flux instantiator.
            if isinstance(item, str):
                item = item.strip()
            if core.rangedNumber.looks_like_rangedNumber(item):
                item = core.rangedNumber.from_string(item)
            else:
                if utils.is_float(item):
                    item = float(item)
            processed.append(item)

        forward  = processed[0]
        backward = processed[1]
        net      = processed[2]
        exchange = processed[3]

        flux = core.flux(for_back_tup=(forward,backward),net_exc_tup=(net,exchange))
        return flux


    def rFromFlux(flux):
        fluxSBML = libsbml.XMLNode()
        
        triple = libsbml.XMLTriple("flux", "", "")

        att = libsbml.XMLAttributes()
        att.add("forward", str(flux.forward)) # Forcing string conversion to invoke __str__ for rangedNumber
        att.add("backward", str(flux.backward))
        att.add("net", str(flux.net))
        att.add("exchange", str(flux.exchange))

        ns = libsbml.XMLNamespaces()
        token = libsbml.XMLToken(triple,att,ns)
        fluxNode = libsbml.XMLNode(token)

        fluxSBML.addChild(fluxNode)
        return [fluxSBML]


    # A dictionary of note headers, providing three things:
    # 1. A function to use for input from SBML line into dictionary key,
    # 2. An inverse function to use for obtaining SBML line from dict element,
    # 3. A function to produce a label for dictionary in the case of multiple entries (e.g. GCMS)
    noteTypes = {'GCMSLabelData'         : (labeling.GCMSLabelData,      rFromLinesInDict,   rGetAbbrev),
                'GCMSLabelData(fit)'     : (labeling.GCMSLabelData,      rFromLinesInDict,   rGetAbbrev),
                'LCMSLabelData'          : (labeling.LCMSLabelData,      rFromLinesInDict,   rGetAbbrev),
                'LCMSLabelData(fit)'     : (labeling.LCMSLabelData,      rFromLinesInDict,   rGetAbbrev),
                'CEMSLabelData'          : (labeling.CEMSLabelData,      rFromLinesInDict,   rGetAbbrev),
                'CEMSLabelData(fit)'     : (labeling.CEMSLabelData,      rFromLinesInDict,   rGetAbbrev),
                'DENS'                   : (rPassthrough,   '',                 ''),
                'CARBON TRANSITIONS TYPE': (rPassthrough,   '',                 ''),
                'CARBON_TRANSITIONS'     : (rToStrippedStr, '',                 ''),
                'LABELING_DEST'          : (rToDestInfo,    rFromDestInfo,      ''),
                'LABELING_SOURCE'        : (rToTrue,        '',                 ''),
                'EXCLUDED'               : (rToTrue,        '',                 ''),
                'FEED'                   : (rPassthrough,   '',                 ''),
                'CARBON NUMBER'          : (rToInt,         rFromNumber,        ''),
                'FLUX VALUE'             : (rToFlux,        rFromFlux,          ''),
                'SUBSYSTEM'              : (rPassthrough,   '',                 ''),     
                'GENE_ASSOCIATION'       : (rPassthrough,   '',                 ''),   
                'PROTEIN_ASSOCIATION'    : (rPassthrough,   '',                 ''),   
                'FORMULA'                : (rPassthrough,   '',                 '')
                }

    # SBML Notes to notesDict
    if notesSBML != 'unset':
        notesDict = {}
        if not notesSBML:
            return notesDict
        notesBody = notesSBML
        if notesSBML.hasChild('body'):
            notesBody = notesSBML.getChild(0)

        for i in range(notesBody.getNumChildren()):
            nElement = notesBody.getChild(i)
            string = nElement.getChild(0).toXMLString()
            try:
                noteHead, noteContent = string.split(':', 1)
            except ValueError:
                continue
            # If there is an additional XML object inside the note element, past the
            # text string, use that as the note content, rather than the remainder of the string.
            if (nElement.getNumChildren() > 1):
                noteContent = nElement.getChild(1);

            header = noteHead.strip()
            if header in noteTypes:
                recordImportFunction, recordExportFunction, labelGeneratorFunc = noteTypes[header]
                # Eliminate header (e.g. GCMS :) from string
                object = recordImportFunction(noteContent)
            else:
                labelGeneratorFunc = ''
                object = noteContent
            # Store
            # This object is meant to return a dictionary, store as a list for the moment
            if labelGeneratorFunc != '' : 
                if header in notesDict:
                    notesDict[header].append(object)                                
                else:
                    notesDict[header] = [object]           
            #    One entry           
            else:
                notesDict[header] = object

        # Convert lists into dictionaries if labelGeneratorFunc has been provided
        for header in notesDict:
            if header in noteTypes:
                recordImportFunction, recordExportFunction, labelGeneratorFunc = noteTypes[header]
                # if there is a function and it is a list (iterable object)
                if labelGeneratorFunc and getattr(notesDict[header],'__iter__',False):
                    dict = {}
                    list = notesDict[header]
                    for object in list:
                        dict[labelGeneratorFunc(object)] = object
                    notesDict[header] = dict

        return notesDict


    # notesDict to SBML Notes    
    elif notesDict != 'unset':
        notesSBML = libsbml.XMLNode()
        
        triple = libsbml.XMLTriple("body", "", "")
        att = libsbml.XMLAttributes()
        ns = libsbml.XMLNamespaces()
        ns.add( "http://www.w3.org/1999/xhtml", "")
        token = libsbml.XMLToken(triple,att,ns)
        notesBodyNode = libsbml.XMLNode(token)

        ns.clear()
        triple = libsbml.XMLTriple("p", "", "")
        token = libsbml.XMLToken(triple,att,ns)
        
        for header in notesDict:
            object = notesDict[header]
            lines = [object]
            if header in noteTypes:
                recordImportFunction, recordExportFunction, labelGeneratorFunc = noteTypes[header]
                if recordExportFunction:
                    lines = recordExportFunction(object)
            # Processing all lines
            for line in lines:
                # If the result of the export function was None, don't bother embedding the note
                if line is None:
                    continue
                node = libsbml.XMLNode(token)
                if isinstance(line, libsbml.XMLNode):
                    tt = libsbml.XMLToken(header + ': ')
                    n1 = libsbml.XMLNode(tt)
                    node.addChild(n1)
                    node.addChild(line)
                else:
                    tt = libsbml.XMLToken(header + ': ' + str(line))
                    n1 = libsbml.XMLNode(tt)
                    node.addChild(n1)
                notesBodyNode.addChild(node)

        notesSBML.addChild(notesBodyNode)
        return notesSBML
    raise Exception('convertSBMLNotes called with no valid input')


##### UNIT TESTS ########

import unittest

# Expected result from tests below
expectedResult =  '153.68'  



class testSBML(unittest.TestCase):  
    """ Testing basic usage of SBML support """

    def setUp(self):
        pass


    def testSBML(self):
        """ Testing SBML support """


    def tearDown(self):
        pass

