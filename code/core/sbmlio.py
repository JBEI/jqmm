# Part of the JBEI Quantitative Metabolic Modeling Library (JQMM)
# Copyright (c) 2016, The Regents of the University of California.
# For licensing details see "license.txt" and "legal.txt".
"""
The sbmlio module provides support for input/output of ReactionNetwork structures in SBML format.
"""

import re
import DB, enhancedLists, labeling
import core
import Genes, Proteins
from SBMLResourceIdentifiers import SBMLResourceIdentifiers
import utilities
import libsbml



class SBMLImporter:
    """
    Helper class to parse the given string as SBML and extract various structures important to a ReactionNetwork.
    """
    def __init__(self, sbmlFileContents=None, exceptionOnError=True):

        self.model_name = None
        self.model_notes = {}

        self.metList      = None
        self.reactionList = None
        self.proteinSet   = None
        self.geneSet      = None

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

        if self.errorCount > 0 and exceptionOnError:
            raise Exception('Errors parsing SBML:\n%s' % ('\n'.join(self.logStrings)))
        else:
            self.parseSBMLDoc(sbmlDoc)



    def parseSBMLDoc(self, sbmlDoc):

        hasFBC = False
        usingFBC = False
        model = sbmlDoc.getModel()

        # Basic check for FBC functionality
        if hasattr(libsbml, 'SBML_FBC_GENEPRODUCT'):
            if libsbml.SBML_FBC_GENEPRODUCT:
                hasFBC = True

        if not hasFBC:
            self.logStrings.append('Warning: libsbml does not appear to have FBC (flux balance constraints) extension compiled in.')

        modelFBCplugin = None
        FBCversion = 0
        FBCstrict = True

        if hasFBC:
            modelFBCplugin = model.getPlugin(str("fbc"))
            if modelFBCplugin is not None:
                usingFBC = True
                FBCversion = modelFBCplugin.getPackageVersion()
                if FBCversion == 2:
                    FBCstrict = modelFBCplugin.getStrict()
                    if not FBCstrict:
                        self.logStrings.append('Warning: This model is using the FBC package in non-strict mode.  Strange things may occur; see the libsbml docs.')

        masterProteinSet = Proteins.ProteinSet()
        masterGeneSet = Genes.GeneSet()

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
        metIdentifiersCount = 0

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

            spFbc = None
            if usingFBC:
                spFbc = sp.getPlugin("fbc")

            if usingFBC and spFbc:
                formula = spFbc.getChemicalFormula()

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

            midentifiers = SBMLResourceIdentifiers()
            if usingFBC and (FBCversion == 2):
                terms = self.readCVTerms(sp)
                if 'is' in terms:
                    metIdentifiersCount += len(terms['is'])
                    midentifiers.add(terms['is'])
                    #self.logStrings.append('  Met identifiers: ' + str(midentifiers))

            met = core.Metabolite(spName, ncarbons, excluded, source, feed, destination, carbonsOut,
                                  full_name=full_name, formula=formula,
                                  extraNotes=spNotes, identifiers=midentifiers)
            mets.append(met)

        self.logStrings.append('Detected %s species, parsed %s metabolites with %s identifiers.' %
                                (str(len(listOfSpecies)), str(len(mets)), str(metIdentifiersCount)))

        # Creating met list    
        metList = enhancedLists.MetaboliteList(mets)

        if (modelFBCplugin is not None) and (FBCversion == 2):
            for gpIndex in range(modelFBCplugin.getNumGeneProducts()):
                oneProduct = modelFBCplugin.getGeneProduct(gpIndex)

                # fbc:geneProduct elements have nultiple sources for names:
                #
                #  fbc:label (example: 'b0584')
                #  fbc:name  (example: 'fepA')
                #  metaid    (example: 'G_b0584')
                #  fbc:id    (example: 'G_b0584')
                #
                # The canonical name is fbc:name, and the name used to resolve against
                # other structures can be metaid or fbc:id, depending on the structure.
                # For our purposes, fbc:id is the relevant one.
                #
                # Some geneProducts do not have an fbc:name, in which case we just use fbc:id.

                pid = oneProduct.getId()
                pid = pid.replace('G_','',1)
                # Since we specify the fbc:name first, it becomes the canonicalName:
                n = oneProduct.getName()
                if n:
                    newProtein = Proteins.Protein([n, pid])
                else:
                    newProtein = Proteins.Protein([pid])
                    newProtein.nonCanonicalName = True

                newProtein.label = oneProduct.getLabel()
                newProtein.sbo = oneProduct.getSBOTermID()
                newGene = Genes.Gene(pid)
                newProtein.encodedBy = [newGene]

                terms = self.readCVTerms(oneProduct)
                if 'is' in terms:
                    newProtein.addIdentifier(terms['is'])
                if 'isEncodedBy' in terms:
                    newGene.addIdentifier(terms['isEncodedBy'])
                masterProteinSet.recastItems(newProtein)
                masterGeneSet.recastItems(newGene)

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
        rxnIdentifiersCount = 0
        for rxn in listOfReactions:
            rxnNameOrig = re.sub('^R_','', rxn.getId())
            rxnNameOrig = re.sub('_e_$','(e)', rxnNameOrig)
            rxnName         = core.ReactionName(rxnNameOrig)
            full_name   = rxn.getName()
            if not isinstance(full_name, str):
                full_name = rxnName
            elif full_name.strip() == '':
                full_name = rxnName

            exchangeReaction = True if rxnName[0:3] == 'EX_' else False 

            if usingFBC:
                reactionFBC = rxn.getPlugin('fbc')

            # Stoichiometry
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
                    
            #  Carbon Transitions and other notes
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
                if 'FLUX_VALUE' in rxnNotes:
                    flux = rxnNotes['FLUX_VALUE']
                if 'FLUX VALUE' in rxnNotes:
                    flux = rxnNotes['FLUX VALUE']
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
            proteinSets = []
            # We will recast the sets we make into identical-looking sets that are composed of
            # already-existing gene/protein objects, if available.
            recastGeneSets = []
            recastProteinSets = []

            # If using FBC package, look for gene product associations in FBC structures.
            if usingFBC:
                gpas = reactionFBC.getGeneProductAssociation()
                if gpas != None:
                    gpasa = gpas.getAssociation()
                    if gpasa != None:
                        self.readGeneProductRefs(gpasa, proteinSets)
                        for subSet in proteinSets:
                            recastProteinSets.append(masterProteinSet.recastSet(subSet, mergeValues=True))
                        # Use the 'encodedBy' information read from the global set
                        # (if it exists) to build a list of correlated genes
                        for ps in recastProteinSets:
                            gs = []
                            for p in ps.itemList():
                                for gn in p.encodedBy:
                                    gs.append(gn)
                            if len(gs) > 0:
                                geneSets.append(Genes.GeneSet(gs))

            if geneString:
                geneString = re.sub('^\s*:','', geneString) # Chop off any bogus colon
                geneString = geneString.strip()
                if geneString.lower() in ["none", "n.a.", "n/a"]:
                    geneString = ""
                if geneString != "":
                    geneSets = Genes.GeneSet.createSetsFromStrings(geneString, geneValueString)
            if proteinString:
                proteinString = re.sub('^\s*:','', proteinString) # Chop off any bogus colon
                proteinString = proteinString.strip()
                if proteinString.lower() in ["none", "n.a.", "n/a"]:
                    proteinString = ""
                if proteinString != "":
                    # The gene names act as alternate names fro the proteins here.
                    # (We assume they are both specified in the same order, and link them up.)
                    # This is done because at least one SBML exporter (it is not known which) erroneously pair the
                    # transcription value with the protein identifier.
                    proteinSets = Proteins.ProteinSet.createSetsFromStrings(proteinString, geneString, proteinValueString)
                    for subSet in proteinSets:
                        recastProteinSets.append(masterProteinSet.recastSet(subSet, mergeValues=True))

            for subSet in geneSets:
                recastGeneSets.append(masterGeneSet.recastSet(subSet, mergeValues=True))

            #  Measured fluxes (MeasFluxes for 13CMFA, UB and LB for 2S-13CMFA or FBA)
            cval = 0
            LB = None
            UB = None

            # If using FBC package, attempt to locate and then do some rudimentary translation
            # on the upper and lower bounds.
            if usingFBC:
                LB = reactionFBC.getLowerFluxBound()
                UB = reactionFBC.getUpperFluxBound()
                bounds_translation = {
                    'cobra_0_bound': '0',
                    'cobra_default_lb': core.DEFAULT_LOWER_BOUND,
                    'cobra_default_ub': core.DEFAULT_UPPER_BOUND
                }
                if LB is not None:
                    if LB in bounds_translation:
                        LB = bounds_translation[LB]
                    # TODO: Right now if it's not a number we don't know what to do with it, so we set it to 0.
                    # This is bad.  It could be something like 'R_ATPM_lower_bound', 'R_EX_cbl1_e_lower_bound', or
                    # anything else defined in the <listOfParameters> section of the SBML document...
                    elif not utilities.is_float(LB):
                        LB = core.DEFAULT_LOWER_BOUND

                if UB is not None:
                    if UB in bounds_translation:
                        UB = bounds_translation[UB]
                    elif not utilities.is_float(UB):
                        LB = core.DEFAULT_UPPER_BOUND

            # If we have a kinetic law annotated with old-style flux bounds, we'll prefer that over the FBC package.
            kLaw = rxn.getKineticLaw()
            if kLaw:
                listOfParameters = kLaw.getListOfParameters()
                for par in listOfParameters:
                    if par.getId() == 'LOWER_BOUND':
                        LB = par.getValue()
                    elif par.getId() == 'UPPER_BOUND':
                        UB = par.getValue()
                    elif par.getId() == 'FLUX_VALUE':
                        if flux is None:    # Prefer the one from the notes
                            flux = core.flux(for_back_tup=('NA','NA'),net_exc_tup=(par.getValue(),'NA'))
                            # Sigh...  Yet another place where we're going to have to remove 'NA' and replace it with something
                            # sane...
                    elif par.getId() == 'OBJECTIVE_COEFFICIENT':
                        cval = par.getValue()
                        Cvec[rxnName] = par.getValue()

            if (LB != None) and (UB != None):
                bounded_fluxes_count += 1
                if (LB == UB):
                    closely_bounded_fluxes_count += 1

            measured,fluxBounds = core.fluxBounds(LB, UB, rxn.getReversible(), exchangeReaction)

            if flux is not None:
                set_fluxes_count += 1

            # If we don't have a determined flux, use the flux bounds.
            if not flux:
                flux = fluxBounds

            ridentifiers = SBMLResourceIdentifiers()
            if usingFBC and (FBCversion == 2):
                terms = self.readCVTerms(rxn)
                if 'is' in terms:
                    rxnIdentifiersCount += len(terms['is'])
                    ridentifiers.add(terms['is'])
                    #self.logStrings.append('  React identifiers: ' + str(ridentifiers))

            # Instantiating reaction
            new_reaction = core.Reaction(rxnName, rxn.getReversible(), False, measured, flux=flux, fluxBounds=fluxBounds,
                                     transitionLine=transitionLine, exchange=exchangeReaction, reactants=reacts, products=prods,
                                     cval=cval, subsystem=subsystem, gene=recastGeneSets, protein=recastProteinSets,
                                     extraNotes=rxnNotes, identifiers=ridentifiers, full_name=full_name)
            reactions.append(new_reaction)

        self.metList      = metList
        self.reactionList = enhancedLists.ReactionList(reactions)
        self.proteinSet   = masterProteinSet
        self.geneSet      = masterGeneSet

        self.initialSetFluxesCount     = set_fluxes_count
        self.initialBoundedFluxesCount = bounded_fluxes_count
        self.initialCloselyBoundedFluxesCount = closely_bounded_fluxes_count

        self.logStrings.append('Detected %s reactions, parsed %s, with %s identifiers.' %
                                (str(len(listOfReactions)), str(len(reactions)), str(rxnIdentifiersCount)))
        self.logStrings.append('Found %s bounded fluxes (with %s closely bounded) and %s set fluxes.' %
                                (str(bounded_fluxes_count), str(closely_bounded_fluxes_count), str(set_fluxes_count)))
        self.logStrings.append('Detected %s genes, %s proteins.' %
                                (str(len(masterGeneSet)), str(len(masterProteinSet))))

    qualifiersToElements = {
        libsbml.BQB_ENCODES         : "encodes",        # 8
        libsbml.BQB_HAS_PART        : "hasPart",        # 1
        libsbml.BQB_HAS_PROPERTY    : "hasProperty",    # 10
        libsbml.BQB_HAS_VERSION     : "hasVersion",     # 4
        libsbml.BQB_IS              : "is",             # 0
        libsbml.BQB_IS_DESCRIBED_BY : "isDescribedBy",  # 6
        libsbml.BQB_IS_ENCODED_BY   : "isEncodedBy",    # 7
        libsbml.BQB_IS_HOMOLOG_TO   : "isHomologTo",    # 5
        libsbml.BQB_IS_PART_OF      : "isPartOf",       # 2
        libsbml.BQB_IS_PROPERTY_OF  : "isPropertyOf",   # 11
        libsbml.BQB_IS_VERSION_OF   : "isVersionOf",    # 3
        libsbml.BQB_OCCURS_IN       : "occursIn",       # 9
        libsbml.BQB_UNKNOWN         : "unknown",        # 12

#        libsbml.BQM_IS              : "isB",             # 0   # Should be in another structure eventually
#        libsbml.BQM_IS_DERIVED_FROM : "isDerivedFrom",   # 2
#        libsbml.BQM_IS_DESCRIBED_BY : "isDescribedByB",  # 1
#        libsbml.BQM_UNKNOWN         : "unknownB"         # 3
    }


    def readCVTerms(self, sbaseObj):
        terms = {}
        if sbaseObj.getNumCVTerms() < 1:
            return terms
        for i in range(sbaseObj.getNumCVTerms()):
            oneTerm = sbaseObj.getCVTerm(i)
            if oneTerm.getQualifierType() == libsbml.BIOLOGICAL_QUALIFIER:
                name = SBMLImporter.qualifiersToElements[oneTerm.getBiologicalQualifierType()]
            else:
                name = SBMLImporter.qualifiersToElements[oneTerm.getModelQualifierType()]
            resourceList = []
            for j in range(oneTerm.getNumResources()):
                uri = oneTerm.getResourceURI(j)
                if name != None and uri != None:
                    resourceList.append(uri)
            terms[name] = tuple(resourceList)
        return terms


    # This function calls itself recursively to walk through sets of sets when
    # indicated by the ibsbml.FbcOr element.
    def readGeneProductRefs(self, gpasa, listContainer):
        if isinstance(gpasa, libsbml.GeneProductRef):
            p = gpasa.getGeneProduct()
            listContainer.append(Proteins.ProteinSet(Proteins.Protein(p.replace('G_','',1))))
        elif isinstance(gpasa, libsbml.FbcAnd):
            subList = []
            for i in range(gpasa.getNumAssociations()):
                gpasaSub = gpasa.getAssociation(i)
                p = gpasaSub.getGeneProduct()
                subList.append(Proteins.Protein(p.replace('G_','',1)))
            listContainer.append(Proteins.ProteinSet(subList))
        elif isinstance(gpasa, libsbml.FbcOr):
            for i in range(gpasa.getNumAssociations()):
                gpasaSub = gpasa.getAssociation(i)
                self.readGeneProductRefs(gpasaSub, listContainer)
        else:
            self.logStrings.append('Strange element structure found while parsing fbc geneProduct associations!')



class SBMLExporter:
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
                if utilities.is_float(item):
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



def test():  
    """ Testing basic usage of SBML support """
    sbmlFileName = "SBMLTestFile.sbml"
    sbmlFileContents = None
    with open(sbmlFileName, 'r') as file:
        sbmlFileContents = file.read()
    file.closed
    assert sbmlFileContents is not None, "Could not load SBML test file %s." % sbmlFileName

    sbmlImporter = SBMLImporter(sbmlFileContents)

    name         = sbmlImporter.model_name
    notes        = sbmlImporter.model_notes
    metList      = sbmlImporter.metList
    reactionList = sbmlImporter.reactionList

    errorCount   = sbmlImporter.errorCount
    logStrings   = sbmlImporter.logStrings

    if errorCount > 0 and exceptionOnError:
        raise Exception('Errors interpreting SBML:\n%s' % ('\n'.join(logStrings)))
    else:
        print '\n'.join(logStrings)


if __name__ == "__main__":
    test()

