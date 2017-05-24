# Part of the JBEI Quantitative Metabolic Modeling Library (JQMM)
# Copyright (c) 2016, The Regents of the University of California.
# For licensing details see "license.txt" and "legal.txt".

import re
import core
import NamedRangedNumber
from SBMLResourceIdentifiers import SBMLResourceIdentifiers


class Protein(NamedRangedNumber.NamedRangedNumber):
    """
    Class for single proteins, and values typically associated with them.
    Typically it is instantiated with a string representing a name, and a value.
       Since proteins can potentially have multiple names due to conflicting standards, the superclass also supports
    receiving a list of names during instantiation, instead of a string.
    The first name in the list will be considered the canonical name when rendering the protein as a string.
       The given value can be an integer, a float, a ranged number, or a string representation of any of these,
    but is kept internally and exported as a ranged number.
    """
    def __init__(self, names, value=None):
        if isinstance(names, list):
            nameList = names
        else:
            nameList = [names]
        for name in nameList:
            if ' ' in name.strip():
                raise ValueError("Protein names cannot contain spaces: '%s'" % name)

        self.label = None
        self.sbo = None
        self.identifiers = SBMLResourceIdentifiers()
        self.encodedBy = []
        # For use when this is an unknown/unnamed protein linking a gene to a reaction.
        # The name is copied from the gene in that case, and this flag is set.
        self.nonCanonicalName = False
        super(Protein, self).__init__(names, value)


    def addName(self, name):
        if ' ' in name.strip():
            raise ValueError("Protein names cannot contain spaces: '%s'" % name)
        super(Protein, self).addName(name)


    def addIdentifier(self, urls):
        self.identifiers.add(urls)


    def mergeValues(self,other):
        "Merge the contents of other into this instance."
        super(Protein, self).mergeValues(other)
        if self.label is None:
            self.label = other.label
        if self.sbo is None:
            self.sbo = other.sbo
        if self.identifiers is None:
            self.identifiers = other.identifiers
        # The sensible behavior for an encodedBy set would be to merge it.
        # We could attempt this, but it would require making assumptions about the Gene objects
        # within that we cannot make.  So instead, we'll just check to see if there is NO set in the
        # target Protein, and carry the old set along if so.
        # (FWIW, this solves our immediate issue with the old vs. new SBML import methods.)
        if len(self.encodedBy) < 1:
            self.encodedBy = other.encodedBy


class ProteinSet(NamedRangedNumber.NamedRangedNumberSet):
    """
    Class for a set of Protein objects, derived from NamedRangedNumberSet.
    """
    def __init__(self, contents=None):
        super(ProteinSet, self).__init__(contents)


    def recastSet(self, victimSet, preferExistingObjects=True, mergeValues=False):
        itemsToRecast = victimSet.contentsList
        recastItems = self.recastItems(itemsToRecast, preferExistingObjects, mergeValues)
        return ProteinSet(recastItems)



    @staticmethod
    def createSetsFromStrings(structureString="", secondaryStructureString=None, valueString=None):
        """
        Construct a list of ProteinSet objects based on the given strings.
          The first string, structureString, determines the number of ProteinSets to create.
          An empty string will return an empty list.
          A string containing anything else will return one or more subsets.
          The number of subsets is determined by the number of times the separator string " or " occurs in the
        structureString.  For example, "(setA) or (setB) or (setC) or proteinD or proteinE or (setF)" will create
        six subsets.  "(setA)", "(setC)", "proteinD", etc are substrings that declare the contents of each set.
          There are two accepted ways to format the substrings:
        Method #1:
            substrings example: "aName=[aLow:aBest:aHigh] and bName=[bLow:bBest:bHigh] and cName=value and dName"
            valueString=None
          In this example, four NamedRangedNumber objects will be created in total:
        aName and bName specify rangedNumber values for NamedRangedNumber, cName specifies just one floating point
        number that is converted to a rangedNumber, and dName crates a NamedRangedNumber with the value set to None.
        valueString can be left out entirely.
        Method #2:
            substrings example: "aName and dName and cName and qName"
            valueString example: "aName=[aLow:aBest:aHigh] bName=[bLow:bBest:bHigh] cName=value dName=value fName=value"
          In this example, four NamedRangedNumber objects will be created, but only two of them will be assigned values
        (the other two will have values of None).  This happens because the structureString declares what items are in
        the set, while the valueString only assigns values.  If a value is given in the second string for a name that is
        not listed in the first, that value is ignored.  No item is created for it.
          While it is possible to supply a mixture of methods 1 and 2, it is not recommended practice.  Values assigned via
        method 1 take precedence over values assigned via method 2, even if the value assigned is "=None".
          Note that Protein objects are re-used from one set to the next.  That is, if the same name is mentioned in two
        different substrings, only one Protein object will be created but it will be placed in two subsets.
          secondaryStructureString:
          This string, if supplied, is meant to supply alternate names for the names given in structureString.
        It should have the same structure - including 'and' and 'or' components and parenthesis - so it can provide a
        1-to-1 correspondence when walked along in linear fashion.
        """

        givenValues = {}
        if valueString is not None:
            pairs = valueString.split()    # Split on whitespace, no need to strip
            for pair in pairs:
                parts = pair.split('=')
                name = parts[0]
                if parts[1:2]:
                    givenValues[name] = parts[1]

        alternateNamesInOrder = []
        if secondaryStructureString:
            if secondaryStructureString != "":
                a = secondaryStructureString.replace('(',' ').replace(')',' ').strip()
                aParts = re.split("\s+or\s+", a)
                for b in aParts:
                    bParts = re.split("\s+and\s+", b)
                    for c in bParts:
                        alternateNamesInOrder.append(c)

        subSets = []
        structureString = structureString.strip() # Stripping initial surrounding whitespace in order to check for a blank entry
        if structureString != "":
            collectionStrings = re.split("\s+or\s+", structureString)
            for collectionStr in collectionStrings:
                items = []

                # Sections of the string are sometimes enclosed in parenthesis.
                # Plus, sometime garbage from badly-done embeds comes up, like so:
                #  <html:p> GENE_ASSOCIATION :( b1901  and  b1900  and  ( b1898  and  b1899 ) ) </html:p>
                collectionStr = collectionStr.replace('(',' ').replace(')',' ').strip()

                itemStrings = re.split("\s+and\s+", collectionStr)
                for itemString in itemStrings:
                    item = Protein.fromString(itemString)
                    if alternateNamesInOrder[0:1]:
                        item.addName(alternateNamesInOrder[0])
                        alternateNamesInOrder = alternateNamesInOrder[1:]
                    for n in item.names:
                        if givenValues.has_key(n):
                            item.set(givenValues[n])
                            break
                    items.append(item)
                if items:
                    subSets.append(ProteinSet(items))
        return subSets



def test():
    try:
        vError = False
        print "Instantiating from illegal string \"test=[2,3,4]\", expecting failure ..."
        a = Protein.fromString("test=[2,3,4]")
    except ValueError:
        vError = True
        print "\tGot ValueError as expected."
        pass
    assert vError, "NamedRangedNumber accepted wrong input."

    print "\nInstantiating from string value \"test=[2:3:4]\" ..."

    a = Protein.fromString("test=[2:3:4]")
    assert a.canonicalName == 'test', "Name wrong"
    b = Protein.fromString("dorks=[0.5:1:1.5]")
    c = a + 3
    print "\t" + str(a) + ' + 3 = ' + str(c)
    d = a + b
    print "\t" + str(a) + ' + ' + str(b) + ' = ' + str(d)
    assert d.value.best == 4.0, "Addition failure, d.value.best should be 4.0."

    print "\nInstantiating a ProteinSet from an invalid string, expecting failure:"
    strA = "(bob fred frank) or (jed and Bill123) and (fred & billyBob) or captainAmerica"
    print "\t" + strA
    try:
        aError = False
        proteinSets = ProteinSet.createSetsFromStrings(strA)
    except ValueError:
        aError = True
        print "\tGot ValueError as expected."
        pass
    assert aError, "ProteinSet.createSetsFromStrings accepted wrong input."

    print "\nInstantiating a ProteinSet from strings:"
    strA = "(bob and fred and frank) or (jed and Bill123) or (fred and billyBob) or captainAmerica"
    strB = "bob=12 fred=45 frank=[1:22:32] jed=13.1 Bill123=100"
    print "\t" + strA
    print "\t" + strB
    subSets = ProteinSet.createSetsFromStrings(strA, None, strB)

    masterSet = ProteinSet()
    newSubSets = []
    print "Master set:"
    for subSet in subSets:
        newSubSets.append(masterSet.recastSet(subSet))
    print "\t" + str(masterSet)
    print "Subsets consolidated, for embedding:"
    print "\t" + ProteinSet.createStringFromSets(newSubSets)

    print "Significance tests against master set:"
    print "\tAgainst 12:" + str(masterSet.testSignificance(12))
    assert masterSet.testSignificance(12), "Significance test fails."
    print "\tAgainst 12 with equals not counting:" + str(masterSet.testSignificance(12, equalIsntEnough=True))
    assert not masterSet.testSignificance(12, equalIsntEnough=True), "Significance test fails."
    print "\tAgainst 14:" + str(masterSet.testSignificance(14))
    assert not masterSet.testSignificance(14), "Significance test fails."

    print "\nInstantiating a ProteinSet from odd input:"
    strA = "(AraF and AraG and AraH)"
    strB = "( b1901  and  b1900  and  ( b1898  and  b1899 ) )"
    strC = "AraG=15 b1898=20"
    print "\tMain string: " + strA
    print "\tSecondary string: " + strB
    print "\tValue string: " + strC
    subSets = ProteinSet.createSetsFromStrings(strA, strB, strC)

    masterSet = ProteinSet()
    newSubSets = []
    for subSet in subSets:
        newSubSets.append(masterSet.recastSet(subSet))
    print "Master set:"
    print "\t" + str(masterSet)
    print "Subsets consolidated, for embedding:"
    print "\t" + ProteinSet.createStringFromSets(newSubSets)
    print "Adding identifiers:"
    a.addIdentifier([ "http://identifiers.org/ncbigi/GI:16130343",
            "http://identifiers.org/ncbigene/946880",
            "http://identifiers.org/asap/ABE-0007971"])
    a.addIdentifier("http://identifiers.org/ecogene/EG10165")
    print "Expecting to find 'ncbigi' identifier with 'GI:16130343':"
    print a.identifiers.dict
    assert 'ncbigi' in a.identifiers.dict, "Key not found."
    assert a.identifiers.dict['ncbigi'] == 'GI:16130343', "Value is incorrect."


if __name__ == "__main__":
    test()

