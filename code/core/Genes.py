from __future__ import print_function
# Part of the JBEI Quantitative Metabolic Modeling Library (JQMM)
# Copyright (c) 2016, The Regents of the University of California.
# For licensing details see "license.txt" and "legal.txt".

from builtins import str
import re
import core
import NamedRangedNumber


class Gene(NamedRangedNumber.NamedRangedNumber):
    """
    Class for single genes, and values typically associated with them.
    Typically it is instantiated with a string representing a name, and a value.
       Since genes can potentially have multiple names due to conflicting standards, the superclass also supports
    receiving a list of names during instantiation, instead of a string.
    The first name in the list will be considered the canonical name when rendering the gene as a string.
       The given value can be an integer, a float, a ranged number, or a string representation of any of these,
    but is kept internally and exported as a ranged number.
    """
    def __init__(self, names, value=None):
        if isinstance(names, list):
            nameList = names
        else:
            nameList = [names]

        for name in nameList:
            assert ' ' not in name.strip(), "Gene names cannot contain spaces: '" + name + "'"
        super(Gene, self).__init__(names, value)


    def addName(self, name):
        assert ' ' not in name.strip(), "Gene names cannot contain spaces: '" + name + "'"
        super(Gene, self).addName(name)



class GeneSet(NamedRangedNumber.NamedRangedNumberSet):
    """
    Class for a set of GeneSet objects, derived from NamedRangedNumberSet.
    """
    def __init__(self, contents=None):
        super(GeneSet, self).__init__(contents)


    def recastSet(self, victimSet, preferExistingObjects=True, preferExistingValues=False):
        itemsToRecast = victimSet.contentsList
        recastItems = self.recastItems(itemsToRecast, preferExistingObjects, preferExistingValues)
        return GeneSet(recastItems)



    @staticmethod
    def createSetsFromStrings(structureString="", valueString=None):
        """
        Construct a list of GeneSet objects based on the given strings.
          The first string, structureString, determines the number of GeneSets to create.
          An empty string will return an empty list.
          A string containing anything else will return one or more subsets.
          The number of subsets is determined by the number of times the separator string " or " occurs in the
        structureString.  For example, "(setA) or (setB) or (setC) or geneD or geneE or (setF)" will create
        six subsets.  "(setA)", "(setC)", "geneD", etc are substrings that declare the contents of each set.
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
        method 2 take precedence over values assigned via method 1, even if the value assigned is "=None".
          Note that Gene objects are re-used from one set to the next.  That is, if the same name is mentioned in two
        different substrings, only one Gene object will be created but it will be placed in two subsets.
        """

        givenValues = {}
        if valueString is not None:
            pairs = valueString.split()    # Split on whitespace, no need to strip
            for pair in pairs:
                parts = pair.split('=')
                name = parts[0]
                if parts[1:2]:
                    givenValues[name] = parts[1]

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
                    item = Gene.fromString(itemString)
                    if item.canonicalName in givenValues:
                        item.set(givenValues[item.canonicalName])
                    items.append(item)
                if items:
                    subSets.append(GeneSet(items))
        return subSets



if __name__ == "__main__":
    test()

def test():

    try:
        vError = False
        print("Instantiating from illegal string \"test=[2,3,4]\", expecting failure ...")
        a = Gene.fromString("test=[2,3,4]")
    except ValueError:
        vError = True
        print("\tGot ValueError as expected.")
        pass
    assert vError, "NamedRangedNumber accepted wrong input."

    print("\nInstantiating from string value \"test=[2:3:4]\" ...")

    a = Gene.fromString("test=[2:3:4]")
    assert a.canonicalName == 'test', "Name wrong"
    b = Gene.fromString("dorks=[0.5:1:1.5]")
    c = a + 3
    print("\t" + str(a) + ' + 3 = ' + str(c))
    d = a + b
    print("\t" + str(a) + ' + ' + str(b) + ' = ' + str(d))
    assert d.value.best == 4.0, "Addition failure, d.value.best should be 4.0."

    print("\nInstantiating a GeneSet from an invalid string, expecting failure:")
    strA = "(bob fred frank) or (jed and Bill123) and (fred & billyBob) or captainAmerica"
    print("\t" + strA)
    try:
        aError = False
        geneSets = GeneSet.createSetsFromStrings(strA)
    except AssertionError:
        aError = True
        print("\tGot AssertionError as expected.")
        pass
    assert aError, "GeneSet.createSetsFromStrings accepted wrong input."

    print("\nInstantiating a GeneSet from strings:")
    strA = "(bob and fred and frank) or (jed and Bill123) or (fred and billyBob) or captainAmerica"
    strB = "bob=12 fred=45 frank=[1:2:3] jed=10.1 Bill123=1"
    print("\t" + strA)
    print("\t" + strB)
    subSets = GeneSet.createSetsFromStrings(strA, strB)
    masterSet = GeneSet()
    newSubSets = []
    print("Master set:")
    for subSet in subSets:
        newSubSets.append(masterSet.recastSet(subSet))
    print("\t" + str(masterSet))
    print("Subsets consolidated, for embedding:")
    print("\t" + GeneSet.createStringFromSets(newSubSets))
    print("Significance test result for master set:" + str(masterSet.testSignificance(12)))


