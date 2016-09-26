# Part of the JBEI Quantitative Metabolic Modeling Library (JQMM)
# Copyright (c) 2016, The Regents of the University of California.
# For licensing details see "license.txt" and "legal.txt".

import re
import core


class NamedRangedNumber(object):
    """
    Class for a rangedNumber, and a set of names associated with it.
       Typically it is instantiated with one string representing a name, but it also supports
    receiving a list of names.  The first name in the list will be considered the canonical name when rendering
    the NamedRangedNumber as a string.
       The given value can be an integer, a float, a rangedNumber, or a string representation of any of these,
    but is kept internally and exported as a rangedNumber.  A special exception is made for the value None,
    which creates a NamedRangedNumber with no rangedNumber component.  When rendering as a string, it will
    only print the name.  Mathematical operations against None will raise an exception.
    """

    def __init__(self, names, value=None):

        if isinstance(names, str):     # Single name
            self.canonicalName = names
            self.names = [names]
        elif isinstance(names, list):  # List of names
            self.names = []
            for n in names:
                if n not in self.names:
                    self.names.append(n)
            self.canonicalName = names[0]
        else:
            raise Exception('Argument 1 must be a string or a list of strings')
        self.set(value)


    def addName(self, name):
        if isinstance(name, str):
            if name not in self.names:
                self.names.append(name)
        else:
            raise Exception('Name to add must be a string.')


    def set(self, value=None):
        if isinstance(value, list):   # Strange that someone would hand us a list, but we might as well deal with it
            value = value[0]
        if value is None:
            self.value = None
        elif isinstance(value, NamedRangedNumber):
            self.value = value.value.copy()   # Clone the RangedNumber within the NamedRangedNumber
        elif isinstance(value, core.rangedNumber):
            self.value = value.copy()   # Clone the RangedNumber
        elif isinstance(value, str):
            self.value = core.rangedNumber.from_string(value)
        else:
            self.value = core.rangedNumber.from_number(value)


    @classmethod
    def fromString(cls, string):
        """
        Construct a new NamedRangedNumber based on a string supplied in the same format as "__str__"
        Format is "name=[low:best:high]" or "name=value" or simply "name".
        In the third case, the value of the NamedRangedNumber is set to 'None'.
        """
        parts = re.split("\s*=\s*", string)
        name = parts[0].strip()

        if not parts[1:2]:
            return cls(name, None)
        return cls(name, parts[1].strip())


    def __str__(self):
        if self.value is None:
            return self.canonicalName
        return self.canonicalName + '=' + str(self.value)   # RangedNumber __str__ method will take care of the rest


    def __add__(self, other):
        "NamedRangedNumber value addition."
        if isinstance(other, NamedRangedNumber):    # (Note: works for derived classes as well, like Protein)
            v = self.value + other.value
        else:
            # We're assuming that rangedNumber can handle any other input we pass it
            v = self.value + other
        return NamedRangedNumber(self.names, v)

    def __sub__(self, other):
        "NamedRangedNumber value subtraction."
        if isinstance(other, NamedRangedNumber):
            v = self.value - other.value
        else:
            v = self.value - other
        return NamedRangedNumber(self.names, v)
        
    
    def __mul__(self,other):
        "NamedRangedNumber value multiplication."
        if isinstance(other, NamedRangedNumber):
            v = self.value * other.value
        else:
            v = self.value * other
        return NamedRangedNumber(self.names, v)


    def __rmul__(self,other):
        "NamedRangedNumber value reverse multiplication."
        if isinstance(other, NamedRangedNumber):
            v = other.value * self.value
        else:
            v = other * self.value
        return NamedRangedNumber(self.names, v)


    def __div__(self,other):
        "NamedRangedNumber value division."
        if isinstance(other, NamedRangedNumber):
            v = self.value / other.value
        else:
            v = self.value / other
        return NamedRangedNumber(self.names, v)


    def __eq__(self,other):
        "NamedRangedNumber value equivalency."
        if isinstance(other, NamedRangedNumber):
            if self.value != other.value:
                return False
        else:
            if self.value != other:
                return False
        return True

    def add(self, other):
        "NamedRangedNumber value addition."
        if isinstance(other, NamedRangedNumber):
            self.value = self.value + other.value
        else:
            # We're assuming that rangedNumber can handle any other input we pass it
            self.value = self.value + other

    def sub(self, other):
        "NamedRangedNumber value subtraction."
        if isinstance(other, NamedRangedNumber):
            self.value = self.value - other.value
        else:
            self.value = self.value - other
        
    
    def mul(self,other):
        "NamedRangedNumber value multiplication."
        if isinstance(other, NamedRangedNumber):
            self.value = self.value * other.value
        else:
            self.value = self.value * other


    def div(self,other):
        "NamedRangedNumber value division."
        if isinstance(other, NamedRangedNumber):
            self.value = self.value / other.value
        else:
            self.value = self.value / other



class NamedRangedNumberSet(object):
    """
    Class for a set of NamedRangedNumber objects, with routines to deal with duplicates and generate a
    string suitable for printing and embedding in SBML notes.
    """

    def __init__(self, contents=None):
        self.contentsList = []
        self.contentsDict = {}
        self.recastItems(contents)


    def has(self, name):
        if name in self.contentsDict:
            return True
        else:
            return False


    def get(self, name):
        if name in self.contentsDict:
            return self.contentsDict[name]
        else:
            return None


    def itemList(self):
        return self.contentsList


    def isEmpty(self):
        if len(self.contentsList) > 0:
            return False
        return True


    def recastItems(self, contents, preferExistingObjects=True, preferExistingValues=False):
        """
        This method is used to add a list of items to the set, and to remake the given list of items using the merged
        contents of the set and the list. Some examples will help.
          Say the set contains four objects, with names A, B, C, and D.  Now we call this method, and
        supply a list of three objects, with names A, C, and E.  We leave the flags set to their defaults.
          This method will update the _value_ in the objects named A and C in the set with the _value_ in the
        objects named A and C from the list, add object E to the set, and then return a new list with objects
        A, C, and E from the set.  The A and C objects given in the list will be discarded.  We have
        effectively recreated the list, with names and values, using objects in the set, and updated the values
        in the set.
          This is very useful for taking many sets of different objects and getting them all re-made using
        one master set, while still preserving their names and values.
        """
        addedOrFoundItems = []
        if contents is None:
            return addedOrFoundItems
        if isinstance(contents, NamedRangedNumber):
            contents = [contents]
        # If it's not a list at this point, bail.
        if not isinstance(contents, list):
            return addedOrFoundItems
        itemsToRemove = []
        itemsToAdd = []
        for i in contents:
            collision = False
            for n in i.names:
                if n in self.contentsDict:
                    collision = self.contentsDict[n]
                    break
            if not collision:
                for n in i.names:
                    self.contentsDict[n] = i
                itemsToAdd.append(i)
                addedOrFoundItems.append(i)
            else:
                mergedNames = {}
                for n in i.names:
                    mergedNames[n] = 1
                for n in collision.names:
                    mergedNames[n] = 1
                mergedNameList = mergedNames.keys()
                if preferExistingObjects:
                    if not preferExistingValues:
                        collision.value = i.value
                    # Make sure that the currently existing NamedRangedNumber can be referenced by
                    # all the names of the colliding one as well.
                    for n in i.names:
                        self.contentsDict[n] = collision
                    addedOrFoundItems.append(collision)
                else:
                    if preferExistingValues:
                        i.value = collision.value
                    # Make sure that the colliding NamedRangedNumber can be referenced by
                    # all its names as well as the names of the one currently existing.
                    for n in mergedNameList:
                        self.contentsDict[n] = i
                    # It's debatable wether we should be updating the names kept by the item itself
                    # to the merged set of names, but for now, we are doing so.
                    i.names = mergedNameList
                    addedOrFoundItems.append(i)
                    itemsToAdd.append(i)
                    itemsToRemove.append(collision)
        if itemsToRemove:
            # If there are items we want to remove from the list, the most efficient way is to
            # build a hash table for those items and then filter the list by the hash table all at once.
            # (It's a LOT faster than walking though the whole list to dump out a single object, for every object.)
            canonicalNameHash = {}  # We'll build a hash by canonical name, since each item is guaranteed to have one.
            for i in itemsToRemove:
                canonicalNameHash[i.canonicalName] = i
            oldList = self.contentsList
            self.contentsList = []
            for i in oldList:
                # If the 'remove' dictionary has an item under the same name,
                # and that item is the same object (by reference) as the current one,
                # don't add the current one to the new list, effectively filtering it out.
                if i.canonicalName in canonicalNameHash:
                    if i is canonicalNameHash[i.canonicalName]:
                        continue
                self.contentsList.append(i)
        # Only add new items after we've removed the old ones, to make the removal process faster.
        for i in itemsToAdd:
            self.contentsList.append(i)

        return addedOrFoundItems


    def testSignificance(self, threshold, equalIsntEnough=False):
        """
        Test all the ranged numbers in the set to see if the set is 'significant'.
        If all the ranged numbers are defined (not 'None') and ALL their 'best' values
        are equal to or above the threshold, the set is considered significant.
        """
        sectionHasASignificantValue = False
        sectionHasAnInsignificantValue = False
        # We need at least one value, and we need all of them to be significant,
        # and none of them to be insignificant.
        for item in self.contentsList:
            if item.value is not None:
                if item.value.best == threshold:
                    if equalIsntEnough:
                        sectionHasAnInsignificantValue = True
                    else:
                        sectionHasASignificantValue = True
                elif item.value.best < threshold:
                    sectionHasAnInsignificantValue = True
                else:
                    sectionHasASignificantValue = True                            
        if sectionHasASignificantValue and not sectionHasAnInsignificantValue:
            return True
        return False


    @staticmethod
    def createStringFromSets(sets):
        nonEmptySets = [i for i in sets if not i.isEmpty()]
        return ' or '.join(map(str, nonEmptySets))


    def __str__(self):
        if len(self.contentsList) > 1:
            return '(' + ' and '.join(map(str, self.contentsList)) + ')'
        elif len(self.contentsList) == 1:
            return str(self.contentsList[0])
        else:
            return ''



if __name__ == "__main__":

    print "TODO: Testing here..."

