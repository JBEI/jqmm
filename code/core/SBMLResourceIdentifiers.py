# Part of the JBEI Quantitative Metabolic Modeling Library (JQMM)
# Copyright (c) 2016, The Regents of the University of California.
# For licensing details see "license.txt" and "legal.txt".

import re
import core


class SBMLResourceIdentifiers(object):
    """
    Parses and stores a set of biology reference resource identifiers, as a dictionary.
    """

    identifiersorg_pattern = re.compile(r'^https?://identifiers.org/([a-zA-Z0-9.\-]+)/([^\s]+)$')

    def __init__(self, urls=None):
        self.dict = {}
        if urls is not None:
            self.add(urls)


    def add(self, urls):
        urlslist = urls
        if isinstance(urls, str):     # Single name
            urlslist = [urls]
        for u in urlslist:
            self.process(u)


    def process(self, url):
        result = SBMLResourceIdentifiers.match(url)
        if result is None:
            return
        if result[0] in self.dict:
            if isinstance(self.dict[result[0]], str):
                self.dict[result[0]] = [self.dict[result[0]]]
            self.dict[result[0]].append(result[1])
        else:
            self.dict[result[0]] = result[1]


    def __str__(self):
        return ' '.join(self.asArray())


    def asArray(self):
        substrs = []
        for k in self.dict:
            v = self.dict[k]
            if isinstance(v, str):
                v = [v]
            for onev in v:
                substrs.append('%s,%s' % (k, onev))
        return substrs


    @staticmethod
    def match(url):
        url_match = SBMLResourceIdentifiers.identifiersorg_pattern.match(url)
        if not url_match:
            return None
        name = url_match.group(1)
        key = url_match.group(2)
        return [name, key]



if __name__ == "__main__":

    a = SBMLResourceIdentifiers()
    a.add([ "http://identifiers.org/ncbigi/GI:16130343",
            "http://identifiers.org/ncbigene/946880",
            "http://identifiers.org/asap/ABE-0007971"])
    a.add("http://identifiers.org/ecogene/EG10165")
    a.add([ "http://identifiers.org/chebi/CHEBI:11304",
            "http://identifiers.org/chebi/CHEBI:15637",
            "http://identifiers.org/chebi/CHEBI:19108"])

    print a.dict

