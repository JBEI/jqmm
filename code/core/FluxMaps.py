# Part of the JBEI Quantitative Metabolic Modeling Library (JQMM)
# Copyright (c) 2016, The Regents of the University of California.
# For licensing details see "license.txt" and "legal.txt".

"""
The FluxMaps module provides flux visualization using svg files. 
It should be eventually phased out in favour of arrowland in:

https://public-arrowland.jbei.org/

"""

import re
import math

import core
from svgpathparser import parse_path
from lxml import etree
from StringIO import StringIO


transcriptomicsHaloColor = '#FFDC60'
proteomicsHaloColor = '#BEF7FF'

class FluxMap:
	"Class for drawing a flux map"
	def __init__(self,inputFilename='default',inputString='default',zoomLevel=3):

		if inputString != 'default':
			self.baseSVGString	= inputString
		elif inputFilename != 'default':
			# Read the file to a string
			with open(inputFilename, 'r') as file:
				self.baseSVGString = file.read()
			file.closed

		self.arrowScaleFactor = 0.2

		self.errorCount = 0
		self.zoomLevel = zoomLevel
		self.logStrings = []
		self._parseSVGTemplate()


	def _parseSVGTemplate(self):
		# Code samples that could be used if adapting to use libxml2 directly, instead of lxml:
		# Also see: http://docs.python.org/3/library/xml.dom.minidom.html#module-xml.dom.minidom
		# doc = libxml2.readDoc(self.basefile, None, None, 0)
		# root = doc.getRootElement()
		# result = doc.xpathEval('//*')
		#	for node in result:
		#	 print node.name
		#	 print node.prop("id")
		# child = root.children	 # the children property returns the FIRST child of a node
		# while child is not None:
		#	if child.type == "element":
		#	  # do something with the child node
		#	  print child.name
		#	child = child.next
		# newNode = libxml2.newNode('bar')
		# root.addChild(newNode)
		# newNode.setProp('attribute', 'the value')

		# Turns out removing blank text also removes the carriage returns that separate multiple lines
		# in text sub-elements - at least, according to Inkscape.  In other words, Inkscape treats them
		# like SPANs, while browsers treat them like DIVs (which is more correct, IMHO).
		#parser = etree.XMLParser(remove_blank_text=True, remove_comments=True)
		#etree.register_namespace('', "http://www.w3.org/2000/svg")
		parser = etree.XMLParser(remove_comments=True, ns_clean=True)
		strIODisguise = StringIO(self.baseSVGString)
		self.svgTree = etree.parse(strIODisguise, parser)
		self.svgRoot = self.svgTree.getroot()

		self.mapElements = self.svgRoot.findall(".//*[@id='map_element_group']//*")
		if not self.mapElements[0:1]:
			self.errorCount += 1
			e = 'SVG input file has missing or empty "map_element_group" group.'
			self.logStrings.append(e)
			raise TypeError(e)

		# build path and text dicts used to reference elements in the svg
		pathDict = {}
		cofactorDict = {}
		textDict = {}  
		for entity in self.mapElements:
			if not entity.attrib.has_key('id'):
				continue
			if not hasattr(entity, 'tag'):
				continue

			id = entity.attrib['id']
			tag = re.sub('\{.*\}','',entity.tag)	# Ignore the XML DTD categorization

			# ungrouped paths, and cofactors
			if tag == "path":
				if id.startswith('reaction-'):
					if id.endswith('-path'):
						rname = id[9:-5]
						rname = re.sub('-part.+$','', rname)	# Drop any sub-part specification
						self.logStrings.append('Found reaction ' + rname + ' path ' + id)
						if rname not in pathDict:
							pathDict[rname] = {}
						if id in pathDict[rname]:
							pathDict[rname][id].append(entity)
						else:
							pathDict[rname][id] = [entity]
					elif id.find('-cofactor-') != -1:	# Example ID: reaction-HEX1-cofactor-ATP-forward-consume
						rbits = id.split("-")
						rname = rbits[1]
						self.logStrings.append('Found cofactor for ' + rname)
						if rname in cofactorDict:
							cofactorDict[rname].append(entity)
						else:
							cofactorDict[rname]=[entity]

			# grouped paths
			if tag == "g":
				if id.startswith('reaction-') and id.endswith('-path'):
					rname = id[9:-5]
					rname = re.sub('-part.+$','', rname)	# Drop any sub-part specification
					# Find all the subitems in the group that are paths
					subPaths = []
					for s in entity.iterchildren():
						if not hasattr(s, 'tag'):
							continue
						if not re.search(r"\}path$", s.tag):
							continue
						subPaths.append(s)
						childID = 'none'
						if s.attrib.has_key('id'):
							childID = s.attrib['id']
						if rname not in pathDict:
							pathDict[rname] = {}
						if childID in pathDict[rname]:
							pathDict[rname][childID].append(entity)
						else:
							pathDict[rname][childID] = [entity]
					if len(subPaths):
						self.logStrings.append('Found reaction' + rname + ' group with ' + str(int(len(subPaths))) + ' paths')

			# get text entry names and numbers
			if tag == "text":
				if id.startswith('reaction-'):
					if id.endswith('-label'):
						name = id[9:-6]
						# Reactions can be represented in multiple places, so may have multiple labels
						name = re.sub('-part.+$','', name) # The particular part doesn't interest us

						if name in textDict:
							textDict[name].append(entity)
						else:
							textDict[name]=[entity]
		
		self.pathDict	= pathDict
		self.cofactorDict = cofactorDict
		self.textDict	= textDict
		
		self.processStyles()


	def processStyles(self):
		"make sure marker-start and marker-end are present"
		"make sure that stroke-width is present"
		"make sure that stroke-dasharray is present"
		if self.svgRoot is None:
			return
		# Locate the 'defs' section at he top of the svg document
		defSection = self.svgRoot.find("{http://www.w3.org/2000/svg}defs")
		# This is odd, but, you cannot FIND an element directly in the global namespace.
		# For example, here, we must specify the explicit svg namespace,
		# even though we're looking in a document with the global namespace declared as the svg namespace.
		# But in the same document, to CREATE an svg element, you must specify the global namespace {}
		# and not the explicit one {http://www.w3.org/2000/svg}, or you'll always get a tag
		# with an explicit namespace.
		# So when finding things, {} is broken, but when creating things, {} is required.
		if defSection is None:
			defSection = etree.Element("{}defs", id="map_definitions")
			self.svgRoot.append(defSection)
		self.svgDefSection = defSection

		fluxInMarker = None
		fluxOutMarker = None

		fluxInTranMarker = None	# Transcriptomics indicator marker designed to hand off the left-half of the flux-in marker
		fluxInProtMarker = None	# Proteomics indicator marker designed to hand off the right-half of the flux-in marker

		fluxOutTranMarker = None	# Transcriptomics indicator marker designed to hand off the right-half of the flux-out marker
		fluxOutProtMarker = None	# Proteomics indicator marker designed to hand off the left-half of the flux-out marker

		for entity in defSection.iterchildren():
			if not entity.attrib.has_key('id'):
				continue
			if not hasattr(entity, 'tag'):
				continue
			id = entity.attrib['id']
			tag = re.sub('\{.*\}','',entity.tag)	# Ignore the XML DTD categorization

			if tag == "marker":
				if id == "fluxInMarker":
					fluxInMarker = entity
					continue
				if id == "fluxOutMarker":
					fluxOutMarker = entity
					continue
				if id == "fluxInTranMarker":
					fluxInTranMarker = entity
					continue
				if id == "fluxInProtMarker":
					fluxInProtMarker = entity
					continue
				if id == "fluxOutTranMarker":
					fluxOutTranMarker = entity
					continue
				if id == "fluxOutProtMarker":
					fluxOutProtMarker = entity
					continue

		if fluxInMarker is None:
			self.logStrings.append('Creating fluxInMarker, missing from defs.')
			fluxInMarker = etree.Element("{}marker", id="fluxInMarker", orient="auto", refY="0", refX="0", style="overflow:visible")
			p = etree.Element("{}path", d="m -5.77,0 8.65,5 0,-10 -8.65,5 z m 1.4,0 6.3,3.65 0,-7.3 -6.3,3.65 z", transform="scale(0.8,0.8)", style="fill:#000000;fill-rule:evenodd;stroke:#000000;stroke-width:0;stroke-linejoin:miter;stroke-dasharray:none;stroke-dashoffset:0;")
			fluxInMarker.append(p)
			p = etree.Element("{}path", d="m -4.37,0 6.3,3.65 0,-7.3 -6.3,3.65 z", transform="scale(0.8,0.8)", style="fill:#ffffff;fill-rule:evenodd;stroke:#000000;stroke-width:0;stroke-linejoin:miter;stroke-dasharray:none;stroke-dashoffset:0;")
			fluxInMarker.append(p)
			defSection.append(fluxInMarker)

		if fluxInTranMarker is None:
			self.logStrings.append('Creating fluxInTranMarker, missing from defs.')
			fluxInTranMarker = etree.Element("{}marker", id="fluxInTranMarker", orient="auto", refY="0", refX="0", style="overflow:visible")
			p = etree.Element("{}path", d="m -5.77,0 9.66,8 0,-8 z", transform="scale(0.85,0.85)", style="fill:"+transcriptomicsHaloColor+";fill-rule:nonzero;stroke:#000000;stroke-width:0;stroke-linejoin:miter;stroke-dasharray:none;stroke-dashoffset:0;")
			fluxInTranMarker.append(p)
			defSection.append(fluxInTranMarker)

		if fluxInProtMarker is None:
			self.logStrings.append('Creating fluxInProtMarker, missing from defs.')
			fluxInProtMarker = etree.Element("{}marker", id="fluxInProtMarker", orient="auto", refY="0", refX="0", style="overflow:visible")
			p = etree.Element("{}path", d="m -5.77,0 9.66,-8 0,8 z", transform="scale(0.85,0.85)", style="fill:"+proteomicsHaloColor+";fill-rule:nonzero;stroke:#000000;stroke-width:0;stroke-linejoin:miter;stroke-dasharray:none;stroke-dashoffset:0;")
			fluxInProtMarker.append(p)
			defSection.append(fluxInProtMarker)

		if fluxOutMarker is None:
			self.logStrings.append('Creating fluxOutMarker, missing from defs.')
			fluxOutMarker = etree.Element("{}marker", id="fluxOutMarker", orient="auto", refY="0", refX="0", style="overflow:visible")
			p = etree.Element("{}path", d="m 5.77,0 -8.65,5 0,-10 8.65,5 z", transform="scale(0.8,0.8)", style="fill-rule:nonzero;stroke:#000000;stroke-width:0;stroke-linejoin:miter;stroke-dasharray:none;stroke-dashoffset:0;")
			fluxOutMarker.append(p)
			defSection.append(fluxOutMarker)

		if fluxOutTranMarker is None:
			self.logStrings.append('Creating fluxOutTranMarker, missing from defs.')
			fluxOutTranMarker = etree.Element("{}marker", id="fluxOutTranMarker", orient="auto", refY="0", refX="0", style="overflow:visible")
			p = etree.Element("{}path", d="m 5.77,0 -9.66,8 0,-8 z", transform="scale(0.85,0.85)", style="fill:"+transcriptomicsHaloColor+";fill-rule:nonzero;stroke:#000000;stroke-width:0;stroke-linejoin:miter;stroke-dasharray:none;stroke-dashoffset:0;")
			fluxOutTranMarker.append(p)
			defSection.append(fluxOutTranMarker)

		if fluxOutProtMarker is None:
			self.logStrings.append('Creating fluxOutProtMarker, missing from defs.')
			fluxOutProtMarker = etree.Element("{}marker", id="fluxOutProtMarker", orient="auto", refY="0", refX="0", style="overflow:visible")
			p = etree.Element("{}path", d="m 5.77,0 -9.66,-8 0,8 z", transform="scale(0.85,0.85)", style="fill:"+proteomicsHaloColor+";fill-rule:nonzero;stroke:#000000;stroke-width:0;stroke-linejoin:miter;stroke-dasharray:none;stroke-dashoffset:0;")
			fluxOutProtMarker.append(p)
			defSection.append(fluxOutProtMarker)


	def fluxToArrowScaler(self, v):
		if not v:
			return 0
		v = abs(v)
		# Above x=9 and we start going logarithmic past about 14
		if v > 9.1:
			v = 13.4 + math.log(v-6.05)
		# Between about x=2/3 and x=9, we go logarithmic on a range of approximately 8 to 14
		elif v >= 0.66:
			v = (8 * math.log((v * 0.33) + 2)) + 1.63
		# Below x=2/3, we approach 8 on a sharp curve that inflects at around x=1/3
		else:
			j = 1 - math.sin((math.pi/2.1)*v)
			v = 8 - (8 * (j*j*j))
		return v


	def changeDrawingParameters(self, arrowScaleFactor=None):
		"change map drawing parameters"
		if arrowScaleFactor is not None:
			self.arrowScaleFactor = arrowScaleFactor


	def resetSVG(self,zoomLevel=3):
		"resets svg to basesvg"
		self.zoomLevel = zoomLevel
		self.pathDict	= {}	# Blanking these beforehand may prevent a memory leak
		self.cofactorDict = {}
		self.textDict	= {}
		self._parseSVGTemplate()


	def finalizeSVG(self):
		"performs various one-time actions to finalize the map for export"
		z = self.zoomLevel
		numberRemoved = 0
		numberAltered = 0

		# Get rid of the reference markers we used
		for p in self.pathDict:
			for q in self.pathDict[p]:
				for r in self.pathDict[p][q]:
					if r.attrib.has_key('hasStartMarker'):
						r.attrib.pop('hasStartMarker')
						numberAltered += 1
					if r.attrib.has_key('hasEndMarker'):
						r.attrib.pop('hasEndMarker')
						numberAltered += 1
		if z < 2:
			# Below zoom level 2, we hide all the cofactor marks
			for name in self.cofactorDict.keys():
				for entity in self.cofactorDict[name]:
					entity.getparent().remove(entity)
					numberRemoved += 1
		if z < 3:
			# Below zoom level 3, we hide all the text labels
			for name in self.textDict.keys():
				for entity in self.textDict[name]:
					entity.getparent().remove(entity)
					numberRemoved += 1
		self.logStrings.append('Finalization removed ' + str(int(numberRemoved)) + ' entities and made ' + str(int(numberAltered)) + ' alterations')


	def writeSVG(self,outfile):
		"writes svg to outfile"
		self.finalizeSVG()
		self.svgTree.write(outfile)


	def getSVG(self):
		"returns svg as a string"
		self.finalizeSVG()
		return etree.tostring(self.svgTree)


	# Modified version of _decorateReactionPath to handle single bidirectional flux arrows
	# value: for unidirectional, handle the sign; for bidirectional, treat as magnitude
	def _decorateReactionPath(self, name, fluxValue, transcriptomicsPresent, proteomicsPresent):
		"private method, changes arrow widths"

		partsDict = None
		try:
			partsDict = self.pathDict[name]
		except KeyError:
			self.logStrings.append("Unable to find paths for reaction " + name + " to change arrow width")
			self.errorCount += 1
			return

		# Only use 'best' flux value for arrow widths
		if fluxValue.__class__.__name__ == 'rangedNumber':
			fluxValue = fluxValue.best

		parts = []
		for p in partsDict:
			for q in partsDict[p]:
				parts.append(q)

		endToPrependFrom = parts[0]

		for q in parts:
			# We need a style attribute and a path definition attribute to work with.  If none is present,
			# we're probably looking at a text node or some other improperly formed element.
			if ((not q.attrib.has_key('style')) or (not q.attrib.has_key('d'))):
				continue
			st = q.attrib['style']
			pathD = q.attrib['d']
			origID = q.attrib['id']

			# The semicolon is optional for the last piece of a style declaration
			# but to make search-and-replace easier we will ensure it exists
			if not st.endswith(';'):
				st = st + ';'

			# Search for these and remove them, and then set another attribute
			# that will persist for future reference
			markstart	= re.search('marker-start:url\([^\)]+\)', st)
			markend		= re.search('marker-end:url\([^\)]+\)', st)
			markstartnone	= 'marker-start:none' in st
			markendnone		= 'marker-end:none' in st
			if markstart and (not markstartnone):
				q.set('hasStartMarker', '1')
			if markend and (not markendnone):
				q.set('hasEndMarker', '1')

			# Rip out just about everything from the style
			st = re.sub('fill-opacity:[0-9\.]+','',st)
			st = re.sub('stroke-opacity:[\d.]+;','',st)
			st = re.sub('stroke-width:[^;]+;','',st)
			st = re.sub('stroke:[^;]+;','',st)
			st = re.sub('fill:[^;]+;','',st)
			st = re.sub('stroke-dasharray:[^;]+;','',st)
			st = re.sub('stroke-dashoffset:[^;]+;','',st)
			st = re.sub('marker-start:url\([^\)]+\);','',st)
			st = re.sub('marker-end:url\([^\)]+\);','',st)

			hasStartMarker = False
			if q.attrib.has_key('hasStartMarker'):
				hasStartMarker = True
			hasEndMarker = False
			if q.attrib.has_key('hasEndMarker'):
				hasEndMarker = True

			# If null fluxValue, mark the line as absent from the map data
			if not isinstance(fluxValue, float):
				st = st + 'stroke:#444444;fill:none;stroke-width:0.5;stroke-dasharray:0.625,1.25;stroke-opacity:0.5;'
				q.set('style', st)
				continue

			# Save this style for later, so we can use it for start and end marker paths.
			capSt = st

			# Parse the path description
			parsedPath = parse_path(pathD)
			# Find the length of the original path, by manually parsing the path description attribute,
			# and building rough estimates at each complex curve.
			# You'd think there'd be an easier way to do this, right?
			pLength = parsedPath.length()

			if fluxValue == 0:
				flatSt = st + 'stroke:#000000;fill:none;stroke-width:0.5;stroke-opacity:0.5;'
				if hasStartMarker:
					flatSt = flatSt + 'marker-start:url(#fluxInMarker);'
				if hasEndMarker:
					flatSt = flatSt + 'marker-end:url(#fluxOutMarker);'
				q.set('style', flatSt)

				wStart = 1
				wEnd = 1
				newPOffset = 0.25
				maxMagnitude = 1

				# Report what we're doing
				self.logStrings.append("Building paths for " + name + " fluxValue 0")
			else:
				# First we convert the target flux fluxValue into a stroke width
				scaledFlux	= self.fluxToArrowScaler(abs(fluxValue))

				# Then we use that to find the start width and the end width
				wStart = min(2.0,scaledFlux/2)
				wEnd = max(scaledFlux,0.5)

				# Report what we're doing
				self.logStrings.append("Building paths for " + name + " fluxValue " + str(fluxValue) + " wStart " + str(wStart) + " wEnd " + str(wEnd))

				# Rebuild the path as the definition for a solid shape
				newPOffset = wEnd/2
				newPStart = wStart/wEnd
				newPEnd = 1
				if fluxValue < 0:
					newPStart = 1
					newPEnd = wStart/wEnd
				offsetPath = parsedPath.simpleOffestPathToShape(newPOffset, 0-newPOffset, newPStart, newPEnd)
				# Embed it back in the original SVG path element
				q.set('d', offsetPath.getPathFragment())
				q.set('style', st + 'stroke:none;fill-opacity:1;fill:#000000;')

				maxMagnitude = max(newPStart, newPEnd)

			prependedPathsList = []
			appendedPathsList = []

			if transcriptomicsPresent is True:
				newSt = st + "stroke:none;fill-opacity:1;fill:"+transcriptomicsHaloColor+";"
				newOffsetPath = parsedPath.simpleOffestPathToShape(0, 0-((newPOffset*1.6)+1.1), maxMagnitude, maxMagnitude)
				pNew = etree.Element("{}path", d=newOffsetPath.getPathFragment(), id=origID + "-tomicsB", style=newSt)
				prependedPathsList.append(pNew)

			if proteomicsPresent is True:
				newSt = st + "stroke:none;fill-opacity:1;fill:"+proteomicsHaloColor+";"
				newOffsetPath = parsedPath.simpleOffestPathToShape((newPOffset*1.6)+1.1, 0, maxMagnitude, maxMagnitude)
				pNew = etree.Element("{}path", d=newOffsetPath.getPathFragment(), id=origID + "-pomicsB", style=newSt)
				prependedPathsList.append(pNew)

			# For the paths that will host the start and end markers, we want a dasharray
			# that will effectively render the path invisible, so we can stick the paths
			# (and therefore the markers) on top without obscuring the other paths.
			# This allows us to control the size of the markers independent of everything else.
			# TODO: Replace this with one or two basic near-zero length lines that reproduce the needed angles,
			# to simplify the geometry.
			capSt = capSt + 'stroke-dasharray:1.0,{1:5.3f};'.format(1,pLength+1)
			capSt = capSt + 'stroke-dashoffset:1.5;stroke-opacity:1;stroke:#000000;fill:none;'

			if hasStartMarker:
				if (fluxValue < 0):
					sw = (wEnd * self.arrowScaleFactor) + 0.4
				else:
					sw = (wStart * self.arrowScaleFactor) + 0.4
				startCapSt = capSt + 'stroke-width:{0:5.3f};'.format(sw)
				newID = origID + "-startcap"
				pNew = etree.Element("{}path", d=pathD, id=newID, style=startCapSt + 'marker-start:url(#fluxInMarker);')
				appendedPathsList.append(pNew)
				if (fluxValue < 0):
					if transcriptomicsPresent is True:
						newID = origID + "-startcaptran"
						pNew = etree.Element("{}path", d=pathD, id=newID, style=startCapSt + 'marker-start:url(#fluxInTranMarker);')
						prependedPathsList.append(pNew)
					if proteomicsPresent is True:
						newID = origID + "-startcapprot"
						pNew = etree.Element("{}path", d=pathD, id=newID, style=startCapSt + 'marker-start:url(#fluxInProtMarker);')
						prependedPathsList.append(pNew)

			if hasEndMarker:
				if (fluxValue < 0):
					sw = (wStart * self.arrowScaleFactor) + 0.4
				else:
					sw = (wEnd * self.arrowScaleFactor) + 0.4
				endCapSt = capSt + 'stroke-width:{0:5.3f};'.format(sw)
				newID = origID + "-endcap"
				pNew = etree.Element("{}path", d=pathD, id=newID, style=endCapSt + 'marker-end:url(#fluxOutMarker);')
				appendedPathsList.append(pNew)
				if (fluxValue > 0):
					if transcriptomicsPresent is True:
						newID = origID + "-endcaptran"
						pNew = etree.Element("{}path", d=pathD, id=newID, style=endCapSt + 'marker-end:url(#fluxOutTranMarker);')
						prependedPathsList.append(pNew)
					if proteomicsPresent is True:
						newID = origID + "-endcapprot"
						pNew = etree.Element("{}path", d=pathD, id=newID, style=endCapSt + 'marker-end:url(#fluxOutProtMarker);')
						prependedPathsList.append(pNew)

			for i in prependedPathsList:
				endToPrependFrom.addprevious(i)
				endToPrependFrom = i

			endToAddFrom = q
			for i in appendedPathsList:
				endToAddFrom.addnext(i)
				endToAddFrom = i


	def _changeFluxText(self, name, fluxValue):
		"private method, changes text labels"

		labels = None
		try:
			labels = self.textDict[name]
		except KeyError:
			self.logStrings.append("Unable to find text label " + name + " to change content")
			self.errorCount += 1
			return

		for label in labels:
			if label.attrib.has_key('style'):
				st = label.attrib['style']
				if not st.endswith(';'):
					st = st + ';'
			else:
				st = ';'

			parts = []
			for s in label.iterchildren():
				if not hasattr(s, 'tag'):
					continue
				if re.search(r"\}tspan$", s.tag):
					parts.append(s)
			if not parts[0:1]:
				# If there isn't even a first item, give up
				continue

			# Formatting values
			stdInsert = ''
			fill = 1
			if isinstance(fluxValue, core.rangedNumber):
				if fluxValue.lo == 0 and fluxValue.hi == 0 and fluxValue.best == 0:
					#fill = 0.5
					txtInsert = '(0.0)'
				else:
					txtInsert = '({0:5.3f})'.format(fluxValue.best)
					if abs(fluxValue.hi-fluxValue.lo) != 0:
						fLow = '{0:5.2f}'.format(fluxValue.lo)  # Float formatting introduces unnecessary spaces that
						fHigh = '{0:5.2f}'.format(fluxValue.hi) # interfere with spacing when the value is positive
						stdInsert = '(' + fLow.strip() + '/' + fHigh.strip() + ')'
			elif isinstance(fluxValue, float):
				if fluxValue == 0:
					txtInsert = '(0.0)'
					#fill = 0.5
				else:
					txtInsert = '({0:4.3G})'.format(fluxValue)
			elif fluxValue == '?':
				txtInsert = ''
				fill = 0.2
			else:
				txtInsert = fluxValue 

			st = re.sub('fill-opacity:[0-9\.]+','',st)
			st = st + 'fill-opacity:{0:5.2f};'.format(fill)
			label.set('style', st)

			# Updating svg with new text
			if parts[1:2]:	  # No need to catch an "IndexError" with this.
				parts[1].text = txtInsert
			if parts[2:3]:
				parts[2].text = stdInsert


	def changeTitle(self,title):
		"Changes title for file"
		#self.svg[self.titleLoc] = title


	def decorateMap(self, rNetwork, fluxDict):
		"changes all arrows and labels based on flux values, protein measurements, and gene transcriptions"
		# First change all present arrows to non existent (to flag non-included reactions)
		for name in self.textDict:
			self._decorateReactionPath(name, '?', None, None)
			self._changeFluxText(name, '?')

		reactDict = rNetwork.reactionList.getReactionDictionary()

		# Only do further decoration for reactions that have a flux value.
		for fluxName in fluxDict.keys():					
			value = fluxDict[fluxName]
			transcriptomicsPresent = False
			proteomicsPresent = False
			if fluxName in reactDict:
				transcriptomicsPresent = reactDict[fluxName].getSignificantTranscriptomicsValue()
				proteomicsPresent = reactDict[fluxName].getSignificantProteomicsValue()
			self._decorateReactionPath(fluxName, value, transcriptomicsPresent, proteomicsPresent)
			self._changeFluxText(fluxName, value)


	def findorphans(self, fluxDict):
		"debug info on orphanfluxes (fluxDict entries without map paths) and orphanpaths (map paths without fluxDict entries)"
		pathset = set(self.pathDict.keys())
		fluxset = set(fluxDict.keys())
		self.orphanfluxes	= fluxset - pathset
		self.orphanpaths	= pathset - fluxset




# Class for drawing flux maps
# TODO: Delete this when all maps are actualized
class FluxMapOLD:
    "class representing a flux map"

    def __init__(self,basefile,strokewidth_multiplier=3,max_strokewidth=5,min_strokewidth=0.2):
        import svgfig
        svg         = svgfig.load(basefile)
        # build path and text dicts used to reference fluxes in the svg
        rawPathList = []
        rawTextDict = {}  
        for index,item in svg:
            
            # Title info 
            try:
                attr = item.attr
                try:
                    id   = attr['id']
                    if 'title' in id:
                        titleLoc = (index,0,0)
                except KeyError:
                    pass
            except AttributeError:
                pass
            
            # The rest
            try:
                # ungrouped paths
                if item.t == "path" and len(index) == 1:
                    rawPathList.append((item.attr['id'],[index]))
                
                # grouped paths
                # if item type is group and the ID isn't a generic "gNNNN"
                if item.t == "g" and not item.attr['id'][1:].isnumeric():
                    # if every subitem in the group is a path
                    pathtest = []
                    for s in item.sub:
                        pathtest.append(s.t == "path")                  
                    if all(pathtest):
                        rawPathList.append((item.attr['id'],[(index[0],u) for u in range(len(item.sub))]))

                # get text entry names and numbers
                if item.t == "text":
                    # set up a dict lookup by name, store top-level index as the dict contents
                    # index,0,0 : the name
                    # index,1,0 : the rate number
                    name = svg[index,0,0].strip()
                    rate = (index[0],1,0)
                    stds = (index[0],2,0)
                    
                    if name in rawTextDict:
                        rawTextDict[name].append((rate,stds))
                    else:
                        rawTextDict[name]=[(rate,stds)]
                    # svgfig silently handles IndexError exception if svg[index,1,0] does not exist
            except:
                pass
        # clean misc paths and cofactors out of the path list
        cleanPathList = []
        badTerms = ["path","atpc","adpc","coac","co2c","nadc","nadhc","nadpc","nadphc"]
        for entry in rawPathList:
            flagBad = False
            for bad in badTerms:
                if bad in entry[0]:
                    flagBad = True
                    break
            if not flagBad:
                cleanPathList.append(entry)
        pathdict = dict(cleanPathList)
        # clean metabolite labels and pure numeric labels out of the text list
        cleanTextDict = {}
        badTerms = ["[c]","[e]"]
        for entry in rawTextDict:
            flagBad = False
            for bad in badTerms:
                if bad in entry:
                    flagBad = True
                    break
            #if not flagBad and entry[0].isnumeric() == False:
            if not flagBad and entry.isnumeric() == False:
                cleanTextDict[entry] = rawTextDict[entry]
        textdict = cleanTextDict
        
        #textdict = dict(cleanTextDict)
        # assign attributes
        self.basefile   = basefile
        self.basesvg    = svg
        self.svg        = svg
        self.pathdict   = pathdict
        self.textdict   = textdict
        self.strokewidth_multiplier     = strokewidth_multiplier
        self.max_strokewidth            = max_strokewidth
        self.min_strokewidth            = min_strokewidth
        self.titleLoc   = titleLoc



    def changestrokewidths(self,strokewidth_multiplier=3,max_strokewidth=5,min_strokewidth=0.2):
        "change strokewidth from initialization settings"
        self.strokewidth_multiplier     = strokewidth_multiplier
        self.max_strokewidth            = max_strokewidth
        self.min_strokewidth            = min_strokewidth       


    def resetsvg(self):
        "resets svg to basesvg"
        self.svg = self.basesvg


    def writesvg(self,outfile):
        "writes svg to outfile"
        self.svg.save(outfile)


    def _changearrowwidth(self,name,value):
        "private method, changes arrow widths"
        # Only use best value for arrow widths
        if value.__class__.__name__ == 'rangedNumber':
            value = value.best
        # If null value, make it zero
        if not isinstance(value, float):
            value = 0

        if value == 0:
            for q in self.pathdict[name]:
                self.svg[q].attr['style'] = re.sub('stroke-width:[\d.]+','stroke-width:1.0',self.svg[q].attr['style'])
                self.svg[q].attr['style'] = re.sub('stroke-dasharray:none','stroke-dasharray:0.625, 1.25',self.svg[q].attr['style'])
                self.svg[q].attr['style'] = re.sub('marker-end:url\(.*\)','marker-end:url(#TriangleOutS)',self.svg[q].attr['style'])
        else:
            for q in self.pathdict[name]:
                # Changing line width
                finalvalue  = min(self.strokewidth_multiplier*abs(value),self.max_strokewidth);
                finalvalue  = max(finalvalue,self.min_strokewidth);
                insert = 'stroke-width:{0:5.3f}'.format(finalvalue)
                self.svg[q].attr['style'] = re.sub('stroke-width:[\d.]+',insert,self.svg[q].attr['style'])
                self.svg[q].attr['style'] = re.sub('stroke-dasharray:[^;]','stroke-dasharray:none',self.svg[q].attr['style'])
                #Arrow point correction
                if 10*abs(value) >= 3:
                    self.svg[q].attr['style'] = re.sub('marker-end:url\(.*\)','marker-end:url(#TriangleOutS)',self.svg[q].attr['style'])
                elif 10*abs(value) >= 0.8 and 10*abs(value) < 3:
                    self.svg[q].attr['style'] = re.sub('marker-end:url\(.*\)','marker-end:url(#TriangleOutM)',self.svg[q].attr['style'])
                else:
                    self.svg[q].attr['style'] = re.sub('marker-end:url\(.*\)','marker-end:url(#TriangleOutL)',self.svg[q].attr['style'])


    def _changetextlabel(self,name,value):
        "private method, changes text labels"
        for group in self.textdict[name]:
            # Getting coordinates in svg tree for text and stds
            if group.__class__.__name__ == 'tuple' and group[0].__class__.__name__ == 'tuple':
                txtCoord,stdCoord = group
            else:
                txtCoord = group
                stdCoord = () 
            
            # Changing text
            if value.__class__.__name__ == 'rangedNumber':
                txtInsert = '{0:5.3f}'.format(value.best)
                stdInsert = '{0:5.2f}-{1:5.2f}'.format(value.lo,value.hi)
            else:
                if isinstance(value, float):
                    txtInsert = '{0:5.3f}'.format(value) 
                    stdInsert = 'NA-NA'
                else:
                    txtInsert = value 
                    stdInsert = value
            
            # Updating svg with new text
            self.svg[txtCoord] = '('+txtInsert+')'
            self.svg[stdCoord] = '('+stdInsert+')'


    def changetitle(self,title):
        "Changes title for file"
        self.svg[self.titleLoc] = title


    def changeflux(self,name,value):
        "changes single flux arrow and label" 
        # Deal with two possible input types for value
        # This should be made more compact
        if value.__class__.__name__ == 'rangedNumber':
            if value.best < 0:
                netname = name+'-r'
                backname = name
                value.best = -value.best
                value.lo   = -value.lo
                value.hi   = -value.hi
            else:
                netname = name
                backname = name+'-r'
        else:
            netname = name
            backname = name+'-r'
            if isinstance(value, float):
                if value < 0:
                    netname = name+'-r'
                    backname = name
                    value = -value          
            
        # Change back arrow and forward arrow    
        try:
            self._changearrowwidth(backname,0)
        except KeyError:
            pass
        try:
            self._changearrowwidth(netname,value)
            self._changetextlabel(name,value)
        except KeyError:
            pass


    def changeallfluxes(self,fluxdict):
        "changes all arrows and labels based on "
        # First change all present arrows to non existent (to flag non-included reactions)
        for name in self.textdict:
            self.changeflux(name,'---')
        
        # Then change flux values
        for name in fluxdict.keys():                    
            value = fluxdict[name]
            self.changeflux(name,value)


    def findorphans(self,fluxdict):
        "debug info on orphanfluxes (fluxdict entries without map paths) and orphanpaths (map paths without fluxdict entries)"
        pathset = set(self.pathdict.keys())
        fluxset = set(fluxdict.keys())
        self.orphanfluxes   = fluxset - pathset
        self.orphanpaths    = pathset - fluxset



############### Tests ##################
if __name__ == "__main__":
	
	print "No tests..."
