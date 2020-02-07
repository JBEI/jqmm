# Part of the JBEI Quantitative Metabolic Modeling Library (JQMM)
# Copyright (c) 2016, The Regents of the University of California.
# For licensing details see "license.txt" and "legal.txt".

# Code heavily adapted from "svg.path", https://github.com/regebro/svg.path (MIT license)
# Copyright (c) 2013-2014 Lennart Regebro
"""
This file contains classes for the different types of SVG path segments as
well as a Path object that manipulates a sequence of path segments.
"""

from __future__ import division
# In older Pythons, if you divide two integer values, you will get an integer:
# 1 / 2.0	--> 0.5
# 1 / 2		--> 0
# In the future, Python will switch to always yielding a real result,
# and to force an integer division operation you use the special "//"
# integer division operator.
from builtins import zip
from builtins import range
from builtins import object
from math import sqrt, cos, sin, acos, degrees, radians, pi
from collections import MutableSequence



class CubicBezier(object):
	def __init__(self, start, control1, control2, end):
		self.start = start
		self.control1 = control1
		self.control2 = control2
		self.end = end
		self._length = None

	
	def __repr__(self):
		return '<CubicBezier start=%s control1=%s control2=%s end=%s>' % (
			   self.start, self.control1, self.control2, self.end)


	def __eq__(self, other):
		if not isinstance(other, CubicBezier):
			return NotImplemented
		return self.start == other.start and self.end == other.end and \
			   self.control1 == other.control1 and self.control2 == other.control2


	def __ne__(self, other):
		if not isinstance(other, CubicBezier):
			return NotImplemented
		return not self == other


	def getPathFragment(self):
		return 'C %s,%s %s,%s %s,%s' % (self.control1.real, self.control1.imag, self.control2.real, self.control2.imag, self.end.real, self.end.imag)


	def point(self, pos=0.0):
		"""Calculate the x,y position at a certain fractional position (0.0 to 1.0) of the path"""
		# Equivalent to findDotAtSegment
		return ((1-pos) ** 3 * self.start) + \
			   (3 * (1-pos) ** 2 * pos * self.control1) + \
			   (3 * (1-pos) * pos ** 2 * self.control2) + \
			   (pos ** 3 * self.end)


	def base3(self, t):
		"""Calculate the derivative at fractional position t"""
		t1 = -3 * self.start + 9 * self.control1 - 9 * self.control2 + 3 * self.end
		t2 = t * t1 + 6 * self.start - 12 * self.control1 + 6 * self.control2
		return t * t2 - 3 * self.start + 3 * self.control1


	def reverse(self):
		temp = self.start
		self.start = self.end
		self.end = temp
		temp = self.control1
		self.control1 = self.control2
		self.control2 = temp
	

	def alternate_length(self):
		"""Calculate the length of the path"""
		if self._length is not None:
			return self._length
		# It's impossible to integrate a Cubic Bezier, so
		# this is a geometric approximation instead.

		current_point = self.point(0)
		# Using 100,000 subdivisions will satisfy assertAlmostEqual on the
		# Arc segment, and we should generally go for the same here,
		# but right now we're tuning for speed rather than accuracy, so
		# we cut it way down to 1000.
		# Over 1,000,000 subdivisions makes no difference in accuracy at all.
		subdivisions = 1000
		lengthb = 0
		delta = 1/subdivisions
		
		for x in range(1, subdivisions+1):
			next_point = self.point(delta*x)
			distance = sqrt((next_point.real - current_point.real)**2 + (next_point.imag - current_point.imag)**2)
			lengthb += distance
			current_point = next_point

		self._length = lengthb
		return lengthb


	def length(self, z=1.0):
		"""Calculate the length of the path up to a certain position (fractional between 0.0 and 1.0)"""
		# It's impossible to integrate a Cubic Bezier, so
		# this is an optimized geometric approximation instead.
		z = max(0, min(1, z))
		z2 = z / 2
		n = 12
		tValues = [-0.1252, 0.1252, -0.3678, 0.3678, -0.5873, 0.5873, -0.7699, 0.7699, -0.9041, 0.9041, -0.9816, 0.9816]
		cValues = [0.2491, 0.2491, 0.2335, 0.2335, 0.2032, 0.2032, 0.1601, 0.1601, 0.1069, 0.1069, 0.0472, 0.0472]
		sum = 0
		for i in range(len(tValues)):
			ct = z2 * tValues[i] + z2
			base = self.base3(ct)
			sum += cValues[i] * sqrt(base.real**2 + base.imag**2)
		return z2 * sum


	def getTatLen(self, ll):
		if ll < 0 or self.length() < ll:
			return
		t = 1
		step = t / 2.0
		t2 = t - step
		e = .01
		l = self.length(t2)
		while abs(l - ll) > e:
			step /= 2
			if l < ll:
				t2 += step
			else:
				t2 -= step
			l = self.length(t2)
		return t2


	def getPointAtSegmentLength(self, length=None):
		if length is None:
			self.point()
		l = self.getTatLen(length)
		return self.point(l)


	def simpleSingleOffsetPath(self, offset, startMagnitude, endMagnitude):
		"""A hideously expensive function that more-or-less properly calculates an offset path
		to the curve, with the offset scaled by a fractional start and end magnitude (0.0 to 1.0).
		The newly constructed curve object is returned.
		To get an offset path on the other side of the curve, pass in a negative value for the offset."""

		# If the start point and first control point match, we'll get mathematical ugliness.
		# So we attempt to address that here by tweaking the points.  Hooray for "good enough" engineering!
		if self.start == self.control1:
			if self.control1 == self.control2:
				if self.control2 == self.end:
					# If every parameter is the same, there's no way to make a meaningful offset. Forget it.
					return CubicBezier(self.start, self.control1, self.control2, self.end)
				self.control1 = (self.start * 0.7) + (self.end * 0.3)
			else:
				self.control1 = (self.start * 0.999) + (self.control2 * 0.001)
		if self.control2 == self.end:
			self.control2 = (self.end * 0.999) + (self.control1 * 0.001)

		# Make appropriately proportional start/end versions of the intended offset
		moveStart = offset * startMagnitude
		moveEnd = offset * endMagnitude

		self._coordForTValue = []
		# Create a lookup table for resolving coordinate -> 't' value
		# This is the computationally expensive part of the whole process
		for idx in range(317):
			self._coordForTValue.append(self.point(idx / 316.0))

		tValues = []
		# Find needed t values
		tValues.append(0)
		tValues.append(self._getPointProjection(self.control1))
		tValues.append(self._getPointProjection(self.control2))
		#tValues.append(self._getPointProjection(self.end))		# Should always be 1?
		tValues.append(1)

		normals = []
		# Compute normals at t values
		ca = cos(pi / -2)
		sa = sin(pi / -2)
		for i in range(4):
			d = self.base3(tValues[i])
			nx = d.real * ca - d.imag * sa
			ny = d.real * sa + d.imag * ca
			# If the bezier is constructed with a forgiving editing tool like Inkscape,
			# a user can end up specifying a curve with control points identical to its start or end points,
			# causing slope calculations to turn out zero or infinity, in turn causing exceptions in this code.
			#  That's one of the reasons why the code elsewhere in this module is careful to calculate novel control
			# points when converting, for example, straight lines to bezier curves.
			dst = sqrt(nx*nx+ny*ny)
			normals.append(nx/dst + (ny/dst) * 1j)

		# Compute ratios
		curveLength = self.length()
		ratio1 = self.length(tValues[1]) / curveLength
		ratio2 = self.length(tValues[2]) / curveLength

		newStart =    self.start +    (moveStart * normals[0])
		newControl1 = self.control1 + (moveStart + (moveEnd - moveStart) * ratio1) * normals[1]
		newControl2 = self.control2 + (moveStart + (moveEnd - moveStart) * ratio2) * normals[2]
		newEnd =      self.end +      (moveEnd * normals[3])

		return CubicBezier(newStart, newControl1, newControl2, newEnd)


	def _getPointProjection(self, pp):
		"""find an approximate t value that acts as the control's
		projection onto the curve, towards the origin."""

		t = 0.5
		mindist = 9999999
		# find a reasonable initial "t"
		for idx in range(317):
			n = pp - self._coordForTValue[idx]
			pdist = sqrt(n.real**2 + n.imag**2)

			if pdist < mindist:
				mindist = pdist
				t = idx / 317.0
		return self._refineProjection(pp, t, mindist, 1.0/(1.01*317))


	def _refineProjection(self, pp, t, distance, precision):
		"""Refine a point projection's [t] value."""
		if precision < 0.0001:
			return t
		# refinement
		prev_t = t - precision
		next_t = t + precision

		if prev_t >= 0:
			prev = self.getPointAtSegmentLength(prev_t)
			p = pp - prev
			prev_distance = sqrt(p.real**2 + p.imag**2)

			if prev_distance < distance:
				return self._refineProjection(pp, prev_t, prev_distance, precision)

		if next_t <= 1:
			next = self.getPointAtSegmentLength(next_t)
			n = pp - next
			next_distance = sqrt(n.real**2 + n.imag**2)

			if next_distance < distance:
				return self._refineProjection(pp, next_t, next_distance, precision)

		return self._refineProjection(pp, t, distance, precision / 2.0)


	def remakeAsCubicBezier(self):
		return CubicBezier(self.start, self.control1, self.control2, self.end)

	
	
class QuadraticBezier(CubicBezier):
	# For Quadratic Bezier we simply subclass the Cubic. This is less efficient
	# and gives more complex calculations, but reuse means fewer bugs.
	# It is possible to calculate the length of a quadratic bezier so a TODO is to
	# replace the geometric approximation here.

	def __init__(self, start, control, end):
		self.start = start
		_13 = 1 / 3.0
		_23 = 2 / 3.0
		self.control1 = _13 * start + _23 * control
		self.control2 = _13 * end + _23 * control
		self.end = end


	def __repr__(self):
		return '<QuadradicBezier start=%s control=%s end=%s>' % (
			   self.start, self.control1, self.end)


	def reverse(self):
		temp = self.start
		self.start = self.end
		self.end = temp
		temp = self.control1
		self.control1 = self.control2
		self.control2 = temp


	def remakeAsCubicBezier(self):
		return CubicBezier(self.start, self.control1, self.control2, self.end)



class Line(object):

	def __init__(self, start, end):
		self.start = start
		self.end = end
	
	
	def __repr__(self):
		return '<Line start=%s end=%s>' % (self.start, self.end)


	def __eq__(self, other):
		if not isinstance(other, Line):
			return NotImplemented
		return self.start == other.start and self.end == other.end


	def __ne__(self, other):
		if not isinstance(other, Line):
			return NotImplemented
		return not self == other


	def getPathFragment(self):
		return 'L %s,%s' % (self.end.real, self.end.imag)

	
	def point(self, pos):
		distance = self.end - self.start
		return self.start + distance * pos


	def reverse(self):
		temp = self.start
		self.start = self.end
		self.end = temp
		
		
	def length(self):
		distance = (self.end - self.start)
		# It's actually possible to write this as "abs(distance)".
		# abs() returns the 'magnitude' of a complex number if it's passed one.
		# But that's a pretty mean trick to play on whoever's maintaining this code, right?
		return sqrt(distance.real**2+distance.imag**2)


	def remakeAsCubicBezier(self):
		return CubicBezier(self.start, self.point(0.33), self.point(0.66), self.end)
		


class Arc(object):

	def __init__(self, start, radius, rotation, arc, sweep, end):
		"""radius is complex, rotation is in degrees, 
		   large and sweep are 1 or 0 (True/False also work)"""

		self.start = start
		self.radius = radius
		self.rotation = rotation
		self.arc = bool(arc)
		self.sweep = bool(sweep)
		self.end = end
		self._length = None

		self._parameterize()


	def __repr__(self):
		return '<Arc start=%s radius=%s rotation=%s arc=%s sweep=%s end=%s>' % (
			   self.start, self.radius, self.rotation, self.arc, self.sweep, self.end)


	def __eq__(self, other):
		if not isinstance(other, Arc):
			return NotImplemented
		return self.start == other.start and self.end == other.end and \
			   self.radius == other.radius and self.rotation == other.rotation and\
			   self.arc == other.arc and self.sweep == other.sweep


	def __ne__(self, other):
		if not isinstance(other, Arc):
			return NotImplemented
		return not self == other


	def getPathFragment(self):
		return 'A %s,%s %s %s,%s %s,%s' % (self.radius.real, self.radius.imag, self.rotation, int(arc), int(sweep), self.end.real, self.end.imag)


	def reverse(self):
		temp = self.start
		self.start = self.end
		self.end = temp
		# Got to reverse the sweep flag
		# Effect observable in http://www.w3.org/TR/SVG11/images/paths/arcs02.svg
		self.sweep = not self.sweep
		self._parameterize()


	def _parameterize(self):
		# Conversion from endpoint to center parameterization
		# http://www.w3.org/TR/SVG/implnote.html#ArcImplementationNotes

		cosr = cos(radians(self.rotation))
		sinr = sin(radians(self.rotation))
		dx = (self.start.real - self.end.real) / 2
		dy = (self.start.imag - self.end.imag) / 2
		x1prim = cosr * dx + sinr * dy
		x1prim_sq = x1prim * x1prim
		y1prim = -sinr * dx + cosr * dy
		y1prim_sq = y1prim * y1prim

		rx = self.radius.real
		rx_sq = rx * rx
		ry = self.radius.imag		 
		ry_sq = ry * ry

		# Correct out of range radii
		radius_check = (x1prim_sq / rx_sq) + (y1prim_sq / ry_sq)
		if radius_check > 1:
			rx *= sqrt(radius_check)
			ry *= sqrt(radius_check)
			rx_sq = rx * rx
			ry_sq = ry * ry

		t1 = rx_sq * y1prim_sq
		t2 = ry_sq * x1prim_sq
		c = sqrt(abs((rx_sq * ry_sq - t1 - t2) / (t1 + t2)))
		
		if self.arc == self.sweep:
			c = -c
		cxprim = c * rx * y1prim / ry
		cyprim = -c * ry * x1prim / rx

		self.center = complex((cosr * cxprim - sinr * cyprim) + 
							  ((self.start.real + self.end.real) / 2),
							  (sinr * cxprim + cosr * cyprim) + 
							  ((self.start.imag + self.end.imag) / 2))

		ux = (x1prim - cxprim) / rx
		uy = (y1prim - cyprim) / ry
		vx = (-x1prim - cxprim) / rx
		vy = (-y1prim - cyprim) / ry
		n = sqrt(ux * ux + uy * uy)
		p = ux
		theta = degrees(acos(p / n))
		if uy < 0:
			theta = -theta
		self.theta = theta % 360

		n = sqrt((ux * ux + uy * uy) * (vx * vx + vy * vy))
		p = ux * vx + uy * vy
		if p == 0:
			delta = degrees(acos(0))
		else:
			delta = degrees(acos(p / n))
		if (ux * vy - uy * vx) < 0:
			delta = -delta
		self.delta = delta % 360
		if not self.sweep:
			self.delta -= 360


	def point(self, pos):
		angle = radians(self.theta + (self.delta * pos))
		cosr = cos(radians(self.rotation))
		sinr = sin(radians(self.rotation))
	
		x = cosr * cos(angle) * self.radius.real - sinr * sin(angle) * self.radius.imag + self.center.real
		y = sinr * cos(angle) * self.radius.real + cosr * sin(angle) * self.radius.imag + self.center.imag
		return complex(x, y)
	
	
	def length(self):
		"""The length of an elliptical arc segment requires numerical
		integration, and in that case it's simpler to just do a geometric
		approximation, as for cubic bezier curves.
		"""

		if self._length is not None:
			return self._length
		
		current_point = self.point(0)
		# Here I need 100,000 subdivisions to satisfy assertAlmostEqual. It's
		# a bit slow, but I'm not in a hurry. Over 1,000,000 subdivisions
		# makes no difference in accuracy at all.
		subdivisions = 1000	  # Note: Cut this down since we really don't need that kind of accuracy for these maps
		lengthb = 0
		delta = 1/subdivisions
		
		for x in range(1, subdivisions+1):
			next_point = self.point(delta*x)
			distance = sqrt((next_point.real - current_point.real)**2 + (next_point.imag - current_point.imag)**2)
			lengthb += distance
			current_point = next_point

		self._length = lengthb
		return lengthb
	

	# This is an ugly conversion and is generally not needed when working with metabolic maps drawn in Inkscape.
	# So for now we're skipping it.  Code that can be adapted for this purpose is in the svgUtilities.js file.
	def remakeAsCubicBezier(self):
		return NotImplemented
		#return CubicBezier(self.start, self.point(0.33), self.point(0.66), self.end)



class Path(MutableSequence):
	"""A Path is a sequence of path segments"""

	def __init__(self, *segments):
		self._segments = list(segments)
		self._length = None
		self._lengths = None
		self._isAllBeziers = False


	def __getitem__(self, index):
		return self._segments[index]


	def __setitem__(self, index, value):
		self._segments[index] = value


	def __delitem__(self, index):
		del self._segments[index]


	def insert(self, index, value):
		self._segments.insert(index, value)


	def __len__(self):
		return len(self._segments)


	def __repr__(self):
		return '<Path %s>' % ', '.join(repr(x) for x in self._segments)


	def __eq__(self, other):
		if not isinstance(other, Path):
			return NotImplemented
		if len(self) != len(other):
			return False
		for s, o in zip(self._segments, other._segments):
			if not s == o:
				return False
		return True
		

	def __ne__(self, other):
		if not isinstance(other, Path):
			return NotImplemented
		return not self == other


	def _calc_lengths(self):
		if self._length is not None:
			return
		
		lengths = [each.length() for each in self._segments]		
		self._length = sum(lengths)
		self._lengths = [each/self._length for each in lengths]


	# Render the parts as a string of SVG path instructions
	# Note that this code assumes the path is meant to be unbroken.
	# One could introduce code to compare .start and .end values and embed M commands if they differ,
	# but that might give rise to unintentionally broken shapes when the user is, for example,
	# creating an offset path as a series of smaller offset fragments (like simpleOffestPathToShape does below)
	def getPathFragment(self):

		l = len(self)
		if l < 1:
			return ''

		parts = []
		first = self._segments[0]
		parts.append('M %s,%s' % (first.start.real, first.start.imag))

		e = None
		for each in self._segments:
			e = each
			parts.append(each.getPathFragment())
		# If the path ends with a line linking back to the start point,
		# assume that the intention is to close the shape, and replace the line with a Z command
		if isinstance(e, Line):
			if first.start == e.end:
				parts.pop()
				parts.append('Z')

		return ' '.join(parts)


	def convertToCubicBeziers(self):
		if self._isAllBeziers:
			return
		l = len(self)
		i = 0
		while i < l:
			e = self._segments[i]
			self._segments[i] = e.remakeAsCubicBezier()
			i += 1
		self._isAllBeziers = True


	def reverse(self):
		for each in self._segments:
			each.reverse()
		self._segments = list(reversed(self._segments))


	def point(self, pos):
		"""Return the point along the path at the given fractional location (0.0 to 1.0)"""
		self._calc_lengths()
		# Find which segment the point we search for is located on:
		segment_start = 0
		for index, segment in enumerate(self._segments):
			segment_end = segment_start + self._lengths[index]
			if segment_end >= pos:
				# This is the segment! How far in on the segment is the point?
				segment_pos = (pos - segment_start) / (segment_end - segment_start)
				break
			segment_start = segment_end
		else:
			# This happens when pos is 1.0, and accumulated errors
			# mean that segment_end of the last segment is not quite 1.0.
			segment_pos = 1.0

		return segment.point(segment_pos)


	def length(self):
		self._calc_lengths()
		return self._length


	def simpleOffestPath(self, offset, startMagnitude, endMagnitude):
		"""This routine takes the current path and creates a path alongside it at the given offset distance,
		magnified gradually by a fractional start and end value (ranging from 0.0 to 1.0)."""

		# The necessary offset function is implemented only for bezier curves right now.
		# ... But they could pretty easily be written for lines if someone cared to make the optimization
		self.convertToCubicBeziers()
		# We should probably truncate the path to the first continuous segment,
		# just in case someone has passed in a path created from svg with "move" commands scattered in.
		# For now, we'll just assume they know what they're doing.
		totalLength = self.length()
		lengthSoFar = 0

		newPathA = []

		for p in self._segments:
			segLength = p.length()

			subStartMag = startMagnitude + (endMagnitude - startMagnitude) * (lengthSoFar / totalLength)
			subEndMag = startMagnitude + (endMagnitude - startMagnitude) * ((lengthSoFar+segLength) / totalLength)

			newPathA.append( p.simpleSingleOffsetPath(offset, subStartMag, subEndMag) )
			lengthSoFar += segLength

		newPath = Path()
		newPath.extend(newPathA)

		return newPath


	def simpleOffestPathToShape(self, offsetA, offsetB, startMagnitude, endMagnitude):
		"""This funky routine takes the current path and creates two paths that follow alongside it,
		at the given offset distances, with both offsets magnified gradually by a fractional start and end value
		(ranging from 0.0 to 1.0).  Then it takes the left-hand path, reverses it, and tacks it onto the
		right-hand path, joining the two with a line command in between, then capping the path with another line
		back to the starting point.
		This effectively creates a path that is suitable for a solid shape, representing a gradually thickening
		line that follows along the original."""

		# The necessary offset function is implemented only for bezier curves right now.
		# ... But they could pretty easily be written for lines if someone cared to make the optimization
		self.convertToCubicBeziers()

		# We should probably truncate the path to the first continuous segment,
		# just in case someone has passed in a path created from svg with "move" commands scattered in.
		# For now, we'll just assume they know what they're doing, even if the resulting "shape" is not a solid one.
		totalLength = self.length()
		lengthSoFar = 0

		newPathA = []
		newPathB = []

		for p in self._segments:
			segLength = p.length()

			subStartMag = startMagnitude + (endMagnitude - startMagnitude) * (lengthSoFar / totalLength)
			subEndMag = startMagnitude + (endMagnitude - startMagnitude) * ((lengthSoFar+segLength) / totalLength)

			newPathA.append( p.simpleSingleOffsetPath(offsetA, subStartMag, subEndMag) )
			newPathB.append( p.simpleSingleOffsetPath(offsetB, subStartMag, subEndMag) )
			lengthSoFar += segLength

		# Reverse the segments in path B individually, then reverse the order of the segments
		for p in newPathB:
			p.reverse()
		newPathB = list(reversed(newPathB))

		startPos = newPathA[-1].end
		# Create and append the connective Lines
		newPathA.append(Line(newPathA[-1].end, newPathB[0].start))
		newPathB.append(Line(newPathB[-1].end, newPathA[0].start))
		# Finally, join the paths

		newPath = Path()
		newPath.extend(newPathA)
		newPath.extend(newPathB)

		return newPath
