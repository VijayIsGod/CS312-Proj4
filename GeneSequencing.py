#!/usr/bin/python3

from cmath import inf
from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		pass
	
	# Generates alignment using the back_pointers dictionary generated while executing either banded or unrestricted.
	def alignment(self, seq1,seq2):
		self.alignment1 = []
		self.alignment2 = []
		current = self.back_pointers[self.min_length1-1, self.min_length2-1]
		i = self.min_length1-1
		j = self.min_length2 -1
		while(current[0] > -1 or current[1] > -1):
			if(current[2] == 'top'):
				self.alignment2.insert(0,'-')
				self.alignment1.insert(0,seq1[i])
			if(current[2] == 'left'):
				self.alignment1.insert(0, '-')
				self.alignment2.insert(0, seq2[j])
			if(current[2] == 'diagonal'):
				self.alignment1.insert(0,seq1[i])
				self.alignment2.insert(0,seq2[j])
			i = current[0]
			j = current[1]
			current = self.back_pointers[current[0], current[1]]
		self.alignment1.insert(0,seq1[i])
		self.alignment2.insert(0,seq2[j])
		if(current[2] == 'top'):
			self.alignment2.insert(0,'-')
		if(current[1] == 'left'):
			self.alignment1.insert(0,'-')
		return

	#Finds the min of the three values and return from which direction the cell would go to current cell.
	# Used in banded or unrestricted
	def getMin(self, a, b, c):
		minimum = min(a,b,c)
		if(minimum == a):
			return 'left'
		if(minimum == b):
			return 'top'
		if(minimum == c):
			return 'diagonal'
	
	# Algorithm for implementing the unrestricted alignment.
	def unrestricted(self, seq1, seq2):
		self.score = 0
		self.matrix[(-1,-1)] = 0
		if(self.MaxCharactersToAlign > len(seq1)):
			self.min_length1 = len(seq1)
		else:
			self.min_length1 = self.MaxCharactersToAlign
		if(self.MaxCharactersToAlign > len(seq2)):
			self.min_length2 = len(seq2)
		else:
			self.min_length2 = self.MaxCharactersToAlign
		for i in range(-1,self.min_length1): #Time complexity of O(m) where m is size of seq1
			for j in range(-1, self.min_length2): #Time complexity of O(n) where n is size of seq2.
				if(i==-1 and j==-1):
					continue
				if(i == -1 and j >-1) :
					self.matrix[-1,j] = self.matrix[-1,j-1] + INDEL
					self.back_pointers[-1,j] = (-1,j-1,'left')
					continue
				if(i > -1 and j == -1):
					self.matrix[i,-1] = self.matrix[i-1, -1] + INDEL
					self.back_pointers[i,-1] = (i-1, -1, 'top')
					continue
				if(seq1[i] == seq2[j]):
					diagonal = self.matrix[i-1,j-1] + MATCH
					direction = self.getMin(self.matrix[i,j-1]+INDEL, self.matrix[i-1,j]+INDEL,diagonal)
				else:
					diagonal = self.matrix[i-1,j-1] + SUB
					direction = self.getMin(self.matrix[i,j-1]+INDEL, self.matrix[i-1,j]+INDEL,diagonal)
				if(direction == 'left'):
					self.matrix[i,j] = self.matrix[i,j-1] + INDEL
					self.back_pointers[i,j] = (i,j-1,'left')
					continue
				if(direction == 'top'):
					self.matrix[i,j] = self.matrix[i-1,j] + INDEL
					self.back_pointers[i,j] = (i-1,j, 'top')
					continue
				if(direction == 'diagonal'):
					self.matrix[i,j] = diagonal
					self.back_pointers[i,j] = (i-1,j-1,'diagonal')
					continue
		self.score = self.matrix[self.min_length1-1,self.min_length2-1]

	# Algorithm for implementing the banded alignment.
	def banded_alignment(self, seq1, seq2):
		self.score = 0
		self.matrix[(-1,-1)] = 0
		if(self.MaxCharactersToAlign > len(seq1)):
			self.min_length1 = len(seq1)
		else:
			self.min_length1 = self.MaxCharactersToAlign
		if(self.MaxCharactersToAlign > len(seq2)):
			self.min_length2 = len(seq2)
		else:
			self.min_length2 = self.MaxCharactersToAlign
		for i in range(-1,self.min_length1): #Time Complexity of O(n).
			if(-1 >= i - MAXINDELS):
				j_start = -1
			else:
				j_start = i - MAXINDELS
			if(self.min_length2 <= i + MAXINDELS):
				j_end = self.min_length2
			else:
				j_end = i + MAXINDELS + 1
			for j in range(j_start, j_end): #Time Complexity of O(K) which is 7 in this case.
				if(i==-1 and j==-1):
					continue
				if(i == -1 and j >-1) :
					self.matrix[-1,j] = self.matrix[-1,j-1] + INDEL
					self.back_pointers[-1,j] = (-1,j-1,'left')
					continue
				if(i > -1 and j == -1):
					self.matrix[i,-1] = self.matrix[i-1, -1] + INDEL
					self.back_pointers[i,-1] = (i-1, -1, 'top')
					continue
				self.matrix.setdefault((i,j-1), float(inf))
				self.matrix.setdefault((i-1,j),float(inf))
				if(seq1[i] == seq2[j]):
					diagonal = self.matrix[i-1,j-1] + MATCH
					direction = self.getMin(self.matrix[i,j-1]+INDEL, self.matrix[i-1,j]+INDEL,diagonal)
				else:
					diagonal = self.matrix[i-1,j-1] + SUB
					direction = self.getMin(self.matrix[i,j-1]+INDEL, self.matrix[i-1,j]+INDEL,diagonal)
				if(direction == 'left'):
					self.matrix[i,j] = self.matrix[i,j-1] + INDEL
					self.back_pointers[i,j] = (i,j-1,'left')
					continue
				if(direction == 'top'):
					self.matrix[i,j] = self.matrix[i-1,j] + INDEL
					self.back_pointers[i,j] = (i-1,j, 'top')
					continue
				if(direction == 'diagonal'):
					self.matrix[i,j] = diagonal
					self.back_pointers[i,j] = (i-1,j-1,'diagonal')
					continue
		self.matrix.setdefault((self.min_length1-1,self.min_length2-1), float(inf))
		self.score = self.matrix[self.min_length1-1,self.min_length2-1]

# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you 
# how many base pairs to use in computing the alignment	
	def align( self, seq1, seq2, banded, align_length):
		if(seq1 == 'polynomial' and seq2 == 'polynomial'):
			self.no_of_calls = 0
		self.no_of_calls += 1
		self.banded = banded
		self.MaxCharactersToAlign = align_length
		self.back_pointers = {}
		self.matrix = {}
		if(self.banded):
			self.banded_alignment(seq1, seq2)
		else:
			self.unrestricted(seq1, seq2)
		if(self.score == float(inf)):
			return {'align_cost':self.score, 'seqi_first100': 'No Alignment Possible', 'seqj_first100':'No Alignment Possible'}
		self.alignment(seq1,seq2)

		# gets the first 100 characters of the alignment string.
		seqi100 = ''.join(self.alignment1[:100])	
		seqj100 = ''.join(self.alignment2[:100])

		# Prints out the first 100 strings of each alignment as the gui was cutting off the values.
		print(str(self.no_of_calls) +': ' + seqi100)
		print(str(self.no_of_calls) + ': ' + seqj100)
		print('')
		return {'align_cost':self.score, 'seqi_first100':seqi100, 'seqj_first100':seqj100}


