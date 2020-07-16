import random
import math
import numpy as np
import itertools
import re 


#given a list of known obstructions, and patterns considered to be maximal, and a max length to consider, finds what patterns are not contained in one of the maximal patterns, that are not known obstructions
#set is the alphabet to generate codes over
#only works for 0 fully specified columns in the sense that it checks all submatrices
#with more, I think what needs to happen is to apply this process to each side ? I guess in the check noncontig substring, it just needs to only work between 0s 
def generateNotContained(obstructions, maximal, set, length):
	for i in range(1, length+1):
		print(i)
		codes = generateCodes(set, i)
		for j in codes:
			isvalid = True 
			#checks for 3 in a row
			for k in set:
				triple = k + k + k
				if triple in j:
					isvalid = False
			
			#checks if contained in a 'maximal' patterns
			for k in maximal:
				if nonContigSubstring(k, j):
					isvalid = False
			
			#checks if it contains an obstruction 
			for k in obstructions:
				if nonContigSubstring(j, k):
					isvalid = False 
			
			if isvalid:
				print(j)

def generateCodes(set, length):
	if length == 0:
		return ['']
	return [s + c for s in generateCodes(set, length-1) for c in set]
	
def nonContigSubstring(bigstring, substring):
	pos = 0
	for x in bigstring:
		if pos < len(substring) and x ==  substring[pos]:
			pos += 1
	return pos == len(substring)

#generateNotContained(['311331'],['11331331131133'], ['1', '3'], 15)

#0 specified columns
#generateNotContained(['213', '132', '231', '312', '2112', '2332', '311331'],['112232233','332212211','332232233','112212211', '11331331131133'], ['1','2', '3'], 15)

#generateNotContained(['311331'],[], ['1', '3'], 15)

#1 specified column on the left, in theory 0xx can be simplified to 0x, see below
#generateNotContained(['213', '132', '231', '312', '2112', '2332', '311331', '1331', '12', '32'],['331131133','2212211','2232233'], ['1','2', '3'], 12)

#1 specified column on the right, in theory, we could reduce xx0 to x0, but the code doesn't work with that
#generateNotContained(['213', '132', '231', '312', '2112', '2332', '311331', '21', '23', '3113'],['1122', '3322', '113313311', '331331133'], ['1','2', '3'], 12)


#1 specified column on each end, still not doing reductions like 0xx -> 0x, or xx0 -> x0

#generateNotContained(['213', '132', '231', '312', '2112', '2332', '311331', '21', '23', '3113', '13', '12', '32', '1331'],['22', '11', '3311'], ['1','2', '3'], 12)


#2 fully specified columns on the left
#generateNotContained(['213', '132', '231', '312', '2112', '2332', '311331', '1331', '12', '32', '31', '21'],['1133', '33', '2232233'], ['1','2', '3'], 12)


#2 fully specified columns on the right 
generateNotContained(['213', '132', '231', '312', '2112', '2332', '311331', '21', '23', '3113', '31', '32'],['33', '22', '1122', '1133'], ['1','2', '3'], 12)
