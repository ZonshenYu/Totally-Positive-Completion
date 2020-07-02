import random
import math
import numpy as np
import itertools

#checks partial TP by treating 0 values as unspecified
#this just checks every single possible minor, so it is not efficient 
def checkPartialTP(arr, debug):
	arr = arr.astype('float64')
	numrows = len(arr)
	numcols = len(arr[0])
	minsize = min(numrows, numcols)
	minorsize = 1
	isPartialTP = True
	while minorsize <= minsize:
		rowSubs = list(itertools.combinations( (range(numrows)), minorsize))
		colSubs = list(itertools.combinations( (range(numcols)), minorsize))
		
		for rows in rowSubs:
			for cols in colSubs:
				submat = arr[np.ix_(rows, cols)]
				if (submat > 0).all():
					det = np.linalg.det(submat)
					if det < 0.000001: #rounding errors means a little leeway is needed, it should be fine as long as we don't input super small numbers
						isPartialTP = False;
						if debug:
							print("Negative minor with rows: ", [x+1 for x in rows], " and columns ", [x+1 for x in cols])
							
		minorsize = minorsize + 1
	return isPartialTP

#checks TP using all initial minors
def checkTP(arr, debug):
	arr = arr.astype('float64')
	numrows = len(arr)
	numcols = len(arr[0])
	minsize = min(numrows, numcols)
	minorsize = 2
	contigStart = 0
	contigEnd = minorsize-1
	isTP = True 
	if (arr < 0).any():
		return False
	
	while minorsize <= minsize:
		#check contiguous rows
		while contigEnd < numrows:
			submat = arr[np.ix_(range(contigStart, contigEnd+1),range(0, minorsize))]
			det = np.linalg.det(submat)
			if det < 0.0000001:
				isTP = False
				if debug:
					print("Negative initial minor with rows: ", contigStart+1, " through ", contigEnd+1, " and columns 1 through ", minorsize)
			contigStart = contigStart + 1
			contigEnd = contigEnd + 1
			
			
		contigStart = 0
		contigEnd = minorsize-1
		
		#check contiguous column
		while contigEnd < numcols:
			submat = arr[np.ix_(range(0, minorsize),range(contigStart, contigEnd+1))]
			det = np.linalg.det(submat)
			if det < 0.000001:
				isTP = False
				if debug:
					print("Negative initial minor with columns: ", contigStart+1, " through ", contigEnd+1, " and rows 1 through ", minorsize)
			contigStart = contigStart + 1
			contigEnd = contigEnd + 1
		
		minorsize = minorsize+1 
		contigStart = 0
		contigEnd = minorsize - 1 
		
	return isTP

#given a partial TP matrix with data for the specified entries, tests for completablility. 
#zeroes are still used for designating unspecified entries
def guessCompletableFilledIn(arr, trials):
	arr = arr.astype('float64') #because we usually only input integers, without this line, python casts the random numbers below to integers too
	unspecList = []
	for i in range(len(arr)):
		for j in range(len(arr[0])):
			if arr[i][j] == 0:
				unspecList.append( (i,j) ) #the 3rd entry will be the lower bound, the 4th the upper bound of values (2^63-1)  0, 9223372036854775807
	
	#use 2 by 2 minors with 3 specified entries to derive better bounds per unspecified entry
	lowerBounds = {} 
	upperBounds = {}
	
	for i in range(len(arr)):
		for j in range(len(arr[0])):
			if arr[i][j] == 0:
				continue 
			for k in range(j+1, len(arr[0])):
				if arr[i][k] == 0:
					continue
				for l in range(i+1, len(arr)):
				#i is the upper row, l the lower row, j the left column, k the right column
					if arr[l][j] > 0 and arr[l][k] == 0:
						#looks like: 
						#x x 
						#x ? 
						newlowerbound = arr[i][k]*arr[l][j]/arr[i][j]
						if (l,k) in lowerBounds:
							lowerBounds[(l,k)] = max(lowerBounds[(l,k)], newlowerbound)
						else:
							lowerBounds[(l,k)] = newlowerbound
					
					if arr[l][k] > 0 and arr[l][j] == 0:
						#looks like: 
						#x x 
						#? x 
						newupperbound = arr[i][j]*arr[l][k]/arr[i][k]
						if (l,j) in upperBounds:
							upperBounds[(l,j)] = min(upperBounds[(l,j)], newupperbound)
						else:
							upperBounds[(l,j)] = newupperbound 
				for l in range(i):
				#i is the lower row, l the upper row, j the left column, k the right column
					if arr[l][j] == 0 and arr[l][k] > 0:
						#looks like: 
						#? x 
						#x x 
						newlowerbound = arr[i][j]*arr[l][k]/arr[i][k]
						if (l,j) in lowerBounds:
							lowerBounds[(l,j)] = max(lowerBounds[(l,j)], newlowerbound)
						else:
							lowerBounds[(l,j)] = newlowerbound
					
					if arr[l][k] == 0 and arr[l][j] > 0:
						#looks like: 
						#x ? 
						#x x 
						newupperbound = arr[i][k]*arr[l][j]/arr[i][j]
						if (l,k) in upperBounds:
							upperBounds[(l,k)] = min(upperBounds[(l,k)], newupperbound)
						else:
							upperBounds[(l,k)] = newupperbound 
	
	#tester for lower and upper bounds, note that these are 0 indexed
#	print("LOWER BOUNDS:")	
#	for x,y in lowerBounds.items():
#		print(x,y)
		
#	print("UPPER BOUNDS:")	
#	for x,y in upperBounds.items():
#		print(x,y)	
	#a random search 
	GeneralUpperBound = np.amax(arr)*2 #this is worth experimenting with
	GeneralLowerBound = .001 #probably not good to go too low, it'll clash with the determinant bound for '0'
	foundComp = False 
	for i in range(trials):
		for entry in unspecList:
			rowPos = entry[0]
			colPos = entry[1]
			lowerbound = GeneralLowerBound
			upperbound = GeneralUpperBound 
			if (rowPos, colPos) in lowerBounds:
				lowerbound = lowerBounds[(rowPos, colPos)]
			
			if (rowPos, colPos) in upperBounds:
				upperbound = upperBounds[(rowPos, colPos)]
			arr[rowPos][colPos] = random.uniform(lowerbound, upperbound)
		isTP = checkTP(arr, False)
		if isTP:
			print("Completion found")
			print(arr)
			foundComp = True 
			break 
	print("Random search finished.")
	if not foundComp:
		print("A random search suggests this is not a completable partial matrix")
	
	#TODO: even better guessing, find immediate contradiction, although I don't think that will happen much, if at all
	
#takes in a pattern in the form of a (0,1)-matrix, 0 being unspecified, 1 being unspecified, and attempts to guess if the pattern is completable 	
#not complete
def guessCompletableNotFilledIn(arr):
	arr = arr.astype('float64') #because we usually only input integers, without this line, python casts the random numbers below to integers too
	unspecList = []
	specList = []
	for i in range(len(arr)):
		for j in range(len(arr[0])):
			if arr[i][j] == 0:
				unspecList.append([i,j])
			else:
				specList.append([i,j])
	
def tester(): 
	#np.set_printoptions(precision=1)
	#print(np.array([1.124125]))
	test = np.array([[9999999999,0, 2], [4, 1, 0], [2, 1, 0]])
#	test = np.array([[1, 0, 0, 1], [1, 1, 1, 2], [0, 1, 2, 0]])
	print(test)
	#print(np.linalg.det(test))
	#isTP = checkTP(test, True)

	isPartialTP = checkPartialTP(test, True)
	if isPartialTP:
		print("The matrix is Partial TP")
	else:
		print("The matrix is not Partial TP")
	
#	isTP = checkTP(test, True)
#	if isTP:
#		print("The matrix is TP")
#	else:
#		print("The matrix is not TP")
	guessCompletableFilledIn(test, 1000)

tester()