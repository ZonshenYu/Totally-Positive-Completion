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
				unspecList.append([i,j])
				
	#first, just a random search 
	upperbound = np.amax(arr)*2 #this is worth experimenthing with
	lowerbound = .001 #probably not good to go too low, it'll clash with the determinant bound for '0'
	foundComp = False 
	for i in range(trials):
		for entry in unspecList:
			arr[entry[0]][entry[1]] = random.uniform(lowerbound, upperbound)
		isTP = checkTP(arr, False)
		if isTP:
			print("Completion found")
			print(arr)
			foundComp = True 
			break 
	print("Random search finished.")
	if not foundComp:
		print("A random search suggests this is not a completable partial matrix")
	
def tester(): 
	#np.set_printoptions(precision=1)
	#print(np.array([1.124125]))
	test = np.array([[1,1,1,0], [0,1,2, 1],[1, 2,0,3]])
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
	guessCompletableFilledIn(test, 10000)

tester()