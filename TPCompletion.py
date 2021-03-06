import random
import math
import numpy as np
import itertools
import re 
import mystic #don't forget to cite this if the code ends up being relevant for the paper

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
						isPartialTP = False
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
			if det < 0.0000001:
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
def guessCompletableFilledIn(arr, trials, debug):
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
	print("LOWER BOUNDS:")	
	for x,y in lowerBounds.items():
		print(x,y)
		
	print("UPPER BOUNDS:")	
	for x,y in upperBounds.items():
		print(x,y)	
	#a random search 
	#calculates a general upper bound, ignoring the (1,1) and (m,n) positions
#	GeneralUpperBound = np.amax(arr)*10 #this is worth experimenting with
	arrminuscorners = arr.copy()
	arrminuscorners[0, 0] = np.nan
	arrminuscorners[-1, -1] = np.nan
	GeneralUpperBoundNoCorners = np.nanmax(arrminuscorners)**2 + 2
	GeneralLowerBound = .0001  # probably not good to go too low, it'll clash with the determinant bound for '0'
	foundComp = False
	for i in range(trials):
		for entry in unspecList:
			rowPos = entry[0]
			colPos = entry[1]
			lowerbound = GeneralLowerBound
			upperbound = GeneralUpperBoundNoCorners
			if (rowPos, colPos) in lowerBounds:
				lowerbound = lowerBounds[(rowPos, colPos)]
				
			if (rowPos, colPos) in upperBounds:
				upperbound = upperBounds[(rowPos, colPos)]
			if upperbound == GeneralUpperBoundNoCorners and upperbound < lowerbound and not lowerbound == GeneralLowerBound: #sometimes this happens with really small numbers, the lower bound gets moved up really high
				upperbound = 2*lowerbound 
			arr[rowPos][colPos] = random.uniform(lowerbound, upperbound)
			
			if rowPos ==0 and colPos ==0:
				arr[rowPos][colPos] = 9223372036854775807 #these can just get really big
			if rowPos == len(arr) - 1 and colPos == len(arr[0])-1:
				arr[rowPos][colPos] = 9223372036854775807 #these can just get really big
		isTP = checkTP(arr, False)
		if isTP:
			if debug:
				print("Completion found")
				print(arr)
			foundComp = True 
			break 
	if debug:
		print("Random search finished.")
	if (not foundComp) and debug:
		print("A random search suggests this is not a completable partial matrix")
	
	return foundComp
	#TODO: even better guessing, find immediate contradiction, although I don't think that will happen much, if at all
	
#takes in a pattern in the form of a (0,1)-matrix, 0 being unspecified, 1 being unspecified, and attempts to guess if the pattern is completable 	
def guessCompletableNotFilledIn(arr):
	arr = arr.astype('float64') #because we usually only input integers, without this line, python casts the random numbers below to integers too
	unspecList = []
	specList = []
	trialsFilledIn = 10000
	isTP = True
	for i in range(len(arr)):
		for j in range(len(arr[0])):
			if arr[i][j] == 0:
				unspecList.append([i,j])
			else:
				specList.append([i,j])
	filledIn = arr
	for i in range(len(arr)):
		for j in range(len(arr[0])):
			filledIn[i][j] = 0

	#first, fill in with 1s and then just increment if 1s are not possible 
	count = 0 
	hardCap = 1000
	for i in range(len(arr)):
		for j in range(len(arr[0])):
			temp = [i,j]
			if temp in specList:
				filledIn[i][j] = 1
				
				while(not checkPartialTP(filledIn, False)) and count < hardCap:
					filledIn[i][j] = filledIn[i][j]+1
					count = count + 1
	if count < hardCap:
		isTP = guessCompletableFilledIn(filledIn, trialsFilledIn, False)
	if not isTP:
		print("This is a obstruction!")
		print("Here is the data for the specified entries")
		print(filledIn)
		return False
	
	#now, random searching 
	trials = 1000
	hardCap = 100
	randomLowerBound = 0
	randomUpperBound = 10 #this was at 100 originally, and the code struggled with completing these matrices. It still struggles here honestly
	for k in range(trials):
		count = 0
		for i in range(len(arr)):
			for j in range(len(arr[0])):
				filledIn[i][j] = 0
		for i in range(len(arr)):
			for j in range(len(arr[0])):
				temp = [i,j]
				if temp in specList:
					#we are free to do scaling to force the first column and row to be 1 (for the specified entries)
					if i == 0 or j == 0:
						filledIn[i][j] = 1 
					else:
						filledIn[i][j] = random.uniform(randomLowerBound, randomUpperBound)
						#this might infinite loop	
						while(not checkPartialTP(filledIn, False)) and count < hardCap:
							if count < hardCap -10:
								filledIn[i][j] = random.uniform(randomLowerBound,randomUpperBound)
							else: #if nothing seems to work, let it go really big
								indices = (-filledIn).argpartition(2, axis =None)[:2] #should find 2 largest values in the array
								x,y = np.unravel_index(indices, filledIn.shape) #gets the actual indices
								largestVal = filledIn[x,y]
								filledIn[i][j] = random.uniform(largestVal[0]*largestVal[1], 2*largestVal[0]*largestVal[1]/np.min(arr[np.nonzero(arr)])) #this is mostly for when the min value is < 1
							count = count + 1
		if count == hardCap:
			continue
		isTP = guessCompletableFilledIn(filledIn, trialsFilledIn, False)
		if not isTP: 
				print("This is a obstruction based on random data!")
				print("Here is the data for the specified entries")
				print(filledIn)
				return False
	print("Likely TP completable")
	return True 

#given a code like 2310, generates a 3 by n pattern ( 0,1 matrix ) such that the only unspecified entry in column i is in row code[i:i+1]-1 (0 corresponds to no unspecified entries)
#123012
def generateMatrixFromCode(code):
	arr = np.ones((3, len(code)))
	for i in range(len(code)):
		curr = int(code[i:i+1]) - 1
		if curr > -1:
			arr[curr][i] = 0
	return arr

def reflection(string):
#empty string which we will join our data to
	s = ''
	#newArray is an array that keeps track of the values in the string.
	stringArray = [int(i) for i in string]
	#newArray is generated so we don't override data when swapping items around
	newArray = [0  for i in range(len(stringArray))]
	#swaps the order of unspecified entries in the given pattern
	for i in range(len(stringArray)):
		if(stringArray[i] == 1):
			newArray[i] = 3
		elif (stringArray[i] == 2):
			newArray[i] = 2
		elif(stringArray[i] == 3):
			newArray[i] = 1
#reverses the columns in the pattern
	for entry in reversed(newArray):
		s += str(entry)
	return(s)
#generates codes for possible minimal 3 by (length) obstructions, avoiding codes in otherObstructions that are contiguously in them 
def generatePossibleMinimal(length, otherObstructions):
	codes = [''.join(i) for i in itertools.product("0123", repeat = length)]
	#removes matrices with 2 fully specified columns in the middle 
	codes = [i for i in codes if re.match(".00.", i) is None]
	#removes those which contain contiguous submatrices that are in otherObstructions 
	codes = [i for i in codes if not any(x in i for x in otherObstructions)]
	#removes those whose 180 degree reflections are also in otherObstructions
	codes = [i for i in codes if not str(reflection(i)) in otherObstructions]
	#expansions of obstructions, things we know are completable like 0-1 unspecified
#	print(codes)
	return codes


def findObstructions(length):
	codes = generatePossibleMinimal(3, [])
	obstructions = []
	for i in codes:
		arr = generateMatrixFromCode(i)
		isComp = guessCompletableNotFilledIn(arr)
		if not isComp:
			obstructions.append(arr)
	
	print("These are the suspected obstructions")
	for i in obstructions:
		print(i)
			

def tester(): 
	#np.set_printoptions(precision=1)
	#print(np.array([1.124125]))
#	test = np.array([[3.46514042,0, 7.10280744], [0, 0.24462701, 0], [8.80677369, 5.07579122, 0]])
	test = np.array([[1, 1, 1], [1, 0, 0], [0, 9.13371718,  9.35585903]]) #at the very least, the matrices we care about will have many specified entries, hence providing many bounds
	testPattern = np.array([[0,1, 1], [1, 1, 1], [1, 1, 0]])
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
	guessCompletableFilledIn(test, 10000, True)
#	guessCompletableNotFilledIn(testPattern)
#	findObstructions(3)
	
#todo scale 1s which seems to let our lower bound be 1, cases for bound?, don't check cases that are obviously obstructions, get more 1s scaled, use bounds with ab > c or ab < c (a,b unspecified, c constant)
tester()