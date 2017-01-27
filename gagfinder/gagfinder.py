#!/usr/bin/python

import time # for timing functions
start_time = time.time()

#################################
# Step 0: imports and functions #
#################################

print "Importing modules...",

import re # for converting chemical formulae to dictionaries and vice versa
import sys # for getting user inputs
import argparse # for getting user inputs
import sqlite3 as sq # for accessing the database
debug = False # variable for debugging

if debug:
	import timeit
	start = timeit.default_timer()

# import pymzml
try:
	import pymzml
except:
	print "You need to install the pymzML module to use this script. Please install and try again."
	sys.exit()

# import numpy
try:
	import numpy as np
except:
	print "You need to install the numpy module to use this script. Please install and try again."
	sys.exit()

# import brainpy
try:
	from brainpy import isotopic_variants as iv
except:
	print "You need to install the brainpy module to use this script. Please install and try again."
	sys.exit()

# get local packages
#from species  import *

# function for converting dictionary to chemical formula or composition
def dict2fmla (dt, type):
	fs = '' # formula string
	
	# check whether chemical formula or composition
	if type == 'formula':
		letters = 'CHONS'
	elif type == 'composition':
		letters = 'DUXNAS'
	else:
		print "Incorrect type entered. Please enter either 'formula' or 'composition'."
		sys.exit()
	
	for let in letters:
		if let in dt:
			if dt[let] > 0:
				if dt[let] > 1:
					fs += let + str(dt[let])
				else:
					fs += let
	
	# return
	return fs

# function for converting chemical formula or composition to dictionary
def fmla2dict (fm, type):
	# check whether chemical formula or composition
	if type == 'formula':
		dt = {'C':0, 'H':0, 'O':0, 'N':0, 'S':0}
	elif type == 'composition':
		dt = {'D':0, 'U':0, 'X':0, 'N':0, 'A':0, 'S':0}
	else:
		print "Incorrect type entered. Please enter either 'formula' or 'composition'."
		sys.exit()
	
	parts = re.findall(r'([A-Z][a-z]*)(\d*)', fm.upper()) # split formula by symbol
	for q in parts:
		if q[0] not in dt: # invalid symbol entered
			print "Invalid chemical formula entered."
			sys.exit()
		else:
			if q[1] == '': # only one of this atom
				dt[q[0]] += 1
			else:
				dt[q[0]] += int(q[1])
	
	# return
	return dt

################################
# Step 1: check user arguments #
################################

print "Done!"
print "Checking user arguments...",

# initiate parser
parser = argparse.ArgumentParser(description='Find peaks in a GAG tandem mass spectrum.')

# add arguments
parser.add_argument('-c', required=True, help='GAG class')
parser.add_argument('-i', required=True, help='Input mzML file')
parser.add_argument('-o', required=True, help='Output file')
parser.add_argument('-m', type=float, required=True, help='Precursor m/z value')
parser.add_argument('-z', type=int, required=True, help='Precursor charge')
parser.add_argument('-r', required=False, help='Reducing end derivatization')
parser.add_argument('-n', type=int, required=False, help='Number of top results to return')
parser.add_argument('-p', type=float, required=False, help='Percentage of top results to return')

# parse arguments
args = parser.parse_args()

# get arguments into proper variables
gClass = args.c
dFile  = args.i
oFile  = args.o
pre_mz = args.m
pre_z  = args.z
fmla   = args.r
top_n  = args.n
top_p  = args.p

# get values ready for reducing end derivatization
df    = {'C':0, 'H':0, 'O':0, 'N':0, 'S':0}
atoms = ['C','H','O','N','S']

# check to make sure a proper GAG class was added
if gClass not in ['HS', 'CS', 'KS']:
	print "You must denote a GAG class, either HS, CS, or KS. Try 'python gagfinder.py --help'"

# parse the reducing end derivatization formula
if fmla:
	parts = re.findall(r'([A-Z][a-z]*)(\d*)', fmla.upper()) # split formula by symbol
	for q in parts:
		if q[0] not in atoms: # invalid symbol entered
			print "Invalid chemical formula entered. Please enter only CHONS. Try 'python gagfinder.py --help'"
		else:
			if q[1] == '': # only one of this atom
				df[q[0]] += 1
			else:
				df[q[0]] += int(q[1])

# check to make sure user didn't input both a top n and a top percentile
if top_n and top_p:
	print "You must select either a number or a percentile to return. Try 'python gagfinder.py -h'"
	sys.exit()

# check to make sure user did input at least one of top n or top percentile
if not top_n and not top_p:
	print "You must select a number or a percentile to return. Try 'python gagfinder.py -h'"
	sys.exit()

# print the system arguments back out to the user
if debug:
	print "class: %s" % (gClass)
	print "mzML file: %s" % (dFile)
	print "output file: %s" % (oFile)
	print "precursor m/z: %f" % (pre_mz)
	print "precursor charge: %i" % (pre_z)
	if 'top_n' in globals():
		print "top %i fragments requested" % (top_n)
	else:
		print "top %.2f%% percent of fragments requested" % (top_p)
	
	formula = ''
	for key in df:
		val = df[key]
		if val > 0:
			formula += key
			if val > 1:
				formula += str(val)
	print "atoms in reducing end derivatization: %s" % (formula)

#############################
# Step 2: Get data in order #
#############################

print "Done!"
print "Getting data in order...",

# from user, under construction
data = pymzml.run.Reader(dFile) # get mzML file into object

# UNDER CONSTRUCTION
# iterate through spectrums in mzML file and add them all together
#for spectrum in data:
	#s += spectrum
	#print len(spectrum.mz), len(spectrum.i)

s = pymzml.spec.Spectrum(measuredPrecision=2e-5) # initialize new Spectrum object
for t in data:
	s = t

s = s.removeNoise('median') # remove peaks below noise threshold

# use % max as intensity
#s = s/max(s.i)
#s = s*100

# get minimum and maximum m/z values
min_mz = min(s.mz)
max_mz = max(s.mz)

print "Minimum m/z: " + str(min_mz)
print "Maximum m/z: " + str(max_mz)
print "Range: " + str(max_mz - min_mz)

# get derivatization weight/formula
wt = {'C':12.0, 'H':1.00782504, 'O':15.99491461956, 'N':14.0030740048, 'S':31.97207100}
dw = 0
for q in df:
	dw += df[q] * wt[q]

# get precursor mass and charge (from user, under construction)
#pre_mz   = 375.7306
#pre_z    = 4
pre_mass = (pre_mz*abs(pre_z)) - (pre_z*1.00782504)

# pick a proper class number
if gClass == 'HS':
	cNum = 3
elif gClass == 'CS':
	cNum = 1
else:
	cNum = 4

# connect to GAGfragDB
conn = sq.connect('../lib/GAGfragDB.db')
c    = conn.cursor()

######################################
# Step 3: Find precursor composition #
######################################

print "Done!"
print "Finding precursor composition...",

# get precursor info
c.execute('''SELECT   cpm.id, f.value, p.value
             FROM     Precursors p, ClassPrecursorMap cpm, Formulae f
             WHERE    p.id = cpm.pId 
             AND      f.id = p.fmId
             AND      cpm.cId = ?
             ORDER BY ABS(f.monoMass - ?) ASC
             LIMIT 1;''', (cNum, pre_mass-dw))
row = c.fetchone()

# place precursor info into variables
id    = row[0]
pFmla = row[1]
pComp = row[2]

#########################################################
# Step 4: Get reducing end/non-reducing end information #
#########################################################

print "Done!"
print "Determining reducing end/non-reducing end information...",

# convert composition into a dictionary
pDict = fmla2dict(pComp, 'composition')
n_pre = pDict['D'] + pDict['U'] + pDict['X'] + pDict['N']

# check if dHexA exists (has to be NR end)
if pDict['D'] > 0:
	NR = 'D'
	
	if pDict['U'] == pDict['N']: # we know that the reducing end is HexA because (HexA+dHexA > HexN)
		RE = 'U'
	else: # we know that the reducing end is HexN because (HexA+dHexA == HexN)
		RE = 'N'
else:
	if cNum == 4: # KS
		if pDict['N'] > pDict['X']:
			NR = 'N'
			RE = 'N'
		elif pDict['N'] < pDict['X']:
			NR = 'X'
			RE = 'X'
		else:
			NR = '?'
			RE = '?'
	else: # HS or CS
		if pDict['N'] > pDict['U']: # we know that both ends are HexN
			NR = 'N'
			RE = 'N'
		elif pDict['N'] < pDict['U']: # we know that both ends are HexA
			NR = 'U'
			RE = 'U'
		else: # we cannot know which end is which just yet
			NR = '?'
			RE = '?'

######################################################################
# Step 5: Get potential fragment ions based on precursor composition #
######################################################################

print "Done!"
print "Retrieving all potential fragment ions for precursor with composition " + pComp + "...",

# generate all IDs
all_IDs   = {}
all_forms = {}

# get all child fragments
c.execute('''SELECT   fr.value, fm.value
             FROM     ChildFragments cf, Formulae fm, Fragments fr, Precursors p, ClassPrecursorMap cp
             WHERE    cf.frId = fr.id
             AND      cf.cpId = cp.id 
             AND      cp.pId = p.id
             AND      fr.fmId = fm.id
             AND      cf.cpId = ?;''', (id,))

# go through the rows
for row in c.fetchall():
	# get formula and convert it to dictionary
	fComp  = row[0]
	fFmla  = row[1]
	fcDict = fmla2dict(fComp.rsplit('+', 1)[0], 'composition')
	ffDict = fmla2dict(fFmla, 'formula')
	
	# get number of FULL monosaccharides in fragment
	n_frag = fcDict['D'] + fcDict['U'] + fcDict['X'] + fcDict['N']
	
	# determine if cross-ring cleavage
	if '+' in fComp:
		xr = 1
	else:
		xr = 0
	
	# determine if reducing end info is needed
	if dw == 0: # no reducing end derivatization
		reFrag = 0
	elif fcDict['D'] == 1: # presence of dHexA means non-reducing end fragment
		reFrag = 0
	elif 'NR' in fComp: # non-reducing end fragment
		reFrag = 0
	else: # easy reasons to reject reducing end all failed
		if RE == '?': # don't know reducing end monosaccharide, so doesn't matter
			reFrag = 1
		else: # do know reducing end monosaccharide, need to dig further
			if n_frag % 2 == 0: # equal number of Hex/HexA and HexN
				if RE+'RE' in fComp: # correct monosaccharide to cut into
					reFrag = 1
				elif RE == 'U' and 'DRE' in fComp: # reducing end HexA
					reFrag = 1
				elif not xr: # not x-ring...since equal number, will be RE
					reFrag = 1
				else: # incorrect monosaccharide to cut into
					reFrag = 0
			else: # unequal number of Hex/HexA and HexN
				sm = n_frag/2     # smaller number
				bg = (n_frag/2)+1 # bigger number
				
				if fcDict[RE] == bg: # if uneven number, RE must be larger number
					reFrag = 1
				else:
					reFrag = 0
	
	## ready to test fragments
	for i in range(2-xr): # determine how much water to lose
		for j in range(reFrag+1): # determine whether reducing end fragment or not
			for k in range(3): # cycle through potential hydrogen losses
				thDict = dict(ffDict) # copy directory of fragment
				
				# H2O loss
				thDict['H'] -= 2*i
				thDict['O'] -= i
				
				# reducing end
				for q in df:
					thDict[q] += j*df[q]
				
				# hydrogen loss
				thDict['H'] -= k
				
				# turn new formula dict into string
				thFmla = dict2fmla(thDict, 'formula')
				
				# append alterations onto fComp
				thComp = fComp
				
				# add RE info
				if j == 1:
					thComp += '+RE'
				
				# add water loss info
				if i == 1:
					thComp += '-H2O'
				elif i == 2:
					thComp += -'2H2O'
				
				# add hydrogen loss info
				if k == 1:
					thComp += '-H'
				elif k == 2:
					thComp += '-2H'
				
				# check if this formula already exists
				if thFmla not in all_forms:
					for z in range((pre_z+1), 0):
						dist                 = iv(thDict, charge=z)
						all_IDs[(thFmla, z)] = dist
					
					all_forms[thFmla] = [thComp]
				else:
					all_forms[thFmla].append(thComp)

print "Done!"
print "Scoring all potential fragment ions for precursor with composition " + pComp + "...",

# score each fragment
found_IDs  = {}
wherecache = {}
mono_int   = {}

# cycle through fragments
for j in all_IDs:
	TID = []
	EID = []
	
	kill = False
	
	for peak in all_IDs[j]:
		if peak.mz < s.mz[0] or peak.mz > s.mz[len(s.mz)-1]: # out of the m/z range
			kill = True
			break
		
		if peak.mz in wherecache:
			EID.append(wherecache[peak.mz])
		else:
			near = s.hasPeak(peak.mz)
			
			val = 0
			for k in near:
				if k[1] > val:
					val = k[1]
			
			if val == 0:
				val = 1e-6
			EID.append(val)
			wherecache[peak.mz] = val
		
		TID.append(peak.intensity)
	
	if not kill: # make sure ID is all in m/z range
		EI = np.array(EID)/sum(EID)
		TI = np.array(TID)
		g  = 2 * np.sum(EI * np.log(EI/TI))
		
		found_IDs[j] = g
		mono_int[j]  = EID[0]

print "Tested " + str(len(found_IDs)) + " fragments"
print str(len(all_IDs)) + " total fragments"

# rank the fragments by G-score
rank = sorted(found_IDs.items(), key=lambda x: x[1], reverse=False)

# get top hits
if top_n is not None: # user wants top n hits
	top = rank[:top_n]
else: # user wants top percentile hits
	pct = top_p/100.0
	ct  = int(pct * len(found_IDs))
	
	top = rank[:ct]

print "Done!"
print "Printing output to file...",

# write to file
f  = open(oFile, 'w')
f.write("m/z\tIntensity\tCharge\tFragments\tG-score\n")

for q in top:
	out = str(all_IDs[q[0]][0].mz)+'\t'+str(mono_int[q[0]])+'\t'+str(all_IDs[q[0]][0].charge)+'\t'
	for ions in all_forms[q[0][0]]:
		out += ions + ', '
	
	out = out[:-2] + '\t'
	out += str(found_IDs[q[0]]) + '\n'
	f.write(out)

f.close()

f = open(dFile[:-5]+'.txt', 'w')
for q in found_IDs:
	f.write(str(found_IDs[q]) + '\n')

f.close()

print "Done!"
print "Finished!"
print time.time() - start_time