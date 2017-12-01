#!/usr/bin/python

#################################
# Step 0: imports and functions #
#################################

print "Importing modules...",

import time # for timing functions
start_time = time.time()

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
	from brainpy import mass_dict
	pt    = mass_dict.nist_mass
	elems = pt.keys()
except:
	print "You need to install the brainpy module to use this script. Please install and try again."
	sys.exit()

# get local packages
from species import *

### CLASSES AND FUNCTIONS ###

# from Josh
class TheoreticalIsotopicPattern(object):
    def __init__(self, base_tid, truncated_tid=None):
        if truncated_tid is None:
            truncated_tid = base_tid
        self.base_tid = base_tid
        self.truncated_tid = truncated_tid
	
    def __getitem__(self, i):
        return self.truncated_tid[i]
	
    def __iter__(self):
        return iter(self.truncated_tid)
	
    def __len__(self):
        return len(self.truncated_tid)
	
    def get(self, i):
        return self.truncated_tid[i]
	
    def clone(self):
        return self.__class__([p.clone() for p in self.base_tid],
                              [p.clone() for p in self.truncated_tid])
	
    @property
    def monoisotopic_mz(self):
        return self.base_tid[0].mz
	
    def shift(self, mz, truncated=True):
        first_peak = self.base_tid[0]
        for peak in self.base_tid[1:]:
            delta = peak.mz - first_peak.mz
            peak.mz = mz + delta
        first_peak.mz = mz
		
        if truncated:
            first_peak = self.truncated_tid[0]
            for peak in self.truncated_tid[1:]:
                delta = peak.mz - first_peak.mz
                peak.mz = mz + delta
            first_peak.mz = mz
		
        return self
	
    def truncate_after(self, truncate_after=0.0):
        cumsum = 0
        result = []
        for peak in self.base_tid:
            cumsum += peak.intensity
            result.append(peak)
            if cumsum >= truncate_after:
                break
        for peak in result:
            peak.intensity *= 1. / cumsum
        self.truncated_tid = result
        return self
	
    def ignore_below(self, ignore_below=0.0):
        total = 0
        kept_tid = []
        for i, p in enumerate(self.truncated_tid):
            if p.intensity < ignore_below and i > 1:
                continue
            else:
                total += p.intensity
                kept_tid.append(p.clone())
        for p in kept_tid:
            p.intensity /= total
        self.truncated_tid = kept_tid
        return self
	
    def scale(self, experimental_distribution, method='sum'):
        if method == 'sum':
            total_abundance = sum(
                p.intensity for p in experimental_distribution)
            for peak in self:
                peak.intensity *= total_abundance
        elif method == 'max':
            i, peak = max(enumerate(self),
                          key=lambda x: x[1].intensity)
            scale_factor = experimental_distribution[
                i].intensity / peak.intensity
            for peak in self:
                peak.intensity *= scale_factor
        return self
	
    def __repr__(self):
        return "TheoreticalIsotopicPattern(%0.4f, charge=%d, (%s))" % (
            self.base_tid[0].mz,
            self.base_tid[0].charge,
            ', '.join("%0.3f" % p.intensity for p in self.truncated_tid))

# function for converting dictionary to chemical formula or composition
def dict2fmla (dt, type):
	fs = '' # formula string
	
	# check whether chemical formula or composition
	if type == 'formula':
		symbols = elems
	elif type == 'composition':
		symbols = ['D', 'U', 'X', 'N', 'A', 'S']
	else:
		print "Incorrect type entered. Please enter either 'formula' or 'composition'."
		sys.exit()
	
	for sym in symbols:
		if sym in dt:
			if dt[sym] > 0:
				if dt[sym] > 1:
					fs += sym + str(dt[sym])
				else:
					fs += sym
	
	# return
	return fs

# function for converting chemical formula or composition to dictionary
def fmla2dict (fm, type):
	# check whether chemical formula or composition
	if type == 'formula':
		symbols = elems
		dt = {'C':0, 'H':0, 'O':0, 'N':0, 'S':0} #, 'Na':0, 'K':0, 'Li':0, 'Mg':0, 'Ca':0}
	elif type == 'composition':
		symbols = ['D', 'U', 'X', 'N', 'A', 'S']
		dt = {'D':0, 'U':0, 'X':0, 'N':0, 'A':0, 'S':0}
	else:
		print "Incorrect type entered. Please enter either 'formula' or 'composition'."
		sys.exit()
	
	parts = re.findall(r'([A-Z][a-z]*)(\d*)', fm.upper()) # split formula by symbol
	for q in parts:
		if q[0] not in dt: # invalid symbol entered
			if q[0] not in symbols:
				print "Invalid chemical formula entered."
				sys.exit()
			else:
				dt[q[0]] = 0
		
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
parser = argparse.ArgumentParser(description='Find isotopic clusters in GAG tandem mass spectra.')

# add arguments
parser.add_argument('-c', required=True, help='GAG class (required)')
parser.add_argument('-i', required=True, help='Input mzML file (required)')
parser.add_argument('-r', required=False, help='Reducing end derivatization (optional)')
parser.add_argument('-n', type=int, required=False, help='Number of top results to return (either)')
parser.add_argument('-p', type=float, required=False, help='Percentage of top results to return (either)')
parser.add_argument('-a', required=False, help='Which metal is adducted (optional)')
parser.add_argument('-t', type=int, required=False, help='Number of adducted metals (required if -a TRUE)')
parser.add_argument('-g', required=False, help='Reagent used (optional)')
parser.add_argument('-m', type=float, required=False, help='Precursor m/z (optional, but must be in mzML file)')
parser.add_argument('-z', type=int, required=False, help='Precursor charge (optional, but must be in mzML file)')
parser.add_argument('-s', type=int, required=False, help='Number of sulfate losses to consider (optional, default 0)')
parser.add_argument('-x', required=False, help='Has the noise already been removed? (y/n)')

# parse arguments
args = parser.parse_args()

# get arguments into proper variables
gClass  = args.c
dFile   = args.i
fmla    = args.r
top_n   = args.n
top_p   = args.p
adduct  = args.a
nMetal  = args.t
reag    = args.g
pre_mz  = args.m
pre_z   = args.z
s_loss  = args.s
removed = args.x

# get values ready for reducing end derivatization and reagent
df    = {'C':0, 'H':0, 'O':0, 'N':0, 'S':0}
rf    = {'C':0, 'H':0, 'O':0, 'N':0, 'S':0}
atoms = ['C','H','O','N','S']

# check to make sure a proper GAG class was added
if gClass not in ['HS', 'CS', 'KS']:
	print "You must denote a GAG class, either HS, CS, or KS. Try 'python gagfinder.py --help'"
	sys.exit()

# pick a proper class number
if gClass == 'HS':
	cNum = 3
elif gClass == 'CS':
	cNum = 1
else:
	cNum = 4

# parse the reducing end derivatization formula
if fmla:
	parts = re.findall(r'([A-Z][a-z]*)(\d*)', fmla.upper()) # split formula by symbol
	for q in parts:
		if q[0] not in atoms: # invalid symbol entered
			print "Invalid chemical formula entered. Please enter only CHONS. Try 'python gagfinder.py --help'"
			sys.exit()
		else:
			if q[1] == '': # only one of this atom
				df[q[0]] += 1
			else:
				df[q[0]] += int(q[1])

# parse the reagent formula
if reag:
	parts = re.findall(r'([A-Z][a-z]*)(\d*)', reag.upper()) # split formula by symbol
	for q in parts:
		if q[0] not in atoms: # invalid symbol entered
			print "Invalid chemical formula entered. Please enter only CHONS. Try 'python gagfinder.py --help'"
			sys.exit()
		else:
			if q[1] == '': # only one of this atom
				rf[q[0]] += 1
			else:
				rf[q[0]] += int(q[1])

# get derivatization weight
wt = {'C':  12.0,
      'H':  1.0078250322,
      'O':  15.994914620,
      'N':  14.003074004,
      'S':  31.972071174,
      'Na': 22.98976928,
      'K':  38.96370649,
      'Li': 7.01600344,
      'Ca': 39.9625909,
      'Mg': 23.98504170}
dw = 0
for q in df:
	dw += df[q] * wt[q]

# get reagent weight
rw = 0
for q in rf:
	rw += rf[q] * wt[q]

# check to make sure that metal adduct is appropriate
if adduct:
	if adduct not in ['Na', 'K', 'Li', 'Ca', 'Mg']:
		print "\nInvalid metal adduct entered. Only Na, K, Li, Ca, and Mg are accepted. Try 'python gagfinder.py --help'"
		sys.exit()
	
	if nMetal is None:
		nMetal = 0
	elif nMetal < 1:
		print "\nYou must enter a positive integer for the number of adducted metals. Try 'python gagfinder.py --help'"
		sys.exit()

# check to make sure user didn't input both a top n and a top percentile
if top_n and top_p:
	print "You must select either a number or a percentile to return, not both. Try 'python gagfinder.py -h'"
	sys.exit()

# check to see if the user wants to consider sulfate loss
if not s_loss:
	s_loss = 0

# check to make sure user did input at least one of top n or top percentile
if not top_n and not top_p:
	print "You must select a number or a percentile to return. Try 'python gagfinder.py -h'"
	sys.exit()

# check to see if the user removed the noise to start with
if not removed:
	removed = False
else:
	# check to make sure the entry was correct
	if removed == 'y' or removed == 'n':
		if removed == 'y':
			removed = True
		else:
			removed = False
	else:
		print "You must enter either 'y' or 'n' for whether the noise has been removed or not. Try 'python gagfinder.py -h'"
		sys.exit()

# print the system arguments back out to the user
if debug:
	print "class: %s" % (gClass)
	print "mzML file: %s" % (dFile)
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

############################################
# Step 2: Load mzml file and connect to DB #
############################################

print "Done!"
print "Loading mzml file and connecting to GAGfragDB...",

# from user, under construction
data = pymzml.run.Reader(dFile, MSn_Precision=5e-06) # get mzML file into object
s    = pymzml.spec.Spectrum(measuredPrecision=20e-06) # initialize new Spectrum object

# connect to GAGfragDB
conn = sq.connect('../lib/GAGfragDB.db')
c    = conn.cursor()

######################
# Step 3a: Sum scans #
######################

print "Done!"
print "Summing scans...",

counter = 0.0

'''for t in data:
	counter += 1
	if t['ms level'] == 2:
		t2 = t.deRef()
		if 'total ion current' in t.keys():
			s += t2/t2['total ion current']
		else:
			s += t2/sum(t2.i)'''

for t in data:
	counter += 1
	if t['ms level'] == 2:
		if 'total ion current' in t.keys():
			s += t/t['total ion current']
		else:
			s += t/sum(t.i)

s = s/counter

#mn = sum([i for mz, i in s.peaks])/float(len(s.peaks))
md = s._median([i for mz, i in s.peaks])

# noise not already removed by the user
if not removed:
	s = s.removeNoise(mode='median')

# get minimum and maximum m/z values
ends   = s.extremeValues('mz')
min_mz = ends[0]
max_mz = ends[1]

if debug:
	print "Minimum m/z: " + str(min_mz)
	print "Maximum m/z: " + str(max_mz)
	print "Range: " + str(max_mz - min_mz)

#######################################
# Step 3b: Find precursor composition #
#######################################

print "Done!"
print "Finding precursor composition...",

# calculate precursor mass
pre_mass = (pre_mz*abs(pre_z)) - (pre_z*wt['H'])

# get proper precursor mass to test
if not nMetal:
	nMetal = 0
	test_mass = pre_mass - dw
else:
	test_mass = pre_mass - dw - (nMetal*wt[adduct]) + (nMetal*wt['H'])
	
# get precursor info
c.execute('''SELECT   cpm.id, f.value, p.value, f.monoMass
   	         FROM     Precursors p, ClassPrecursorMap cpm, Formulae f
       	     WHERE    p.id = cpm.pId 
           	 AND      f.id = p.fmId
             AND      cpm.cId = ?
   	         ORDER BY ABS(f.monoMass - ?) ASC
       	     LIMIT 1;''', (cNum, test_mass))
row = c.fetchone()

# place precursor info into variables
id    = row[0]
pFmla = row[1]
pComp = row[2]
tMass = row[3] # theoretical mass, for precision calculation

#########################################################
# Step 3c: Get reducing end/non-reducing end information #
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

#######################################################################
# Step 3d: Get potential fragment ions based on precursor composition #
#######################################################################

print "Done!"
print "Retrieving all potential fragment ions for precursor with composition " + pComp + "...",

# generate all IDs
all_IDs   = {}
all_forms = {}
all_comps = []

## look through precursor-based fragments
ffDict = fmla2dict(pFmla, 'formula')

# add metal adduct info
if adduct is not None:
	ffDict[adduct] = nMetal
	ffDict['H']   -= nMetal

# ready to test fragments
for i in range(2): # determine how much water to lose
	for j in range(2): # determine whether reducing end fragment or not
		for k in range(3): # cycle through potential hydrogen losses
			for l in range(pDict['U'] + pDict['D'] + 1): # cycle through potential CO2 losses
				for m in range(min(s_loss, ffDict['S'])+1): # cycle through potential SO3 losses
					for a in range(2): # cycle through reagent possibilities
						thDict = dict(ffDict) # copy directory of fragment
						
						# add reagent
						for q in rf:
							thDict[q] += a*rf[q]
						
						# H2O loss
						thDict['H'] -= 2*i
						thDict['O'] -= i
						
						# reducing end
						for q in df:
							thDict[q] += j*df[q]
						
						# hydrogen loss
						thDict['H'] -= k
						
						# CO2 loss
						thDict['C'] -= l
						thDict['O'] -= 2*l
						
						# SO3 loss
						thDict['S'] -= m
						thDict['O'] -= 3*m
						
						# turn new formula dict into string
						thFmla = dict2fmla(thDict, 'formula')
						
						# append alterations onto fComp
						thComp = 'M'
						
						# add RE info
						if j == 0:
							thComp += '-RE'
						
						# add water loss info
						if i == 1:
							thComp += '-H2O'
						elif i == 2:
							thComp += '-2H2O'
						
						# add hydrogen loss info
						if k == 1:
							thComp += '-H'
						elif k == 2:
							thComp += '-2H'
						
						# add CO2 loss info
						if l == 1:
							thComp += '-CO2'
						elif l > 1:
							thComp += '-' + str(l) + 'CO2'
						
						# add SO3 loss info
						if m == 1:
							thComp += '-SO3'
						elif m > 1:
							thComp += '-' + str(m) + 'SO3'
						
						# add reagent info
						if a == 1:
							thComp += '+A'
						
						# check if this formula already exists
						if thFmla not in all_forms:
							for z in range((pre_z+1), 0):
								dist                 = TheoreticalIsotopicPattern(iv(thDict, charge=z)).truncate_after(0.95)
								all_IDs[(thFmla, z)] = dist
							
							all_forms[thFmla] = [thComp]
						else:
							all_forms[thFmla].append(thComp)

## get all child fragments
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
	
	if '+' in fComp: # cross-ring fragment
		fcDict  = fmla2dict(fComp.rsplit('+', 1)[0], 'composition')
		xr_info = fComp.rsplit('+', 1)[1]
	else: # glycosidic fragment
		fcDict  = fmla2dict(fComp, 'composition')
		xr_info = ''
	
	ffDict = fmla2dict(fFmla, 'formula')
	
	# get number of FULL monosaccharides in fragment
	n_frag = fcDict['D'] + fcDict['U'] + fcDict['X'] + fcDict['N']
	
	# calculate values for searching adducted metals
	atleast = nMetal - (pDict['S'] + pDict['U'] + pDict['D']) + fcDict['S'] + fcDict['U'] + fcDict['D']
	atmost  = fcDict['S'] + fcDict['U'] + fcDict['D']
	
	# calculate number of CO2 here
	nCO2 = fcDict['U'] + fcDict['D']
	
	# determine if cross-ring cleavage
	if '+' in fComp:
		xr    = 1
		
		# get info about the cross-ring monosaccharide
		x_det = fComp.rsplit('+')[1]
		x_ms  = x_det[0]
		x_end = x_det[1:3]
		x_clv = x_det[3:len(x_det)]
		
		# transcribe from single letters to three-letter codes
		if x_ms == 'D':
			x_ms = 'dHexA'
		elif x_ms == 'U':
			x_ms = 'HexA'
		elif x_ms == 'X':
			x_ms = 'Hex'
		else:
			x_ms = 'HexN'
		
		atleast += xmod[gClass][x_ms][x_end][x_clv]['COOH']
		atmost  += xmod[gClass][x_ms][x_end][x_clv]['COOH']
		
		nCO2 += xmod[gClass][x_ms][x_end][x_clv]['COOH']
	else:
		xr = 0
	
	# determine if reducing end info is needed
	if fcDict['D'] == 1: # presence of dHexA means non-reducing end fragment
		reFrag = [0]
	elif 'NR' in fComp and 'NRE' not in fComp: # non-reducing end fragment
		reFrag = [0]
	else: # easy reasons to reject reducing end all failed
		if RE == '?': # don't know reducing end monosaccharide...this will be tricky
			if not xr: # not cross-ring fragment...could be anything
				reFrag = [0,1]
			else: # cross-ring fragment
				if n_frag == n_pre - 1: # has to be RE fragment
					reFrag = [1]
				else:
					reFrag = [0,1]
		else: # do know reducing end monosaccharide
			if n_frag % 2 == 0: # equal number of Hex/HexA and HexN
				if RE+'RE' in fComp: # correct monosaccharide to cut into for cross-ring fragment on reducing end
					if n_frag < n_pre - 2: # doesn't HAVE to be on the end
						reFrag = [0,1]
					else: # HAS to be on the end
						reFrag = [1]
				elif RE == 'U' and 'DRE' in fComp: # reducing end HexA
					if n_frag == n_pre - 1: # HAS to be on the end
						reFrag = [1]
					else:
						reFrag = [0,1]
				elif not xr: # not x-ring...since equal number, can be either RE or not
					reFrag = [0,1]
				else: # incorrect monosaccharide to cut into
					reFrag = [0]
			else: # unequal number of Hex/HexA and HexN
				sm = n_frag/2     # smaller number
				bg = (n_frag/2)+1 # bigger number
				
				if fcDict[RE] == bg: # if uneven number, RE must be larger number to possibly be a RE fragment
					if xr: # cross-ring fragment
						if n_frag < n_pre - 2: # doesn't HAVE to be on the end
							reFrag = [0,1]
						else: # HAS to be on the end
							reFrag = [1]
					else: # glycosidic fragment
						if n_frag < n_pre - 1: # doesn't HAVE to be on the end
							reFrag = [0,1]
						else: # HAS to be on the end
							reFrag = [1]
				else:
					reFrag = [0]
	
	## ready to test fragments
	for i in range(2-xr): # determine how much water to lose
		for j in reFrag: # determine whether reducing end fragment or not
			for k in range(-2, 1): # cycle through potential hydrogen losses/gains
				for c in range(min(nCO2, 0) + 1):
					for m in range(max(0, atleast), (min(nMetal, atmost)+1)): # cycle through potential metal adduction
						for so3 in range(min(s_loss, ffDict['S'])+1): # cycle through potential SO3 losses
							thffDict = dict(ffDict) # copy directory of fragment formula
							thfcDict = dict(fcDict) # copy directory of fragment composition
							
							# H2O loss
							thffDict['H'] -= 2*i
							thffDict['O'] -= i
							
							# reducing end
							if j == 1:
								for q in df:
									thffDict[q] += j*df[q]
							
							# hydrogen loss/gain
							thffDict['H'] += k
							
							# CO2 loss
							thffDict['C'] -= c
							thffDict['O'] -= c*2
							
							# metal adduction
							if adduct is not None:
								thffDict[adduct] = m
								thffDict['H']   -= m
							
							# sulfate loss
							thffDict['S'] -= so3
							thffDict['O'] -= 3*so3
							
							thfcDict['S'] -= so3
							
							# turn new formula dict into string
							thffFmla = dict2fmla(thffDict, 'formula')
							
							# append alterations onto composition
							thComp = dict2fmla(thfcDict, 'composition')
							if xr_info != '':
								thComp += '+' + xr_info
							
							# add RE info
							if j == 1:
								thComp += '+RE'
							
							# add water loss info
							if i == 1:
								thComp += '-H2O'
							elif i == 2:
								thComp += '-2H2O'
							
							# add hydrogen gain/loss info
							sign = ''
							if abs(k) == 1:
								if k < 0:
									sign = '-'
								else:
									sign = '+'
								
								thComp += sign + 'H'
							elif abs(k) == 2:
								if k < 0:
									sign = '-'
								else:
									sign = '+'
								
								thComp += sign + '2H'
							
							# add CO2 loss info
							if c == 1:
								thComp += '-CO2'
							elif c > 1:
								thComp += '-' + str(c) + 'CO2'
							
							# add metal adduction info
							if m == 1:
								thComp += '+' + adduct
							elif m > 1:
								thComp += '+' + str(m) + adduct
							
							# check to see if this one has already been searched
							if thComp not in all_comps:
								all_comps.append(thComp)
							
								# check if this formula already exists
								if thffFmla not in all_forms:
									for z in range((pre_z+1), 0):
										dist                   = TheoreticalIsotopicPattern(iv(thffDict, charge=z)).truncate_after(0.95)
										all_IDs[(thffFmla, z)] = dist
									
									all_forms[thffFmla] = [thComp]
								else:
									all_forms[thffFmla].append(thComp)

##############################################
# Step 3e: Score each potential fragment ion #
##############################################

print "Done!"
print "Scoring all potential fragment ions for precursor with composition " + pComp + "...",

# dictionaries to store info
found_IDs = {}
mono_int  = {}
mono_mz   = {}
errors    = {}
'''wherecache = {}

# cycle through fragments
for j in all_IDs:
	TID = []
	EID = []
	
	kill = False # boolean that states if the searched m/z is in the spectrum
	
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
'''

# cycle through fragments
for j in all_IDs:
	TID = []
	EID = []
	
	kill = False # variable to determine whether to kill or not
	
	# get first peak to determine precision
	peak = all_IDs[j][0]
	if peak.mz < s.mz[0] or peak.mz > s.mz[len(s.mz)-1]: # out of the m/z range
		continue
	else: # first peak is in m/z range
		near = s.hasPeak(peak.mz)
		
		val  = 0
		for k in near:
			if k[1] > val:
				loc = k[0]
				val = k[1]
		
		if val == 0:
			kill = True
		else:
			diff = peak.mz - loc
			ppm  = 1e6 * ((loc - peak.mz)/peak.mz)
		
		TID.append(peak.intensity)
		EID.append(val)
	
	if kill:
		continue
	
	# cycle through remaining peaks
	for peak in all_IDs[j][1:len(all_IDs[j])]:
		if peak.mz < s.mz[0] or peak.mz > s.mz[len(s.mz)-1]: # out of the m/z range
			kill = True
			break
			
		near = s.hasPeak(peak.mz)
		
		val  = 0
		prec = 1000000
		for k in near:
			if peak.mz - k[0] - diff < prec:
				val = k[1]
		
		if val == 0:
			val = 1e-100
		
		EID.append(val)
		TID.append(peak.intensity)
	
	if not kill:
		EI = np.array(EID)/sum(EID)
		TI = np.array(TID)
		
		g  = 2 * np.sum(EI * np.log(EI/TI))
		
		found_IDs[j] = g
		mono_int[j]  = EID[0]
		mono_mz[j]   = loc
		errors[j]    = ppm

print "Done!"
print "Tested " + str(len(found_IDs)) + " out of " + str(len(all_IDs)) + " fragments"

# rank the fragments by G-score
rank = sorted(found_IDs.items(), key=lambda x: x[1], reverse=False)

# get top hits
if top_n is not None: # user wants top n hits
	top = rank[:top_n]
else: # user wants top percentile hits
	pct = top_p/100.0
	ct  = int(pct * len(found_IDs))
	
	top = rank[:ct]

print "Printing output to file...",

# write to file
oFile = dFile[:-5] + '.txt'
f  = open(oFile, 'w')
f.write("m/z\tIntensity\tCharge\tFragments\tG-score\tError (ppm)\n")

for q in top:
	out = str(mono_mz[q[0]])+'\t'+str(mono_int[q[0]])+'\t'+str(all_IDs[q[0]][0].charge)+'\t'
	for ions in all_forms[q[0][0]]:
		out += ions + '; '
	
	out = out[:-2] + '\t'
	out += str(found_IDs[q[0]]) + '\t' + str(errors[q[0]]) + '\n'
	f.write(out)

f.close()

print "Finished!"
print time.time() - start_time