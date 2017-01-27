#!/usr/bin/python

import re # for converting chemical formulae to dictionaries and vice versa
import sqlite3 as sq # for accessing the database
from species import *

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

def choose_iter(elements, length):
	for i in xrange(len(elements)):
		if length == 1:
			yield (elements[i],)
		else:
			for next in choose_iter(elements[i+1:len(elements)], length-1):
				yield (elements[i],) + next

def choose(l, k):
	return list(choose_iter(l, k))

# sulfate positions
xloc = {}

# HexA
xloc['U'] = {}

# initialize dictionaries for non-reducing end and reducing end
xloc['U']['NR'] = {}
xloc['U']['RE'] = {}

# add modification possibilities
xloc['U']['NR']['0,2'] = []
xloc['U']['NR']['1,5'] = [2]
xloc['U']['NR']['2,4'] = []
xloc['U']['NR']['3,5'] = []
xloc['U']['NR']['0,3'] = []
xloc['U']['NR']['1,4'] = [2]
xloc['U']['NR']['2,5'] = []
xloc['U']['RE']['0,2'] = [2]
xloc['U']['RE']['1,5'] = []
xloc['U']['RE']['2,4'] = [2]
xloc['U']['RE']['3,5'] = [2]
xloc['U']['RE']['0,3'] = [2]
xloc['U']['RE']['1,4'] = []
xloc['U']['RE']['2,5'] = [2]

# N
xloc['N'] = {}

# initialize dictionaries for non-reducing end and reducing end
xloc['N']['NR'] = {}
xloc['N']['RE'] = {}

# add modification possibilities
xloc['N']['NR']['0,2'] = [3,6]
xloc['N']['NR']['1,5'] = [2,3,6]
xloc['N']['NR']['2,4'] = [2]
xloc['N']['NR']['3,5'] = [6]
xloc['N']['NR']['0,3'] = [6]
xloc['N']['NR']['1,4'] = [2,3]
xloc['N']['NR']['2,5'] = [3,6]
xloc['N']['RE']['0,2'] = [2]
xloc['N']['RE']['1,5'] = []
xloc['N']['RE']['2,4'] = [3,6]
xloc['N']['RE']['3,5'] = [2,3]
xloc['N']['RE']['0,3'] = [2,3]
xloc['N']['RE']['1,4'] = [6]
xloc['N']['RE']['2,5'] = [2]

# load G-scores
f = open('/Users/jdhogan5/Desktop/test.txt', 'r')
f.readline()

# dictionary to store G-scores
gdict = {}

# go through lines
for line in f.readlines():
	cells = line.strip('\n').split('\t')
	
	G   = float(cells[4])
	fgs = cells[3].split(', ')
	
	if fgs in gdict:
		if G < gdict[fgs]:
			gdict[fgs] = G
	else:
		gdict[fgs] = G

gsc = []

for key in gdict:
	gsc.append((key, G))

dt  = np.dtype([('frags', 'object'), ('G', 'float')])
gsc = np.array(gsc, dtype=dt)

# normalize G-scores
gsc['G'] = (np.mean(gsc['G'])-gsc['G'])/np.std(gsc['G'])

mol      = 'U2N3S8'
backbone = {'1':'N','2':'U','3':'N','4':'U','5':'N'}
all_mods = ['1-2','1-3','1-6','2-2','3-2','3-3','3-6','4-2','5-2','5-3','5-6']

nA = 0
nS = 8
loss = 0

frags = {}

# function for walking
def walk(scores, f_set):
	run    = 0
	maxrun = 0.0
	
	allfrags = []
	for one in f_set:
		allfrags.append(one)
	
	for row in scores:
		frags = row[0]
		
		ct = 0
		
		for f in frags:
			if f in allfrags:
				ct += 1
		
		if ct > 0:
			#run += np.sqrt((13053-len(f_set))/len(f_set))
			run += 1
		else:
			#run -= np.sqrt(len(f_set)/(13053-len(f_set)))
			run -= 1
		
		if run > maxrun:
			maxrun = run
	
	return maxrun

sf_perms  = choose(all_mods, nS)
en_scores = []

for gag in sf_perms:
	mod_spot = {'1':{'2':0,'3':0,'6':0}, '2':{'2':0}, '3':{'2':0,'3':0,'6':0}, '4':{'2':0}, '5':{'2':0,'3':0,'6':0}}
	for spot in gag:
		unit, loc = spot.split('-')
		mod_spot[unit][loc] = 1
	
	nn = 1
	monos = ['1','2','3','4','5']
	
	frags[gag] = []
	
	while nn < 5:
		for i in range(len(monos)-nn+1):
			j = monos[i:i+nn]
			
			# get the complete monosaccharides
			fml = {'D':0, 'U':0, 'X':0, 'N':0, 'A':0, 'S':0}
			for k in j:
				fml[backbone[k]] += 1
				
				for p in mod_spot[k]:
					fml['S'] += mod_spot[k][p]
			
			## glycosidic fragments
			# before considering sulfate losses
			sf = dict2fmla(fml, 'composition')
			if i != (len(monos)-nn):
				if sf not in frags[gag]:
					frags[gag].append(sf)
			else:
				sf1 = sf + '+RE'
				if (sf1, z) not in frags[gag]:
					frags[gag].append(sf1)
				
			sf2 = sf+'-H2O'
			if i != (len(monos)-nn):
				if sf2 not in frags[gag]:
					frags[gag].append(sf2)
			else:
				sf3 = sf + '+RE-H2O'
				if sf3 not in frags[gag]:
					frags[gag].append(sf3)
			
			# consider sulfate losses
			fml1 = dict(fml)
			for s in range(1, (loss+1)):
				fml1['S'] = max(0, fml1['S']-1)
			
				sf = dict2fmla(fml1, 'composition')
				if i != (len(monos)-nn):
					if sf not in frags[gag]:
						frags[gag].append(sf)
				else:
					sf1 = sf + '+RE'
					if sf1 not in frags[gag]:
						frags[gag].append(sf1)
				
				sf2 = sf+'-H2O'
				if i != (len(monos)-nn):
					if sf2 not in frags[gag]:
						frags[gag].append(sf2)
				else:
					sf3 = sf + '+RE-H2O'
					if sf3 not in frags[gag]:
						frags[gag].append(sf3)
			
			## crossring fragments
			if monos[i] == '1': # non-reducing end
				addto = monos[nn]
				msch  = backbone[addto]
				
				for x in xloc[msch]['NR']:
					fml2 = dict(fml)
					
					for y in xloc[msch]['NR'][x]:
						fml2['S'] += mod_spot[addto][str(y)]
					
					# before considering sulfate losses
					sf = dict2fmla(fml2, 'composition') + '+' + msch + 'NR' + x
					if sf not in frags[gag]:
						frags[gag].append(sf)
					
					# consider sulfate losses
					for s in range(1, (loss+1)):
						fml2['S'] = max(0, fml2['S']-1)
						
						sf = dict2fmla(fml2, 'composition') + '+' + msch + 'NR' + x
						if sf not in frags[gag]:
							frags[gag].append(sf)
			elif i == (len(monos) - nn): # reducing end
				addto = monos[i-1]
				msch  = backbone[addto]
				
				for x in xloc[msch]['RE']:
					fml2 = dict(fml)
					
					for y in xloc[msch]['RE'][x]:
						fml2['S'] += mod_spot[addto][str(y)]
					
					# before considering sulfate losses
					sf = dict2fmla(fml2, 'composition') + '+' + msch + 'RE' + x + '+RE'
					if sf not in frags[gag]:
						frags[gag].append(sf)
					
					# consider sulfate losses
					for s in range(1, (loss+1)):
						fml2['S'] = max(0, fml2['S']-1)
						
						sf = dict2fmla(fml2, 'composition') + '+' + msch + 'RE' + x + '+RE'
						if sf not in frags[gag]:
							frags[gag].append(sf)
		
		nn += 1
	
	en_scores.append((gag, walk(gsc, frags[gag])))

en_scores = np.array(en_scores, dtype=[('GAG', 'object'), ('ES', 'float')])
en_scores = np.sort(en_scores, order='ES')
#en_scores[0:len(np.where(en_scores['ES'] == max(en_scores['ES']))[0])]
en_scores