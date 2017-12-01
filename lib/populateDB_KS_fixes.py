#!/usr/bin/python

# imports
import sqlite3 as sq
import re
import sys
from species import *

# function for making dictionary of formula/composition
def makeDict (type, D=0, U=0, X=0, N=0, A=0, S=0):
	if type == 'formula':
		dc      = {'C':0, 'H':0, 'O':0, 'N':0, 'S':0}
		dc['C'] = D*fm['dHexA']['C'] + U*fm['HexA']['C'] + X*fm['Hex']['C'] + N*fm['HexN']['C'] + A*fm['Ac']['C']
		dc['H'] = D*fm['dHexA']['H'] + U*fm['HexA']['H'] + X*fm['Hex']['H'] + N*fm['HexN']['H'] + A*fm['Ac']['H'] - ((D+U+X+N-1)*fm['H2O']['H'])
		dc['O'] = D*fm['dHexA']['O'] + U*fm['HexA']['O'] + X*fm['Hex']['O'] + N*fm['HexN']['O'] + A*fm['Ac']['O'] + S*fm['SO3']['O'] - ((D+U+X+N-1)*fm['H2O']['O'])
		dc['N'] = N*fm['HexN']['N']
		dc['S'] = S*fm['SO3']['S']
	else: # type is composition
		dc = {'D':D, 'U':U, 'X':X, 'N':N, 'A':A, 'S':S}
	
	# return
	return dc

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

# connect to GAGfragDB
conn = sq.connect('GAGfragDB.db')
c    = conn.cursor()

#################
# ADD FRAGMENTS #
#################

## add KS fragments
# get KS precursors
c.execute('''SELECT cpm.id, p.value
             FROM   ClassPrecursorMap cpm, Precursors p
             WHERE  cpm.pId = p.id
             AND    cpm.cId = 4;''')

# loop through precursors
for row in c.fetchall():
	# get map ID and formula
	mapId = row[0]
	fmla  = row[1]
	
	# only use ids that haven't been added yet
	#if mapId < 3389:
		#continue
	
	# convert formula to dictionary
	cp = fmla2dict(fmla, 'composition')
	
	# get precursor info
	n_pre   = cp['X'] + cp['N']
	sys.stdout.write(fmla + '\n')
	
	# initialize variable that denotes length of fragment
	nn = 1
	
	########################
	# GLYCOSIDIC FRAGMENTS #
	########################
	
	# loop through all fragment sizes smaller than full molecule
	while nn < n_pre:
		if nn % 2 == 0: # equal number of Hex and HexN
			eq = nn/2
			
			j = eq # each HexN will be acetylated
			for k in range(max(0, cp['S'] - (cp['N'] + cp['X']) + nn),(min(cp['S'],nn)+1)): # go through SO3 possibilities
				# chemical formula
				f = makeDict('formula', D=0, U=0, X=eq, N=eq, A=j, S=k)
				
				# convert formula to string
				fstr = (dict2fmla(f, 'formula'),)
				
				# get formula's id value
				c.execute('''SELECT id
				             FROM   Formulae
				             WHERE  value = ?;''', fstr)
				fmId = c.fetchone()[0]
				
				# composition
				cp2 = makeDict('composition', D=0, U=0, X=eq, N=eq, A=j, S=k)
				
				# convert composition to string
				cstr = (dict2fmla(cp2, 'composition'),)
				
				# check if this fragment already exists in the database
				c.execute('''SELECT COUNT(*)
				             FROM   Fragments
				             WHERE  value = ?;''', cstr)
				
				# this fragment does not exist in the database
				if c.fetchone()[0] == 0:
					# insert into Fragments table
					c.execute('''INSERT INTO Fragments (value, fmId)
					             VALUES (?,?);''', (cstr[0], fmId))
					conn.commit()
				
				# get the fragment id
				c.execute('''SELECT id
				             FROM   Fragments
				             WHERE  value = ?;''', cstr)
				frId = c.fetchone()[0]
				
				# add fragment as a child
				c.execute('''INSERT INTO ChildFragments (cpId, frId)
				             VALUES (?,?);''', (mapId, frId))
				conn.commit()
		else: # unequal number of Hex and HexN
			sm = nn/2     # smaller number
			bg = (nn/2)+1 # bigger number
			
			# HexN > Hex
			j = bg # each HexN will be acetylated
			for k in range(max(0, cp['S'] - (cp['N'] + cp['X']) + nn),(min(cp['S'],nn)+1)): # go through SO3 possibilities
				# chemical formula
				f = makeDict('formula', D=0, U=0, X=sm, N=bg, A=j, S=k)
				
				# convert formula to string
				fstr = (dict2fmla(f, 'formula'),)
				
				# get formula's id value
				c.execute('''SELECT id
				             FROM   Formulae
				             WHERE  value = ?;''', fstr)
				fmId = c.fetchone()[0]
				
				# composition
				cp2 = makeDict('composition', D=0, U=0, X=sm, N=bg, A=j, S=k)
				
				# convert composition to string
				cstr = (dict2fmla(cp2, 'composition'),)
				
				# check if this fragment already exists in the database
				c.execute('''SELECT COUNT(*)
				             FROM   Fragments
				             WHERE  value = ?;''', cstr)
				
				# this fragment does not exist in the database
				if c.fetchone()[0] == 0:
					# insert into Fragments table
					c.execute('''INSERT INTO Fragments (value, fmId)
					             VALUES (?,?);''', (cstr[0], fmId))
					conn.commit()
				
				# get the fragment id
				c.execute('''SELECT id
				             FROM   Fragments
				             WHERE  value = ?;''', cstr)
				frId = c.fetchone()[0]
				
				# add fragment as a child
				c.execute('''INSERT INTO ChildFragments (cpId, frId)
				             VALUES (?,?);''', (mapId, frId))
				conn.commit()
			
			# Hex > HexN
			j = sm  # each HexN will be acetylated
			for k in range(max(0, cp['S'] - (cp['N'] + cp['X']) + nn),(min(cp['S'],nn)+1)): # go through SO3 possibilities
				# chemical formula
				f = makeDict('formula', D=0, U=0, X=bg, N=sm, A=j, S=k)
				
				# convert formula to string
				fstr = (dict2fmla(f, 'formula'),)
				
				# get formula's id value
				c.execute('''SELECT id
				             FROM   Formulae
				             WHERE  value = ?;''', fstr)
				fmId = c.fetchone()[0]
				
				# composition
				cp2 = makeDict('composition', D=0, U=0, X=bg, N=sm, A=j, S=k)
				
				# convert composition to string
				cstr = (dict2fmla(cp2, 'composition'),)
				
				# check if this fragment already exists in the database
				c.execute('''SELECT COUNT(*)
				             FROM   Fragments
				             WHERE  value = ?;''', cstr)
				
				# this fragment does not exist in the database
				if c.fetchone()[0] == 0:
					# insert into Fragments table
					c.execute('''INSERT INTO Fragments (value, fmId)
					             VALUES (?,?);''', (cstr[0], fmId))
					conn.commit()
				
				# get the fragment id
				c.execute('''SELECT id
				             FROM   Fragments
				             WHERE  value = ?;''', cstr)
				frId = c.fetchone()[0]
				
				# add fragment as a child
				c.execute('''INSERT INTO ChildFragments (cpId, frId)
				             VALUES (?,?);''', (mapId, frId))
				conn.commit()
		
		nn += 1

	########################
	# CROSS-RING FRAGMENTS #
	########################
	nn = 1
	
	# loop through all fragment sizes smaller than full molecule
	while nn < n_pre:
		if nn % 2 == 0: # equal number of Hex and HexN in fragment
			eq = nn/2
			
			if n_pre % 2 == 0: # equal number of Hex and HexN in precursor
				for x in ['Hex', 'HexN']: # cycle through the monosaccharides to be added to
					# get letters for the formula
					if x == 'Hex':
						fval = 'X'
					else:
						fval = 'N'
					
					for y in ['NR', 'RE']: # cycle through the ends
						for z in xmod['KS'][x][y]: # cycle through the cross-ring cleavage possibilities
							ext = xmod['KS'][x][y][z]
							j = eq + ext['Ac'] # each HexN will be acetylated
							for k in range(max(0, cp['S'] - (cp['N'] + cp['X']) + nn + ext['SO3']), min(cp['S'], nn + ext['SO3'])+1): # go through SO3 possibilities
								# chemical formula
								f = makeDict('formula', D=0, U=0, X=eq, N=eq, A=j, S=k)
								
								# add necessary x-ring values
								for sym in 'CHONS':
									f[sym] += xfm[x][y][z][sym]
								
								# remove water
								f['H'] -= 2
								f['O'] -= 1
								
								# convert formula to string
								fstr = (dict2fmla(f, 'formula'),)
								
								# check if formula exists
								c.execute('''SELECT COUNT(*)
								             FROM   Formulae
								             WHERE  value = ?;''', fstr)
								
								# this formula does not exist in the database
								if c.fetchone()[0] == 0:
									# monoisotopic weight
									w = eq*wt['monoHex'] + eq*wt['monoHexN'] + xwt[x][y][z] + j*wt['monoAc'] + k*wt['monoSO3'] - (nn*wt['monoH2O'])
									
									# insert into Formulae table
									c.execute('''INSERT INTO Formulae (value, monoMass)
				            					 VALUES (?,?);''', (fstr[0], w))
									conn.commit()
								
								# get formula's id value
								c.execute('''SELECT id
								             FROM   Formulae
								             WHERE  value = ?;''', fstr)
								fmId = c.fetchone()[0]
								
								# composition
								cp2 = makeDict('composition', D=0, U=0, X=eq, N=eq, A=j, S=k)
								
								# convert composition to string
								cstr = (dict2fmla(cp2, 'composition') + "+" + fval + y + z,)
								
								# check if this fragment already exists in the database
								c.execute('''SELECT COUNT(*)
								             FROM   Fragments
								             WHERE  value = ?;''', cstr)
								
								# this fragment does not exist in the database
								if c.fetchone()[0] == 0:
									# insert into Fragments table
									c.execute('''INSERT INTO Fragments (value, fmId)
										         VALUES (?,?);''', (cstr[0], fmId))
									conn.commit()
								
								# get the fragment id
								c.execute('''SELECT id
						    			     FROM   Fragments
								             WHERE  value = ?;''', cstr)
								frId = c.fetchone()[0]
								
								# add fragment as a child
								c.execute('''INSERT INTO ChildFragments (cpId, frId)
						    			     VALUES (?,?);''', (mapId, frId))
								conn.commit()
			else: # unequal number of Hex and HexN in precursor
				# get ends
				if cp['X'] > cp['N']:
					x    = 'Hex'
					fval = 'X'
				else:
					x    = 'HexN'
					fval = 'N'
				
				# see if we're only looking at fragments that cut into ends or not
				if nn < (n_pre-1): # fragments don't cut into ends just yet
					for x in ['Hex', 'HexN']: # cycle through the monosaccharides to be added to
						# get letters for the formula
						if x == 'Hex':
							fval = 'X'
						else:
							fval = 'N'
							
						for y in ['NR', 'RE']: # cycle through the ends
							for z in xmod['KS'][x][y]: # cycle through the cross-ring cleavage possibilities
								ext = xmod['KS'][x][y][z]
								j = eq + ext['Ac'] # each HexN will be acetylated
								for k in range(max(0, cp['S'] - (cp['N'] + cp['X']) + nn + ext['SO3']), min(cp['S'], nn + ext['SO3'])+1): # go through SO3 possibilities
									# chemical formula
									f = makeDict('formula', D=0, U=0, X=eq, N=eq, A=j, S=k)
									
									# add necessary x-ring values
									for sym in 'CHONS':
										f[sym] += xfm[x][y][z][sym]
									
									# remove water
									f['H'] -= 2
									f['O'] -= 1
									
									# convert formula to string
									fstr = (dict2fmla(f, 'formula'),)
									
									# check if formula exists
									c.execute('''SELECT COUNT(*)
									             FROM   Formulae
									             WHERE  value = ?;''', fstr)
									
									# this formula does not exist in the database
									if c.fetchone()[0] == 0:
										# monoisotopic weight
										w = eq*wt['monoHex'] + eq*wt['monoHexN'] + xwt[x][y][z] + j*wt['monoAc'] + k*wt['monoSO3'] - (nn*wt['monoH2O'])
									
										# insert into Formulae table
										c.execute('''INSERT INTO Formulae (value, monoMass)
					            					 VALUES (?,?);''', (fstr[0], w))
										conn.commit()
									
									# get formula's id value
									c.execute('''SELECT id
									             FROM   Formulae
									             WHERE  value = ?;''', fstr)
									fmId = c.fetchone()[0]
									
									# composition
									cp2 = makeDict('composition', D=0, U=0, X=eq, N=eq, A=j, S=k)
									
									# convert composition to string
									cstr = (dict2fmla(cp2, 'composition') + "+" + fval + y + z,)
									
									# check if this fragment already exists in the database
									c.execute('''SELECT COUNT(*)
									             FROM   Fragments
									             WHERE  value = ?;''', cstr)
									
									# this fragment does not exist in the database
									if c.fetchone()[0] == 0:
										# insert into Fragments table
										c.execute('''INSERT INTO Fragments (value, fmId)
											         VALUES (?,?);''', (cstr[0], fmId))
										conn.commit()
									
									# get the fragment id
									c.execute('''SELECT id
						        			     FROM   Fragments
									             WHERE  value = ?;''', cstr)
									frId = c.fetchone()[0]
									
									# add fragment as a child
									c.execute('''INSERT INTO ChildFragments (cpId, frId)
							    			     VALUES (?,?);''', (mapId, frId))
									conn.commit()
				else: # fragments cut into ends
					# get letters for the formula
					if cp['X'] > cp['N']:
						x    = 'Hex'
						fval = 'X'
					else:
						x    = 'HexN'
						fval = 'N'
					
					for y in ['NR', 'RE']: # cycle through the ends
						for z in xmod['KS'][x][y]: # cycle through the cross-ring cleavage possibilities
							ext = xmod['KS'][x][y][z]
							j = eq + ext['Ac'] # each HexN will be acetylated
							for k in range(max(0, cp['S'] - (cp['N'] + cp['X']) + nn + ext['SO3']), min(cp['S'], nn + ext['SO3'])+1): # go through SO3 possibilities
								# chemical formula
								f = makeDict('formula', D=0, U=0, X=eq, N=eq, A=j, S=k)
								
								# add necessary x-ring values
								for sym in 'CHONS':
									f[sym] += xfm[x][y][z][sym]
								
								# remove water
								f['H'] -= 2
								f['O'] -= 1
								
								# convert formula to string
								fstr = (dict2fmla(f, 'formula'),)
								
								# check if formula exists
								c.execute('''SELECT COUNT(*)
								             FROM   Formulae
								             WHERE  value = ?;''', fstr)
								
								# this formula does not exist in the database
								if c.fetchone()[0] == 0:
									# monoisotopic weight
									w = eq*wt['monoHex'] + eq*wt['monoHexN'] + xwt[x][y][z] + j*wt['monoAc'] + k*wt['monoSO3'] - (nn*wt['monoH2O'])
									
									# insert into Formulae table
									c.execute('''INSERT INTO Formulae (value, monoMass)
				            					 VALUES (?,?);''', (fstr[0], w))
									conn.commit()
								
								# get formula's id value
								c.execute('''SELECT id
								             FROM   Formulae
								             WHERE  value = ?;''', fstr)
								fmId = c.fetchone()[0]
								
								# composition
								cp2 = makeDict('composition', D=0, U=0, X=eq, N=eq, A=j, S=k)
								
								# convert composition to string
								cstr = (dict2fmla(cp2, 'composition') + "+" + fval + y + z,)
								
								# check if this fragment already exists in the database
								c.execute('''SELECT COUNT(*)
								             FROM   Fragments
								             WHERE  value = ?;''', cstr)
								
								# this fragment does not exist in the database
								if c.fetchone()[0] == 0:
									# insert into Fragments table
									c.execute('''INSERT INTO Fragments (value, fmId)
										         VALUES (?,?);''', (cstr[0], fmId))
									conn.commit()
								
								# get the fragment id
								c.execute('''SELECT id
						    			     FROM   Fragments
								             WHERE  value = ?;''', cstr)
								frId = c.fetchone()[0]
								
								# add fragment as a child
								c.execute('''INSERT INTO ChildFragments (cpId, frId)
						    			     VALUES (?,?);''', (mapId, frId))
								conn.commit()
		else: # unequal number of Hex and HexN in fragment
			sm = nn/2     # smaller number
			bg = (nn/2)+1 # bigger number
			
			for mas in ['Hex', 'HexN']: # get which monosaccharide there is more of
				if mas == 'Hex':
					x    = 'HexN'
					fval = 'N'
					nX   = bg
					nN   = sm
				else:
					x    = 'Hex'
					fval = 'X'
					nX   = sm
					nN   = bg
				
				for y in ['NR', 'RE']: # cycle through the ends
					for z in xmod['KS'][x][y]: # cycle through the cross-ring cleavage possibilities
						ext = xmod['KS'][x][y][z]
						j = nN + ext['Ac'] # each HexN will be acetylated
						for k in range(max(0, cp['S'] - (cp['N'] + cp['X']) + nn + ext['SO3']), min(cp['S'], nn + ext['SO3'])+1): # go through SO3 possibilities
							# chemical formula
							f = makeDict('formula', D=0, U=0, X=nX, N=nN, A=j, S=k)
							
							# add necessary x-ring values
							for sym in 'CHONS':
								f[sym] += xfm[x][y][z][sym]
							
							# remove water
							f['H'] -= 2
							f['O'] -= 1
							
							# convert formula to string
							fstr = (dict2fmla(f, 'formula'),)
							
							# check if formula exists
							c.execute('''SELECT COUNT(*)
							             FROM   Formulae
							             WHERE  value = ?;''', fstr)
							
							# this formula does not exist in the database
							if c.fetchone()[0] == 0:
								# monoisotopic weight
								w = nX*wt['monoHex'] + nN*wt['monoHexN'] + xwt[x][y][z] + j*wt['monoAc'] + k*wt['monoSO3'] - (nn*wt['monoH2O'])
								
								# insert into Formulae table
								c.execute('''INSERT INTO Formulae (value, monoMass)
				        					 VALUES (?,?);''', (fstr[0], w))
								conn.commit()
							
							# get formula's id value
							c.execute('''SELECT id
							             FROM   Formulae
							             WHERE  value = ?;''', fstr)
							fmId = c.fetchone()[0]
							
							# composition
							cp2 = makeDict('composition', D=0, U=0, X=nX, N=nN, A=j, S=k)
							
							# convert composition to string
							cstr = (dict2fmla(cp2, 'composition') + "+" + fval + y + z,)
							
							# check if this fragment already exists in the database
							c.execute('''SELECT COUNT(*)
							             FROM   Fragments
							             WHERE  value = ?;''', cstr)
							
							# this fragment does not exist in the database
							if c.fetchone()[0] == 0:
								# insert into Fragments table
								c.execute('''INSERT INTO Fragments (value, fmId)
									         VALUES (?,?);''', (cstr[0], fmId))
								conn.commit()
							
							# get the fragment id
							c.execute('''SELECT id
									     FROM   Fragments
							             WHERE  value = ?;''', cstr)
							frId = c.fetchone()[0]
							
							# add fragment as a child
							c.execute('''INSERT INTO ChildFragments (cpId, frId)
					    			     VALUES (?,?);''', (mapId, frId))
							conn.commit()
		
		nn += 1
