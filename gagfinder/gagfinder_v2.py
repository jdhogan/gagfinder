#!/usr/bin/python

#################################
# Step 0: imports and functions #
#################################

print "Importing modules and functions...",

# initialize dictionaries for weights and formulae
wt = {}
fm = {}

# monoisotopic element weights
wt['monoH'] = 1.00782504
wt['monoC'] = 12.0
wt['monoO'] = 15.99491461956
wt['monoN'] = 14.0030740048
wt['monoS'] = 31.97207100

# monoisotopic compound weights
wt['monoHex']   = 6*wt['monoC'] + 12*wt['monoH'] + 6*wt['monoO']
wt['monoHexA']  = 6*wt['monoC'] + 10*wt['monoH'] + 7*wt['monoO']
wt['monoHexN']  = 6*wt['monoC'] + 13*wt['monoH'] + wt['monoN'] + 5*wt['monoO']
wt['monodHexA'] = 6*wt['monoC'] + 8*wt['monoH'] + 6*wt['monoO']
wt['monoH2O']   = 2*wt['monoH'] + wt['monoO']
wt['monoSO3']   = wt['monoS'] + 3*wt['monoO']
wt['monoAc']    = 2*wt['monoC'] + 2*wt['monoH'] + wt['monoO']

# formulae
fm['Hex']   = {'C':6, 'H':12, 'O':6, 'N':0, 'S':0}
fm['HexA']  = {'C':6, 'H':10, 'O':7, 'N':0, 'S':0}
fm['HexN']  = {'C':6, 'H':13, 'O':5, 'N':1, 'S':0}
fm['dHexA'] = {'C':6, 'H':8, 'O':6, 'N':0, 'S':0}
fm['H2O']   = {'C':0, 'H':2, 'O':1, 'N':0, 'S':0}
fm['SO3']   = {'C':0, 'H':0, 'O':3, 'N':0, 'S':1}
fm['Ac']    = {'C':2, 'H':2, 'O':1, 'N':0, 'S':0}

### modification locations
# initialize dictionary
modlocs       = {}
modlocs['HS'] = {'D':{}, 'U':{}, 'N':{}}
modlocs['CS'] = {'D':{}, 'U':{}, 'N':{}}
modlocs['KS'] = {'X':{}, 'N':{}}

# Acetyl
modlocs['HS']['N'] = {'Ac': [2]}
modlocs['CS']['N'] = {'Ac': [2]}
modlocs['KS']['N'] = {'Ac': [2]}

# Sulfate
modlocs['HS']['D']['SO3'] = [2]
modlocs['HS']['U']['SO3'] = [2]
modlocs['HS']['N']['SO3'] = [2,3,6]
modlocs['CS']['D']['SO3'] = [2]
modlocs['CS']['U']['SO3'] = [2]
modlocs['CS']['N']['SO3'] = [4,6]
modlocs['KS']['X']['SO3'] = [6]
modlocs['KS']['N']['SO3'] = [6]

### cross-ring fragments ###
# initialize dictionaries for weights and formulae
xwt  = {}
xfm  = {}
xmod = {}

### cross-ring formulae ###
## HexA
xfm['HexA'] = {}

# initialize dictionaries for non-reducing end and reducing end
xfm['HexA']['NR'] = {}
xfm['HexA']['RE'] = {}

# add fragments
xfm['HexA']['NR']['0,2'] = {'C':4, 'H':6, 'O':5, 'N':0, 'S':0}
xfm['HexA']['NR']['1,3'] = {'C':2, 'H':4, 'O':2, 'N':0, 'S':0}
xfm['HexA']['NR']['1,5'] = {'C':5, 'H':8, 'O':5, 'N':0, 'S':0}
xfm['HexA']['NR']['2,4'] = {'C':2, 'H':4, 'O':2, 'N':0, 'S':0}
xfm['HexA']['NR']['3,5'] = {'C':3, 'H':4, 'O':3, 'N':0, 'S':0}
xfm['HexA']['NR']['0,3'] = {'C':3, 'H':4, 'O':4, 'N':0, 'S':0}
xfm['HexA']['NR']['1,4'] = {'C':3, 'H':6, 'O':3, 'N':0, 'S':0}
xfm['HexA']['NR']['2,5'] = {'C':4, 'H':6, 'O':4, 'N':0, 'S':0}
xfm['HexA']['RE']['0,2'] = {'C':2, 'H':4, 'O':2, 'N':0, 'S':0}
xfm['HexA']['RE']['1,3'] = {'C':4, 'H':6, 'O':5, 'N':0, 'S':0}
xfm['HexA']['RE']['1,5'] = {'C':1, 'H':2, 'O':2, 'N':0, 'S':0}
xfm['HexA']['RE']['2,4'] = {'C':4, 'H':6, 'O':5, 'N':0, 'S':0}
xfm['HexA']['RE']['3,5'] = {'C':3, 'H':6, 'O':4, 'N':0, 'S':0}
xfm['HexA']['RE']['0,3'] = {'C':3, 'H':6, 'O':3, 'N':0, 'S':0}
xfm['HexA']['RE']['1,4'] = {'C':3, 'H':4, 'O':4, 'N':0, 'S':0}
xfm['HexA']['RE']['2,5'] = {'C':2, 'H':4, 'O':3, 'N':0, 'S':0}

## HexA
xfm['dHexA'] = {}

# initialize dictionaries for non-reducing end and reducing end
xfm['dHexA']['RE'] = {}

# add fragments
xfm['dHexA']['RE']['0,2'] = {'C':2, 'H':4, 'O':2, 'N':0, 'S':0}
xfm['dHexA']['RE']['1,3'] = {'C':4, 'H':4, 'O':4, 'N':0, 'S':0}
xfm['dHexA']['RE']['1,5'] = {'C':1, 'H':2, 'O':2, 'N':0, 'S':0}
xfm['dHexA']['RE']['2,4'] = {'C':4, 'H':5, 'O':5, 'N':0, 'S':0}
xfm['dHexA']['RE']['3,5'] = {'C':3, 'H':6, 'O':4, 'N':0, 'S':0}
xfm['dHexA']['RE']['0,3'] = {'C':3, 'H':6, 'O':3, 'N':0, 'S':0}
xfm['dHexA']['RE']['1,4'] = {'C':3, 'H':3, 'O':4, 'N':0, 'S':0}
xfm['dHexA']['RE']['2,5'] = {'C':2, 'H':4, 'O':3, 'N':0, 'S':0}

## HexN
xfm['HexN'] = {}

# initialize dictionaries for non-reducing end and reducing end
xfm['HexN']['NR'] = {}
xfm['HexN']['RE'] = {}

# add fragments
xfm['HexN']['NR']['0,2'] = {'C':4, 'H':8, 'O':4, 'N':0, 'S':0}
xfm['HexN']['NR']['1,3'] = {'C':2, 'H':5, 'O':1, 'N':1, 'S':0}
xfm['HexN']['NR']['1,5'] = {'C':5, 'H':11, 'O':3, 'N':1, 'S':0}
xfm['HexN']['NR']['2,4'] = {'C':2, 'H':4, 'O':2, 'N':0, 'S':0}
xfm['HexN']['NR']['3,5'] = {'C':3, 'H':6, 'O':2, 'N':0, 'S':0}
xfm['HexN']['NR']['0,3'] = {'C':3, 'H':7, 'O':2, 'N':1, 'S':0}
xfm['HexN']['NR']['1,4'] = {'C':3, 'H':7, 'O':2, 'N':1, 'S':0}
xfm['HexN']['NR']['2,5'] = {'C':4, 'H':5, 'O':3, 'N':0, 'S':0}
xfm['HexN']['RE']['0,2'] = {'C':2, 'H':5, 'O':1, 'N':1, 'S':0}
xfm['HexN']['RE']['1,3'] = {'C':4, 'H':8, 'O':4, 'N':0, 'S':0}
xfm['HexN']['RE']['1,5'] = {'C':1, 'H':2, 'O':2, 'N':0, 'S':0}
xfm['HexN']['RE']['2,4'] = {'C':4, 'H':9, 'O':3, 'N':1, 'S':0}
xfm['HexN']['RE']['3,5'] = {'C':3, 'H':7, 'O':3, 'N':1, 'S':0}
xfm['HexN']['RE']['0,3'] = {'C':3, 'H':6, 'O':3, 'N':0, 'S':0}
xfm['HexN']['RE']['1,4'] = {'C':3, 'H':6, 'O':3, 'N':0, 'S':0}
xfm['HexN']['RE']['2,5'] = {'C':2, 'H':5, 'O':2, 'N':1, 'S':0}

## Hex
xfm['Hex'] = {}

# initialize dictionaries for non-reducing end and reducing end
xfm['Hex']['NR'] = {}
xfm['Hex']['RE'] = {}

# add fragments
xfm['Hex']['NR']['0,2'] = {'C':4, 'H':8, 'O':4, 'N':0, 'S':0}
xfm['Hex']['NR']['1,3'] = {'C':2, 'H':4, 'O':2, 'N':0, 'S':0}
xfm['Hex']['NR']['1,5'] = {'C':5, 'H':10, 'O':4, 'N':0, 'S':0}
xfm['Hex']['NR']['2,4'] = {'C':2, 'H':4, 'O':2, 'N':0, 'S':0}
xfm['Hex']['NR']['1,4'] = {'C':3, 'H':6, 'O':3, 'N':0, 'S':0}
xfm['Hex']['NR']['2,5'] = {'C':4, 'H':8, 'O':3, 'N':0, 'S':0}
xfm['Hex']['RE']['0,2'] = {'C':2, 'H':4, 'O':2, 'N':0, 'S':0}
xfm['Hex']['RE']['1,3'] = {'C':4, 'H':8, 'O':4, 'N':0, 'S':0}
xfm['Hex']['RE']['1,5'] = {'C':1, 'H':2, 'O':2, 'N':0, 'S':0}
xfm['Hex']['RE']['2,4'] = {'C':4, 'H':8, 'O':4, 'N':0, 'S':0}
xfm['Hex']['RE']['1,4'] = {'C':3, 'H':6, 'O':3, 'N':0, 'S':0}
xfm['Hex']['RE']['2,5'] = {'C':2, 'H':4, 'O':3, 'N':0, 'S':0}



### cross-ring weights ###
## HexA
xwt['HexA'] = {}

# initialize dictionaries for non-reducing end and reducing end
xwt['HexA']['NR'] = {}
xwt['HexA']['RE'] = {}

# add weights
xwt['HexA']['NR']['0,2'] = 4*wt['monoC'] + 6*wt['monoH'] + 5*wt['monoO']
xwt['HexA']['NR']['1,3'] = 2*wt['monoC'] + 4*wt['monoH'] + 2*wt['monoO']
xwt['HexA']['NR']['1,5'] = 5*wt['monoC'] + 8*wt['monoH'] + 5*wt['monoO']
xwt['HexA']['NR']['2,4'] = 2*wt['monoC'] + 4*wt['monoH'] + 2*wt['monoO']
xwt['HexA']['NR']['3,5'] = 3*wt['monoC'] + 4*wt['monoH'] + 3*wt['monoO']
xwt['HexA']['NR']['0,3'] = 3*wt['monoC'] + 4*wt['monoH'] + 4*wt['monoO']
xwt['HexA']['NR']['1,4'] = 3*wt['monoC'] + 6*wt['monoH'] + 3*wt['monoO']
xwt['HexA']['NR']['2,5'] = 4*wt['monoC'] + 6*wt['monoH'] + 4*wt['monoO']
xwt['HexA']['RE']['0,2'] = 2*wt['monoC'] + 4*wt['monoH'] + 2*wt['monoO']
xwt['HexA']['RE']['1,3'] = 4*wt['monoC'] + 6*wt['monoH'] + 5*wt['monoO']
xwt['HexA']['RE']['1,5'] = 1*wt['monoC'] + 2*wt['monoH'] + 2*wt['monoO']
xwt['HexA']['RE']['2,4'] = 4*wt['monoC'] + 6*wt['monoH'] + 5*wt['monoO']
xwt['HexA']['RE']['3,5'] = 3*wt['monoC'] + 6*wt['monoH'] + 4*wt['monoO']
xwt['HexA']['RE']['0,3'] = 3*wt['monoC'] + 6*wt['monoH'] + 3*wt['monoO']
xwt['HexA']['RE']['1,4'] = 3*wt['monoC'] + 4*wt['monoH'] + 4*wt['monoO']
xwt['HexA']['RE']['2,5'] = 2*wt['monoC'] + 4*wt['monoH'] + 3*wt['monoO']

## dHexA
xwt['dHexA'] = {}

# initialize dictionary for non-reducing end and reducing end
xwt['dHexA']['RE']  = {}

# add weights
xwt['dHexA']['RE']['0,2'] = 2*wt['monoC'] + 4*wt['monoH'] + 2*wt['monoO']
xwt['dHexA']['RE']['1,3'] = 4*wt['monoC'] + 4*wt['monoH'] + 4*wt['monoO']
xwt['dHexA']['RE']['1,5'] = wt['monoC'] + 2*wt['monoH'] + 2*wt['monoO']
xwt['dHexA']['RE']['2,4'] = 4*wt['monoC'] + 5*wt['monoH'] + 5*wt['monoO']
xwt['dHexA']['RE']['3,5'] = 3*wt['monoC'] + 6*wt['monoH'] + 4*wt['monoO']
xwt['dHexA']['RE']['0,3'] = 3*wt['monoC'] + 6*wt['monoH'] + 3*wt['monoO']
xwt['dHexA']['RE']['1,4'] = 3*wt['monoC'] + 3*wt['monoH'] + 4*wt['monoO']
xwt['dHexA']['RE']['2,5'] = 2*wt['monoC'] + 4*wt['monoH'] + 3*wt['monoO']

## HexN
xwt['HexN'] = {}

# initialize dictionaries for non-reducing end and reducing end
xwt['HexN']['NR'] = {}
xwt['HexN']['RE']  = {}

# add weights
xwt['HexN']['NR']['0,2'] = 4*wt['monoC'] + 8*wt['monoH'] + 4*wt['monoO']
xwt['HexN']['NR']['1,3'] = 2*wt['monoC'] + 5*wt['monoH'] + wt['monoO'] + wt['monoN']
xwt['HexN']['NR']['1,5'] = 5*wt['monoC'] + 11*wt['monoH'] + 3*wt['monoO'] + wt['monoN']
xwt['HexN']['NR']['2,4'] = 2*wt['monoC'] + 4*wt['monoH'] + 2*wt['monoO']
xwt['HexN']['NR']['3,5'] = 3*wt['monoC'] + 6*wt['monoH'] + 2*wt['monoO']
xwt['HexN']['NR']['0,3'] = 3*wt['monoC'] + 7*wt['monoH'] + 2*wt['monoO'] + wt['monoN']
xwt['HexN']['NR']['1,4'] = 3*wt['monoC'] + 7*wt['monoH'] + 2*wt['monoO'] + wt['monoN']
xwt['HexN']['NR']['2,5'] = 4*wt['monoC'] + 5*wt['monoH'] + 3*wt['monoO']
xwt['HexN']['RE']['0,2'] = 2*wt['monoC'] + 5*wt['monoH'] + 1*wt['monoO'] + wt['monoN']
xwt['HexN']['RE']['1,3'] = 4*wt['monoC'] + 8*wt['monoH'] + 4*wt['monoO']
xwt['HexN']['RE']['1,5'] = 1*wt['monoC'] + 2*wt['monoH'] + 2*wt['monoO']
xwt['HexN']['RE']['2,4'] = 4*wt['monoC'] + 9*wt['monoH'] + 3*wt['monoO'] + wt['monoN']
xwt['HexN']['RE']['3,5'] = 3*wt['monoC'] + 7*wt['monoH'] + 3*wt['monoO'] + wt['monoN']
xwt['HexN']['RE']['0,3'] = 3*wt['monoC'] + 6*wt['monoH'] + 3*wt['monoO']
xwt['HexN']['RE']['1,4'] = 3*wt['monoC'] + 6*wt['monoH'] + 3*wt['monoO']
xwt['HexN']['RE']['2,5'] = 2*wt['monoC'] + 5*wt['monoH'] + 2*wt['monoO'] + wt['monoN']

## Hex
xwt['Hex'] = {}

# initialize dictionaries for non-reducing end and reducing end
xwt['Hex']['NR'] = {}
xwt['Hex']['RE'] = {}

# add fragments
xwt['Hex']['NR']['0,2'] = 4*wt['monoC'] + 8*wt['monoH'] + 4*wt['monoO']
xwt['Hex']['NR']['1,3'] = 2*wt['monoC'] + 4*wt['monoH'] + 2*wt['monoO']
xwt['Hex']['NR']['1,5'] = 5*wt['monoC'] + 10*wt['monoH'] + 4*wt['monoO']
xwt['Hex']['NR']['2,4'] = 2*wt['monoC'] + 4*wt['monoH'] + 2*wt['monoO']
xwt['Hex']['NR']['1,4'] = 3*wt['monoC'] + 6*wt['monoH'] + 3*wt['monoO']
xwt['Hex']['NR']['2,5'] = 4*wt['monoC'] + 8*wt['monoH'] + 3*wt['monoO']
xwt['Hex']['RE']['0,2'] = 2*wt['monoC'] + 4*wt['monoH'] + 2*wt['monoO']
xwt['Hex']['RE']['1,3'] = 4*wt['monoC'] + 8*wt['monoH'] + 4*wt['monoO']
xwt['Hex']['RE']['1,5'] = wt['monoC'] + 2*wt['monoH'] + 2*wt['monoO']
xwt['Hex']['RE']['2,4'] = 4*wt['monoC'] + 8*wt['monoH'] + 4*wt['monoO']
xwt['Hex']['RE']['1,4'] = 3*wt['monoC'] + 6*wt['monoH'] + 3*wt['monoO']
xwt['Hex']['RE']['2,5'] = 2*wt['monoC'] + 4*wt['monoH'] + 3*wt['monoO']



### cross-ring sulfation/acetylation/adduct possibilities, and COOH places
## HS/Heparin
xmod['HS'] = {}

# HexA
xmod['HS']['HexA'] = {}

# initialize dictionaries for non-reducing end and reducing end
xmod['HS']['HexA']['NR'] = {}
xmod['HS']['HexA']['RE'] = {}

# add modification possibilities
xmod['HS']['HexA']['NR']['0,2'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['HS']['HexA']['NR']['1,5'] = {'SO3':1, 'Ac':0, 'COOH':1}
xmod['HS']['HexA']['NR']['2,4'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['HS']['HexA']['NR']['3,5'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['HS']['HexA']['NR']['0,3'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['HS']['HexA']['NR']['1,4'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['HexA']['NR']['2,5'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['HS']['HexA']['RE']['0,2'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['HexA']['RE']['1,5'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['HS']['HexA']['RE']['2,4'] = {'SO3':1, 'Ac':0, 'COOH':1}
xmod['HS']['HexA']['RE']['3,5'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['HexA']['RE']['0,3'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['HexA']['RE']['1,4'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['HS']['HexA']['RE']['2,5'] = {'SO3':1, 'Ac':0, 'COOH':0}

# dHexA
xmod['HS']['dHexA'] = {}

# initialize dictionary for non-reducing end and reducing end
xmod['HS']['dHexA']['RE']  = {}

# add weights
xmod['HS']['dHexA']['RE']['0,2'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['dHexA']['RE']['1,5'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['HS']['dHexA']['RE']['2,4'] = {'SO3':1, 'Ac':0, 'COOH':1}
xmod['HS']['dHexA']['RE']['3,5'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['dHexA']['RE']['0,3'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['dHexA']['RE']['1,4'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['HS']['dHexA']['RE']['2,5'] = {'SO3':1, 'Ac':0, 'COOH':0}

# HexN
xmod['HS']['HexN'] = {}

# initialize dictionaries for non-reducing end and reducing end
xmod['HS']['HexN']['NR'] = {}
xmod['HS']['HexN']['RE'] = {}

# add modification possibilities
xmod['HS']['HexN']['NR']['0,2'] = {'SO3':2, 'Ac':0, 'COOH':0}
xmod['HS']['HexN']['NR']['1,5'] = {'SO3':3, 'Ac':1, 'COOH':0}
xmod['HS']['HexN']['NR']['2,4'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['HexN']['NR']['3,5'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['HexN']['NR']['0,3'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['HexN']['NR']['1,4'] = {'SO3':2, 'Ac':1, 'COOH':0}
xmod['HS']['HexN']['NR']['2,5'] = {'SO3':2, 'Ac':0, 'COOH':0}
xmod['HS']['HexN']['RE']['0,2'] = {'SO3':1, 'Ac':1, 'COOH':0}
xmod['HS']['HexN']['RE']['1,5'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['HS']['HexN']['RE']['2,4'] = {'SO3':2, 'Ac':1, 'COOH':0}
xmod['HS']['HexN']['RE']['3,5'] = {'SO3':2, 'Ac':1, 'COOH':0}
xmod['HS']['HexN']['RE']['0,3'] = {'SO3':2, 'Ac':1, 'COOH':0}
xmod['HS']['HexN']['RE']['1,4'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['HS']['HexN']['RE']['2,5'] = {'SO3':1, 'Ac':1, 'COOH':0}


## CS/DS
xmod['CS'] = {}

# HexA
xmod['CS']['HexA'] = {}

# initialize dictionaries for non-reducing end and reducing end
xmod['CS']['HexA']['NR'] = {}
xmod['CS']['HexA']['RE'] = {}

# add modification possibilities
xmod['CS']['HexA']['NR']['0,2'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['CS']['HexA']['NR']['1,5'] = {'SO3':1, 'Ac':0, 'COOH':1}
xmod['CS']['HexA']['NR']['2,4'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['CS']['HexA']['NR']['3,5'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['CS']['HexA']['NR']['0,3'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['CS']['HexA']['NR']['1,4'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['CS']['HexA']['NR']['2,5'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['CS']['HexA']['RE']['0,2'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['CS']['HexA']['RE']['1,5'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['CS']['HexA']['RE']['2,4'] = {'SO3':1, 'Ac':0, 'COOH':1}
xmod['CS']['HexA']['RE']['3,5'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['CS']['HexA']['RE']['0,3'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['CS']['HexA']['RE']['1,4'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['CS']['HexA']['RE']['2,5'] = {'SO3':1, 'Ac':0, 'COOH':0}

# dHexA
xmod['CS']['dHexA'] = {}

# initialize dictionary for non-reducing end and reducing end
xmod['CS']['dHexA']['RE']  = {}

# add weights
xmod['CS']['dHexA']['RE']['0,2'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['CS']['dHexA']['RE']['1,5'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['CS']['dHexA']['RE']['2,4'] = {'SO3':1, 'Ac':0, 'COOH':1}
xmod['CS']['dHexA']['RE']['3,5'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['CS']['dHexA']['RE']['0,3'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['CS']['dHexA']['RE']['1,4'] = {'SO3':0, 'Ac':0, 'COOH':1}
xmod['CS']['dHexA']['RE']['2,5'] = {'SO3':1, 'Ac':0, 'COOH':0}

# HexN
xmod['CS']['HexN'] = {}

# initialize dictionaries for non-reducing end and reducing end
xmod['CS']['HexN']['NR'] = {}
xmod['CS']['HexN']['RE'] = {}

# add modification possibilities
xmod['CS']['HexN']['NR']['0,2'] = {'SO3':2, 'Ac':0, 'COOH':0}
xmod['CS']['HexN']['NR']['1,3'] = {'SO3':2, 'Ac':0, 'COOH':0}
xmod['CS']['HexN']['NR']['1,5'] = {'SO3':2, 'Ac':1, 'COOH':0}
xmod['CS']['HexN']['NR']['2,4'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['CS']['HexN']['NR']['1,4'] = {'SO3':1, 'Ac':1, 'COOH':0}
xmod['CS']['HexN']['NR']['2,5'] = {'SO3':2, 'Ac':0, 'COOH':0}
xmod['CS']['HexN']['RE']['0,2'] = {'SO3':0, 'Ac':1, 'COOH':0}
xmod['CS']['HexN']['RE']['1,3'] = {'SO3':0, 'Ac':1, 'COOH':0}
xmod['CS']['HexN']['RE']['1,5'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['CS']['HexN']['RE']['2,4'] = {'SO3':1, 'Ac':1, 'COOH':0}
xmod['CS']['HexN']['RE']['1,4'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['CS']['HexN']['RE']['2,5'] = {'SO3':0, 'Ac':1, 'COOH':0}


## ks
xmod['KS'] = {}

# Hex
xmod['KS']['Hex'] = {}

# initialize dictionaries for non-reducing end and reducing end
xmod['KS']['Hex']['NR'] = {}
xmod['KS']['Hex']['RE'] = {}

# add modification possibilities
xmod['KS']['Hex']['NR']['0,2'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['NR']['1,3'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['NR']['1,5'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['NR']['2,4'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['NR']['1,4'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['NR']['2,5'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['RE']['0,2'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['RE']['1,3'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['RE']['1,5'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['RE']['2,4'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['RE']['1,4'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['Hex']['RE']['2,5'] = {'SO3':0, 'Ac':0, 'COOH':0}

# HexN
xmod['KS']['HexN'] = {}

# initialize dictionaries for non-reducing end and reducing end
xmod['KS']['HexN']['NR'] = {}
xmod['KS']['HexN']['RE'] = {}

# add modification possibilities
xmod['KS']['HexN']['NR']['0,2'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['HexN']['NR']['1,5'] = {'SO3':1, 'Ac':1, 'COOH':0}
xmod['KS']['HexN']['NR']['2,4'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['KS']['HexN']['NR']['3,5'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['HexN']['NR']['0,3'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['HexN']['NR']['1,4'] = {'SO3':0, 'Ac':1, 'COOH':0}
xmod['KS']['HexN']['NR']['2,5'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['HexN']['RE']['0,2'] = {'SO3':0, 'Ac':1, 'COOH':0}
xmod['KS']['HexN']['RE']['1,5'] = {'SO3':0, 'Ac':0, 'COOH':0}
xmod['KS']['HexN']['RE']['2,4'] = {'SO3':1, 'Ac':1, 'COOH':0}
xmod['KS']['HexN']['RE']['3,5'] = {'SO3':0, 'Ac':1, 'COOH':0}
xmod['KS']['HexN']['NR']['0,3'] = {'SO3':0, 'Ac':1, 'COOH':0}
xmod['KS']['HexN']['NR']['1,4'] = {'SO3':1, 'Ac':0, 'COOH':0}
xmod['KS']['HexN']['NR']['2,5'] = {'SO3':0, 'Ac':1, 'COOH':0}

import time # for timing functions
start_time = time.time()

import re # for converting chemical formulae to dictionaries and vice versa
import sys # for getting user inputs
import argparse # for getting user inputs
import sqlite3 as sq # for accessing the database
import pymzml # for handling MS data
import numpy as np # for handling numerical operations
import brainpy as bp # for generating theoretical isotopic distribution

# get individual classes from brainpy
iv    = bp.isotopic_variants
pt    = bp.mass_dict.nist_mass
elems = pt.keys()

debug = False # variable for debugging

# get local package
#from gagfinder.species import *

print "Done!"

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

def binsearch(array, x, hint=None):
    n = len(array)
    lo = 0
    hi = n

    while hi != lo:
        mid = (hi + lo) / 2
        y = array[mid]
        err = y - x
        if hi - lo == 1:
            return mid
        elif err > 0:
            hi = mid
        else:
            lo = mid
    return 0

class Window(object):

    def __init__(self, mz_array, intensity_array, start_index=None, end_index=None,
                 center_mz=None):
        self.mz_array = mz_array
        self.intensity_array = intensity_array
        self.start_index = start_index
        self.end_index = end_index
        self.center_mz = center_mz

    def deduct_intensity(self, value):
        self.intensity_array -= value
        self.intensity_array.clip(min=0, out=self.intensity_array)

    @property
    def mean_intensity(self):
        return self.intensity_array.mean()

    def truncated_mean(self, threshold=0.95):
        try:
            hist_count, hist_values = np.histogram(self.intensity_array)
        except ValueError:
            return 1e-6
        mask = (hist_count.max() * (1 - threshold)) < hist_count
        mean = (hist_values[1:][mask] * hist_count[mask]
                ).sum() / hist_count[mask].sum()
        return mean

    def __repr__(self):
        return "Window(start_index=%d, end_index=%d, mean_intensity=%f, center_mz=%f)" % (
            self.start_index, self.end_index, self.mean_intensity, self.center_mz)

def windowed_spectrum(mz_array, intensity_array, window_size=1.):
    mz_min = mz_array.min()
    mz_max = mz_array.max()

    step_size = window_size / 2.
    center_mz = mz_min + step_size
    center_i = binsearch(mz_array, center_mz, 1)

    windows = []

    niter = 0

    while center_mz < mz_max:
        lo_mz = center_mz - step_size
        hi_mz = center_mz + step_size

        lo_i = binsearch(mz_array, lo_mz, center_i)
        hi_i = binsearch(mz_array, hi_mz, center_i)
        mz_window = mz_array[lo_i:hi_i + 1]
        if lo_mz < mz_window.mean() < hi_mz:
            win = Window(mz_window,
                         intensity_array[lo_i:hi_i + 1], lo_i, hi_i)
        else:
            win = Window(np.array([center_mz]), np.array([0.0]), lo_i, hi_i)
        win.center_mz = center_mz
        windows.append(win)

        center_mz = center_mz + window_size
        center_i = binsearch(mz_array, center_mz, center_i)

        niter += 1

    return windows

class NoiseRegion(object):

    def __init__(self, windows, width=10):
        self.windows = windows
        self.width = width
        self.start_index = windows[0].start_index
        self.end_index = windows[-1].end_index

    def __repr__(self):
        return "NoiseRegion(start_index=%d, end_index=%d)" % (
            self.start_index, self.end_index)

    def noise_window(self):
        return min(self.windows, key=lambda x: x.mean_intensity)

    def arrays(self):
        mz_array = np.concatenate([w.mz_array for w in self.windows])
        intensity_array = np.concatenate(
            [w.intensity_array for w in self.windows])
        return mz_array, intensity_array

    def denoise(self, scale=5, maxiter=10):
        if scale == 0:
            return 0
        noise_mean = self.noise_window().truncated_mean() * scale
        for window in self.windows:
            window.deduct_intensity(noise_mean)
        last_mean = noise_mean
        noise_mean = self.noise_window().truncated_mean()
        niter = 1
        while abs(last_mean - noise_mean) > 1e-3 and niter < maxiter:
            niter += 1
            noise_mean = self.noise_window().truncated_mean() * scale
            for window in self.windows:
                window.deduct_intensity(noise_mean)
            last_mean = noise_mean

        return last_mean - noise_mean

def group_windows(windows, width=10):
    step = int(width / 2)
    regions = []
    i = step
    while i < len(windows):
        lo = i - step
        hi = i + step
        reg = NoiseRegion(windows[lo:hi])
        regions.append(reg)
        i += step * 2
    return regions

class FTICRScan(object):

    def __init__(self, mz_array, intensity_array):
        self.mz_array = mz_array.copy()
        self.intensity_array = intensity_array.copy()

        self.windows = None
        self.regions = None

    def __iter__(self):
        yield self.mz_array
        yield self.intensity_array

    def denoise(self, window_size=1., region_width=10, scale=5):
        self.windows = windowed_spectrum(
            self.mz_array, self.intensity_array, window_size=window_size)
        self.regions = group_windows(self.windows, region_width)
        mz = []
        intensity = []
        for region in self.regions:
            region.denoise(scale)
            mz_, intensity_ = region.arrays()
            mz.extend(mz_)
            intensity.extend(intensity_)
        return FTICRScan(np.array(mz), np.array(intensity))

try:
	_FTICRScan = FTICRScan
	has_c = True
	from ms_peak_picker._c.fticr_denoising import FTICRScan
except ImportError:
	has_c = False

def denoise(mz_array, intensity_array, window_size=1., region_width=10, scale=5):
    scan = FTICRScan(mz_array, intensity_array)
    denoised = scan.denoise(window_size, region_width, scale)
    return list(denoised)

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

# function for summing the scans
def sum_scans(data, spec, noise_removed):
	counter = 0.0 # count the number of scans in the data file
	
	for t in data:
		counter += 1
		
		if t['ms level'] == 2:
			spec += t/sum(t.i)
	
	spec = spec/counter
	
	# noise not already removed by the user, we need to remove, using MasSPIKE
	if not noise_removed:
		backup  = spec.copy() # create a copy of the current metadata
		cleaned = denoise(np.array(spec.mz), 1.e6*np.array(spec.i), 1.) # remove the baseline, after boosting signal
		spec.clear() # remove all metadata and state
		cleaned[1] = cleaned[1]/1.e6 # take signal back to where it was
		
		# set the signal and peak information
		spec.mz, spec.i = map(list, cleaned)
		spec._peaks = list(zip(*map(list, cleaned)))
		
		spec.update(backup) # copy back the metadata
	
	return spec # return

# function for getting information about the precursor
def get_precursor(charge, mz, cursor, deriv_wt, weights, gag_class, adct=None, n_adct=None):
	mass = (mz*abs(charge)) - (charge*weights['H'])
	
	# get proper precursor mass to test
	if not adct:
		n_adct    = 0
		test_mass = mass - deriv_wt
	else:
		test_mass = mass - deriv_wt - (n_adct*weights[adct]) + (n_adct*weights['H'])
	
	# get precursor info from database
	cursor.execute('''SELECT   cpm.id, f.value, p.value, f.monoMass
   	       	 		  FROM     Precursors p, ClassPrecursorMap cpm, Formulae f
	 	       	      WHERE    p.id = cpm.pId 
    	       		  AND      f.id = p.fmId
        	     	  AND      cpm.cId = ?
   	        	 	  ORDER BY ABS(f.monoMass - ?) ASC
	       	     	  LIMIT 1;''', (gag_class, test_mass))
	row = cursor.fetchone()
	
	# return
	return row

# function for getting the reducing and non-reducing end info
def get_ends(pd, gag_class):
	n = pd['D'] + pd['U'] + pd['X'] + pd['N'] # length of GAG
	
	# check if dHexA exists (has to be NR end)
	if pd['D'] > 0:
		nonred = 'D'
		
		if pd['U'] == pd['N']: # we know that the reducing end is HexA because (HexA+dHexA > HexN)
			redend = 'U'
		else: # we know that the reducing end is HexN because (HexA+dHexA == HexN)
			redend = 'N'
	else:
		if gag_class == 4: # KS
			if pd['N'] > pd['X']:
				nonred = 'N'
				redend = 'N'
			elif pd['N'] < pd['X']:
				nonred = 'X'
				redend = 'X'
			else:
				nonred = '?'
				redend = '?'
		else: # HS or CS
			if pd['N'] > pd['U']: # we know that both ends are HexN
				nonred = 'N'
				redend = 'N'
			elif pd['N'] < pd['U']: # we know that both ends are HexA
				nonred = 'U'
				redend = 'U'
			else: # we cannot know which end is which just yet
				nonred = '?'
				redend = '?'
	
	# return
	return [nonred, redend, n]

# function for retrieving all fragment ions
def get_frags(pf, pd, sl, reag, d_dict, d_wt, chrg, n_mono, cursor, cpid, crossmods, redend, adct=None, n_adct=None):
	# variables to store fragment info
	IDs   = {}
	forms = {}
	comps = []
	
	## first look through precursor-based fragments
	ffDict = fmla2dict(pf, 'formula')
	
	# add metal adduct info
	if adct:
		ffDict[adct] = n_adct
		ffDict['H'] -= n_adct
	
	# ready to test fragments
	for i in range(2): # determine how much water to lose
		for j in range(2): # determine whether reducing end fragment or not
			for k in range(3): # cycle through potential hydrogen losses
				for l in range(pd['U'] + pd['D'] + 1): # cycle through potential CO2 losses
					for m in range(min(sl, ffDict['S'])+1): # cycle through potential SO3 losses
						for a in range(2): # cycle through reagent possibilities
							thDict = dict(ffDict) # copy directory of fragment
							
							# add reagent
							for q in reag:
								thDict[q] += a*reag[q]
							
							# H2O loss
							thDict['H'] -= 2*i
							thDict['O'] -= i
							
							# reducing end
							for q in d_dict:
								thDict[q] += j*d_dict[q]
							
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
							if thFmla not in forms:
								for z in range((chrg+1), 0):
									dist             = TheoreticalIsotopicPattern(iv(thDict, charge=z)).truncate_after(0.95)
									IDs[(thFmla, z)] = dist
								
								forms[thFmla] = [thComp]
							else:
								forms[thFmla].append(thComp)
	
	## get all child fragments
	cursor.execute('''SELECT   fr.value, fm.value
	   	    	      FROM     ChildFragments cf, Formulae fm, Fragments fr, Precursors p, ClassPrecursorMap cp
	       	     	  WHERE    cf.frId = fr.id
	       	     	  AND      cf.cpId = cp.id
	       	     	  AND      cp.pId = p.id
	       	     	  AND      fr.fmId = fm.id
	       	     	  AND      cf.cpId = ?;''', (cpid,))
	
	# go through the rows
	for row in cursor.fetchall():
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
		if adct:
			atleast = n_adct - (pd['S'] + pd['U'] + pd['D']) + fcDict['S'] + fcDict['U'] + fcDict['D']
			atmost  = fcDict['S'] + fcDict['U'] + fcDict['D']
		else:
			atleast = 0
			atmost  = 0
			n_adct  = 0
		
		# determine if cross-ring cleavage
		if '+' in fComp:
			xr = 1
			
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
			
			# for adducted metals
			if adct:
				atleast += xmod[gClass][x_ms][x_end][x_clv]['COOH']
				atmost  += xmod[gClass][x_ms][x_end][x_clv]['COOH']
		else:
			xr = 0
		
		# determine if reducing end info is needed
		if d_wt == 0:
			reFrag = [0]
		elif xr: # cross-ring fragment
			if 'NR' in fComp and 'NRE' not in fComp: # non-reducing end fragment
				reFrag = [0]
			else: # reducing end fragment
				reFrag = [1]
		else: # glycosidic fragment
			if fcDict['D'] == 1: # presence of dHexA means non-reducing end fragment
				reFrag = [0]
			else:
				if redend == '?': # don't know reducing end monosaccharide
					reFrag = [0,1] # could be anything
				else: # do know reducing end monosaccharide
					if n_frag % 2 == 0: # equal number, can be either RE or not
						reFrag = [0,1]
					else: # unequal number of Hex/HexA and HexN
						sm = n_frag/2     # smaller number
						bg = (n_frag/2)+1 # bigger number
						
						if fcDict[redend] == bg: # if uneven number, RE must be larger number to possibly be a RE fragment
							if n_frag < n_mono - 1: # doesn't HAVE to be on the end
								reFrag = [0,1]
							else: # HAS to be on the end
								reFrag = [1]
						else: # if RE is smaller number, it cannot be a RE fragment
							reFrag = [0]
		
		## ready to test fragments
		for i in range(2-xr): # determine how much water to lose
			for j in reFrag: # determine whether reducing end fragment or not
				for k in range(-2, 1): # cycle through potential hydrogen losses/gains
					for m in range(max(0, atleast), (min(n_adct, atmost)+1)): # cycle through potential metal adduction
						for so3 in range(min(sl, ffDict['S'])+1): # cycle through potential SO3 losses
							thffDict = dict(ffDict) # copy directory of fragment formula
							thfcDict = dict(fcDict) # copy directory of fragment composition
							
							# H2O loss
							thffDict['H'] -= 2*i
							thffDict['O'] -= i
							
							# reducing end
							if j == 1:
								for q in d_dict:
									thffDict[q] += j*d_dict[q]
							
							# hydrogen loss/gain
							thffDict['H'] += k
							
							# metal adduction
							if adct is not None:
								thffDict[adct] = m
								thffDict['H'] -= m
							
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
							if k < 0:
								if k == -1:
									thComp += '-H'
								elif k == -2:
									thComp += '-2H'
							
							# add metal adduction info
							if m == 1:
								thComp += '+' + adct
							elif m > 1:
								thComp += '+' + str(m) + adct
							
							# check to see if this one has already been searched
							if thComp not in comps:
								comps.append(thComp)
							
								# check if this formula already exists
								if thffFmla not in forms:
									for z in range((chrg+1), 0):
										dist               = TheoreticalIsotopicPattern(iv(thffDict, charge=z)).truncate_after(0.95)
										IDs[(thffFmla, z)] = dist
									
									forms[thffFmla] = [thComp]
								else:
									forms[thffFmla].append(thComp)
	
	# return
	return [IDs, forms]

# function for scoring fragments
def score_frags(IDs, spec, error):
	# dictionaries for storing info
	found = {}
	m_int = {}
	m_mz  = {}
	errs  = {}
	
	# cycle through all fragments
	for j in IDs:
		t = []
		e = []
		
		kill = False # variable to determine whether to kill or not
		
		peak = IDs[j][0] # get first peak to determine precision
		exp  = ((error / 1e6) * peak.mz) + peak.mz # get expected peak position based on precursor error
		
		if exp < spec.mz[0] or exp > spec.mz[len(spec.mz)-1]: # out of the m/z range
			continue
		else: # first peak is in m/z range
			near = spec.hasPeak(exp)
			
			val  = 0
			for k in near:
				if k[1] > val:
					loc = k[0]
					val = k[1]
			
			if val == 0:
				kill = True
			else:
				diff = exp - loc
				ppm  = 1e6 * ((loc - exp)/exp)
			
			t.append(peak.intensity)
			e.append(val)
		
		if kill:
			continue
		
		# get second peak to make certain the hit is "real"
		peak = IDs[j][1]
		if peak.mz < spec.mz[0] or peak.mz > spec.mz[len(spec.mz)-1]: # out of the m/z range
			continue
		else: # second peak is in m/z range
			near = spec.hasPeak(peak.mz)
			
			val  = 0
			prec = 1000000
			for k in near:
				if abs(peak.mz - k[0] - diff) < prec:
					val  = k[1]
					prec = abs(peak.mz - k[0] - diff)
			
			if val == 0:
				kill = True
			
			t.append(peak.intensity)
			e.append(val)
		
		if kill:
			continue
		
		# cycle through remaining peaks
		for peak in IDs[j][2:len(IDs[j])]:
			exp = ((error / 1e6) * peak.mz) + peak.mz # get expected peak position based on precursor error
			
			if exp < spec.mz[0] or exp > spec.mz[len(spec.mz)-1]: # out of the m/z range
				kill = True
				break
				
			near = spec.hasPeak(exp)
			
			val  = 0
			prec = 1000000
			for k in near:
				if abs(exp - k[0] - diff) < prec:
					val  = k[1]
					prec = abs(exp - k[0] - diff)
			
			if val == 0:
				val = 1e-100
			
			e.append(val)
			t.append(peak.intensity)
		
		if not kill:
			EI = np.array(e)/sum(e)
			TI = np.array(t)
			
			g  = 2 * np.sum(EI * np.log(EI/TI))
			
			found[j] = g
			m_int[j] = e[0]
			m_mz [j] = loc
			errs[j]  = ppm
	
	return [found, m_int, m_mz, errs]

# function for running the guts of GAGfinder
def find_gags(mzml_path, gag_class, re_form, N, P, metal, metal_ct, reagent, mz, chg, so3loss, precision, noise_gone, db_path='../lib/GAGfragDB.db'):
	# get values ready for reducing end derivatization and reagent
	df    = {'C':0, 'H':0, 'O':0, 'N':0, 'S':0}
	rf    = {'C':0, 'H':0, 'O':0, 'N':0, 'S':0}
	atoms = ['C','H','O','N','S']
	
	# pick a proper class number
	if gag_class == 'HS':
		cNum = 3
	elif gag_class == 'CS':
		cNum = 1
	else:
		cNum = 4
	
	# parse the reducing end derivatization formula
	if re_form:
		parts = re.findall(r'([A-Z][a-z]*)(\d*)', re_form.upper()) # split formula by symbol
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
	if reagent:
		parts = re.findall(r'([A-Z][a-z]*)(\d*)', reagent.upper()) # split formula by symbol
		for q in parts:
			if q[0] not in atoms: # invalid symbol entered
				print "Invalid chemical formula entered. Please enter only CHONS. Try 'python gagfinder.py --help'"
				sys.exit()
			else:
				if q[1] == '': # only one of this atom
					rf[q[0]] += 1
				else:
					rf[q[0]] += int(q[1])
	
	# weights
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
	
	# get derivatization weight
	dw = 0
	for q in df:
		dw += df[q] * wt[q]
	
	# get reagent weight
	rw = 0
	for q in rf:
		rw += rf[q] * wt[q]
	
	############################################
	# Step 2: Load mzml file and connect to DB #
	############################################
	
	print "Loading mzml file and connecting to GAGfragDB...",
	
	# from user, under construction
	d = pymzml.run.Reader(mzml_path, MSn_Precision=5e-6) # get mzML file into object
	s = pymzml.spec.Spectrum(measuredPrecision=precision) # initialize new Spectrum object
	
	# connect to GAGfragDB
	conn = sq.connect(db_path)
	c    = conn.cursor()
	
	print "Done!"
	
	######################
	# Step 3a: Sum scans #
	######################
	
	print "Summing scans...",
	
	s = sum_scans(d, s, noise_gone) # sum scans
	
	print "Done!"

	########################################
	# Step 3b1: Find precursor composition #
	########################################
	
	print "Finding precursor composition...",
	
	id, pFmla, pComp, tMass = get_precursor(chg, mz, c, dw, wt, cNum, metal, metal_ct) # get info about precursor
	
	print "Done!"

	############################
	# Step 3b2: Find m/z error #
	############################
	
	print "Calculating m/z error...",
	
	tmz = ((tMass + dw) + (chg * 1.00782504)) / abs(chg)
	err = 1e6 * ((mz - tmz)/tmz)
	
	print "Done!"
	
	##########################################################
	# Step 3c: Get reducing end/non-reducing end information #
	##########################################################
	
	print "Determining reducing end/non-reducing end information...",
	
	# convert composition into a dictionary
	pDict         = fmla2dict(pComp, 'composition')
	NR, RE, n_pre = get_ends(pDict, cNum)
	
	print "Done!"
	
	#######################################################################
	# Step 3d: Get potential fragment ions based on precursor composition #
	#######################################################################
	
	print "Retrieving all potential fragment ions for precursor with composition " + pComp + "...",
	
	# get all isotopic distributions
	all_IDs, all_forms = get_frags(pFmla, pDict, so3loss, rf, df, dw, chg, n_pre, c, id, xmod, RE, metal, metal_ct)
	
	print "Done!"
	
	##############################################
	# Step 3e: Score each potential fragment ion #
	##############################################
	
	print "Scoring all potential fragment ions for precursor with composition " + pComp + "...",
	
	# score all fragments
	found_IDs, mono_int, mono_mz, errors = score_frags(all_IDs, s, err)
	
	print "Tested " + str(len(found_IDs)) + " out of " + str(len(all_IDs)) + " fragments"
	
	# rank the fragments by G-score
	rank = sorted(found_IDs.items(), key=lambda x: x[1], reverse=False)
	
	# get top hits
	if N: # user wants top n hits
		top = rank[:N]
	else: # user wants top percentile hits
		pct = P/100.
		ct  = int(pct * len(found_IDs))
		
		top = rank[:ct]
	
	# return
	return [top, found_IDs, mono_int, mono_mz, errors, all_IDs, all_forms]

# function for writing to file
def write_result_to_file(mzml_path, scores, fids, m_mz, m_int, all_dist, formulae, errs):
	print "Done!"
	print "Printing output to file...",
	
	# write to file
	oFile = mzml_path[:-5] + '.tsv'
	f     = open(oFile, 'w')
	f.write("m/z\tIntensity\tCharge\tFragments\tG-score\tError (ppm)\n")
	
	for q in scores:
		out = str(m_mz[q[0]])+'\t'+str(m_int[q[0]])+'\t'+str(all_dist[q[0]][0].charge)+'\t'
		for ions in formulae[q[0][0]]:
			out += ions + '; '
		
		out = out[:-2] + '\t'
		out += str(fids[q[0]]) + '\t' + str(errs[q[0]]) + '\n'
		f.write(out)
	
	f.close()

# main function
def main():
	################################
	# Step 1: check user arguments #
	################################
	
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
	parser.add_argument('-e', type=float, required=False, help='Precision, in ppm (optional, default 20)')
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
	mPrec   = args.e
	removed = args.x
	
	# check to make sure a proper GAG class was added
	if gClass not in ['HS', 'CS', 'KS']:
		print "You must denote a GAG class, either HS, CS, or KS. Try 'python gagfinder.py --help'"
		sys.exit()
	
	# check to make sure that metal adduct is appropriate
	if adduct:
		if adduct not in ['Na', 'K', 'Li', 'Ca', 'Mg']:
			print "\nInvalid metal adduct entered. Only Na, K, Li, Ca, and Mg are accepted. Try 'python gagfinder.py --help'"
			sys.exit()
		
		if nMetal is None or nMetal < 1:
			print "\nYou must enter a positive integer for the number of adducted metals. Try 'python gagfinder.py --help'"
			sys.exit()
	else:
		if nMetal:
			print "\nYou did not select a metal adduct, only the number. Try 'python gagfinder.py --help'"
			sys.exit()
	
	# check to make sure user didn't input both a top n and a top percentile
	if top_n and top_p:
		print "You must select either a number or a percentile to return, not both. Try 'python gagfinder.py -h'"
		sys.exit()
	
	# check to see if the user wants to consider sulfate loss
	if not s_loss:
		s_loss = 0
	
	# check to see if the user has a specific product ion mass error
	if not mPrec:
		mPrec = 20.
	
	# change error to correct scale
	mPrec *= 1.e-6
	
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
	
	print "Done!"
	
	# run the guts of GAGfinder
	result = find_gags(dFile, gClass, fmla, top_n, top_p, adduct, nMetal, reag, pre_mz, pre_z, s_loss, mPrec, removed)
	
	# get individual stuff into variables
	output   = result[0]
	f_IDs    = result[1]
	monint   = result[2]
	monmz    = result[3]
	mistakes = result[4]
	a_IDs    = result[5]
	a_forms  = result[6]
	
	# write to file
	write_result_to_file(dFile, output, f_IDs, monmz, monint, a_IDs, a_forms, mistakes)
	
	print "Finished!"
	print time.time() - start_time

# run main
if __name__ == '__main__':
	main()