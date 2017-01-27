# initialize dictionaries for weights and formulae
wt = {}
fm = {}

# average element weights
wt['avgH'] = 1.00794
wt['avgC'] = 12.0107
wt['avgO'] = 15.9994
wt['avgN'] = 14.0067
wt['avgS'] = 32.0650

# monoisotopic element weights
wt['monoH'] = 1.00782504
wt['monoC'] = 12.0
wt['monoO'] = 15.99491461956
wt['monoN'] = 14.0030740048
wt['monoS'] = 31.97207100

# avg compound weights
wt['avgHex']   = 6*wt['avgC'] + 12*wt['avgH'] + 6*wt['avgO']
wt['avgHexA']  = 6*wt['avgC'] + 10*wt['avgH'] + 7*wt['avgO']
wt['avgHexN']  = 6*wt['avgC'] + 13*wt['avgH'] + wt['avgN'] + 5*wt['avgO']
wt['avgdHexA'] = 6*wt['avgC'] + 8*wt['avgH'] + 6*wt['avgO']
wt['avgH2O']   = 2*wt['avgH'] + wt['avgO']
wt['avgSO3']   = wt['avgS'] + 3*wt['avgO']
wt['avgAc']    = 2*wt['avgC'] + 2*wt['avgH'] + wt['avgO']

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



### cross-ring sulfation/acetylation possibilities
## HS/Heparin
xmod['HS'] = {}

# HexA
xmod['HS']['HexA'] = {}

# initialize dictionaries for non-reducing end and reducing end
xmod['HS']['HexA']['NR'] = {}
xmod['HS']['HexA']['RE'] = {}

# add modification possibilities
xmod['HS']['HexA']['NR']['0,2'] = {'SO3':0, 'Ac':0}
xmod['HS']['HexA']['NR']['1,5'] = {'SO3':1, 'Ac':0}
xmod['HS']['HexA']['NR']['2,4'] = {'SO3':0, 'Ac':0}
xmod['HS']['HexA']['NR']['3,5'] = {'SO3':0, 'Ac':0}
xmod['HS']['HexA']['NR']['0,3'] = {'SO3':0, 'Ac':0}
xmod['HS']['HexA']['NR']['1,4'] = {'SO3':1, 'Ac':0}
xmod['HS']['HexA']['NR']['2,5'] = {'SO3':0, 'Ac':0}
xmod['HS']['HexA']['RE']['0,2'] = {'SO3':1, 'Ac':0}
xmod['HS']['HexA']['RE']['1,5'] = {'SO3':0, 'Ac':0}
xmod['HS']['HexA']['RE']['2,4'] = {'SO3':1, 'Ac':0}
xmod['HS']['HexA']['RE']['3,5'] = {'SO3':1, 'Ac':0}
xmod['HS']['HexA']['RE']['0,3'] = {'SO3':1, 'Ac':0}
xmod['HS']['HexA']['RE']['1,4'] = {'SO3':0, 'Ac':0}
xmod['HS']['HexA']['RE']['2,5'] = {'SO3':1, 'Ac':0}

# dHexA
xmod['HS']['dHexA'] = {}

# initialize dictionary for non-reducing end and reducing end
xmod['HS']['dHexA']['RE']  = {}

# add weights
xmod['HS']['dHexA']['RE']['0,2'] = {'SO3':1, 'Ac':0}
xmod['HS']['dHexA']['RE']['1,5'] = {'SO3':0, 'Ac':0}
xmod['HS']['dHexA']['RE']['2,4'] = {'SO3':1, 'Ac':0}
xmod['HS']['dHexA']['RE']['3,5'] = {'SO3':1, 'Ac':0}
xmod['HS']['dHexA']['RE']['0,3'] = {'SO3':1, 'Ac':0}
xmod['HS']['dHexA']['RE']['1,4'] = {'SO3':0, 'Ac':0}
xmod['HS']['dHexA']['RE']['2,5'] = {'SO3':1, 'Ac':0}

# HexN
xmod['HS']['HexN'] = {}

# initialize dictionaries for non-reducing end and reducing end
xmod['HS']['HexN']['NR'] = {}
xmod['HS']['HexN']['RE'] = {}

# add modification possibilities
xmod['HS']['HexN']['NR']['0,2'] = {'SO3':2, 'Ac':0}
xmod['HS']['HexN']['NR']['1,5'] = {'SO3':3, 'Ac':1}
xmod['HS']['HexN']['NR']['2,4'] = {'SO3':1, 'Ac':0}
xmod['HS']['HexN']['NR']['3,5'] = {'SO3':1, 'Ac':0}
xmod['HS']['HexN']['NR']['0,3'] = {'SO3':1, 'Ac':0}
xmod['HS']['HexN']['NR']['1,4'] = {'SO3':2, 'Ac':1}
xmod['HS']['HexN']['NR']['2,5'] = {'SO3':2, 'Ac':0}
xmod['HS']['HexN']['RE']['0,2'] = {'SO3':1, 'Ac':1}
xmod['HS']['HexN']['RE']['1,5'] = {'SO3':0, 'Ac':0}
xmod['HS']['HexN']['RE']['2,4'] = {'SO3':2, 'Ac':1}
xmod['HS']['HexN']['RE']['3,5'] = {'SO3':2, 'Ac':1}
xmod['HS']['HexN']['RE']['0,3'] = {'SO3':2, 'Ac':1}
xmod['HS']['HexN']['RE']['1,4'] = {'SO3':1, 'Ac':0}
xmod['HS']['HexN']['RE']['2,5'] = {'SO3':1, 'Ac':1}


## CS/DS
xmod['CS'] = {}

# HexA
xmod['CS']['HexA'] = {}

# initialize dictionaries for non-reducing end and reducing end
xmod['CS']['HexA']['NR'] = {}
xmod['CS']['HexA']['RE'] = {}

# add modification possibilities
xmod['CS']['HexA']['NR']['0,2'] = {'SO3':0, 'Ac':0}
xmod['CS']['HexA']['NR']['1,5'] = {'SO3':1, 'Ac':0}
xmod['CS']['HexA']['NR']['2,4'] = {'SO3':0, 'Ac':0}
xmod['CS']['HexA']['NR']['3,5'] = {'SO3':0, 'Ac':0}
xmod['CS']['HexA']['NR']['0,3'] = {'SO3':0, 'Ac':0}
xmod['CS']['HexA']['NR']['1,4'] = {'SO3':1, 'Ac':0}
xmod['CS']['HexA']['NR']['2,5'] = {'SO3':0, 'Ac':0}
xmod['CS']['HexA']['RE']['0,2'] = {'SO3':1, 'Ac':0}
xmod['CS']['HexA']['RE']['1,5'] = {'SO3':0, 'Ac':0}
xmod['CS']['HexA']['RE']['2,4'] = {'SO3':1, 'Ac':0}
xmod['CS']['HexA']['RE']['3,5'] = {'SO3':1, 'Ac':0}
xmod['CS']['HexA']['RE']['0,3'] = {'SO3':1, 'Ac':0}
xmod['CS']['HexA']['RE']['1,4'] = {'SO3':0, 'Ac':0}
xmod['CS']['HexA']['RE']['2,5'] = {'SO3':1, 'Ac':0}

# dHexA
xmod['CS']['dHexA'] = {}

# initialize dictionary for non-reducing end and reducing end
xmod['CS']['dHexA']['RE']  = {}

# add weights
xmod['CS']['dHexA']['RE']['0,2'] = {'SO3':1, 'Ac':0}
xmod['CS']['dHexA']['RE']['1,5'] = {'SO3':0, 'Ac':0}
xmod['CS']['dHexA']['RE']['2,4'] = {'SO3':1, 'Ac':0}
xmod['CS']['dHexA']['RE']['3,5'] = {'SO3':1, 'Ac':0}
xmod['CS']['dHexA']['RE']['0,3'] = {'SO3':1, 'Ac':0}
xmod['CS']['dHexA']['RE']['1,4'] = {'SO3':0, 'Ac':0}
xmod['CS']['dHexA']['RE']['2,5'] = {'SO3':1, 'Ac':0}

# HexN
xmod['CS']['HexN'] = {}

# initialize dictionaries for non-reducing end and reducing end
xmod['CS']['HexN']['NR'] = {}
xmod['CS']['HexN']['RE'] = {}

# add modification possibilities
xmod['CS']['HexN']['NR']['0,2'] = {'SO3':2, 'Ac':0}
xmod['CS']['HexN']['NR']['1,3'] = {'SO3':2, 'Ac':0}
xmod['CS']['HexN']['NR']['1,5'] = {'SO3':2, 'Ac':1}
xmod['CS']['HexN']['NR']['2,4'] = {'SO3':1, 'Ac':0}
xmod['CS']['HexN']['NR']['1,4'] = {'SO3':1, 'Ac':1}
xmod['CS']['HexN']['NR']['2,5'] = {'SO3':2, 'Ac':0}
xmod['CS']['HexN']['RE']['0,2'] = {'SO3':0, 'Ac':1}
xmod['CS']['HexN']['RE']['1,3'] = {'SO3':0, 'Ac':1}
xmod['CS']['HexN']['RE']['1,5'] = {'SO3':0, 'Ac':0}
xmod['CS']['HexN']['RE']['2,4'] = {'SO3':1, 'Ac':1}
xmod['CS']['HexN']['RE']['1,4'] = {'SO3':1, 'Ac':0}
xmod['CS']['HexN']['RE']['2,5'] = {'SO3':0, 'Ac':1}


## ks
xmod['KS'] = {}

# Hex
xmod['KS']['Hex'] = {}

# initialize dictionaries for non-reducing end and reducing end
xmod['KS']['Hex']['NR'] = {}
xmod['KS']['Hex']['RE'] = {}

# add modification possibilities
xmod['KS']['Hex']['NR']['0,2'] = {'SO3':1, 'Ac':0}
xmod['KS']['Hex']['NR']['1,3'] = {'SO3':0, 'Ac':0}
xmod['KS']['Hex']['NR']['1,5'] = {'SO3':1, 'Ac':0}
xmod['KS']['Hex']['NR']['2,4'] = {'SO3':0, 'Ac':0}
xmod['KS']['Hex']['NR']['1,4'] = {'SO3':0, 'Ac':0}
xmod['KS']['Hex']['NR']['2,5'] = {'SO3':1, 'Ac':0}
xmod['KS']['Hex']['RE']['0,2'] = {'SO3':0, 'Ac':0}
xmod['KS']['Hex']['RE']['1,3'] = {'SO3':1, 'Ac':0}
xmod['KS']['Hex']['RE']['1,5'] = {'SO3':0, 'Ac':0}
xmod['KS']['Hex']['RE']['2,4'] = {'SO3':1, 'Ac':0}
xmod['KS']['Hex']['RE']['1,4'] = {'SO3':1, 'Ac':0}
xmod['KS']['Hex']['RE']['2,5'] = {'SO3':0, 'Ac':0}

# HexN
xmod['KS']['HexN'] = {}

# initialize dictionaries for non-reducing end and reducing end
xmod['KS']['HexN']['NR'] = {}
xmod['KS']['HexN']['RE'] = {}

# add modification possibilities
xmod['KS']['HexN']['NR']['0,2'] = {'SO3':1, 'Ac':0}
xmod['KS']['HexN']['NR']['1,5'] = {'SO3':1, 'Ac':1}
xmod['KS']['HexN']['NR']['2,4'] = {'SO3':0, 'Ac':0}
xmod['KS']['HexN']['NR']['3,5'] = {'SO3':1, 'Ac':0}
xmod['KS']['HexN']['NR']['0,3'] = {'SO3':1, 'Ac':0}
xmod['KS']['HexN']['NR']['1,4'] = {'SO3':0, 'Ac':1}
xmod['KS']['HexN']['NR']['2,5'] = {'SO3':1, 'Ac':0}
xmod['KS']['HexN']['RE']['0,2'] = {'SO3':0, 'Ac':1}
xmod['KS']['HexN']['RE']['1,5'] = {'SO3':0, 'Ac':0}
xmod['KS']['HexN']['RE']['2,4'] = {'SO3':1, 'Ac':1}
xmod['KS']['HexN']['RE']['3,5'] = {'SO3':0, 'Ac':1}
xmod['KS']['HexN']['NR']['0,3'] = {'SO3':0, 'Ac':1}
xmod['KS']['HexN']['NR']['1,4'] = {'SO3':1, 'Ac':0}
xmod['KS']['HexN']['NR']['2,5'] = {'SO3':0, 'Ac':1}