# per residue sidechain hedra and dihedra definitions, in order of output for
# internal coordinates specification file

pic_data_backbone = {
    
}
pic_data_sidechains = {
    'V': [
        ['CA', 'CB', 'CG1'],
        ['N', 'CA', 'CB', 'CG1', 'chi1'],  # chi1
        ['CA', 'CB', 'CG2'],
        ['N', 'CA', 'CB', 'CG2'],
        # 'chi1': 1
    ],
    'L': [
        ['CA', 'CB', 'CG'],
        ['N', 'CA', 'CB', 'CG', 'chi1'],  # chi1
        ['CB', 'CG', 'CD1'],
        ['CA', 'CB', 'CG', 'CD1', 'chi2'],  # chi2
        ['CB', 'CG', 'CD2'],
        ['CA', 'CB', 'CG', 'CD2'],
        # 'chi1': 1, 'chi2': 3
    ],
    'I': [
        ['CA', 'CB', 'CG1'],
        ['N', 'CA', 'CB', 'CG1', 'chi1'],   # chi1
        ['CB', 'CG1', 'CD1'],
        ['CA', 'CB', 'CG1', 'CD1', 'chi2'],   # chi2
        ['CA', 'CB', 'CG2'],
        ['N', 'CA', 'CB', 'CG2'],
        # 'chi1': 1, 'chi2': 3
    ],
    'M': [
        ['CA', 'CB', 'CG'],
        ['N', 'CA', 'CB', 'CG', 'chi1'],   # chi1
        ['CB', 'CG', 'SD'],
        ['CA', 'CB', 'CG', 'SD', 'chi2'],   # chi2
        ['CG', 'SD', 'CE'],
        ['CB', 'CG', 'SD', 'CE', 'chi3'],   # chi3
        # 'chi1': 2, 'chi2': 4, 'chi3': 6
    ],
    'F': [
        ['CA', 'CB', 'CG'],
        ['N', 'CA', 'CB', 'CG', 'chi1'],   # chi1
        ['CB', 'CG', 'CD1'],
        ['CA', 'CB', 'CG', 'CD1', 'chi2'],   # chi2
        ['CG', 'CD1', 'CE1'],
        ['CB', 'CG', 'CD1', 'CE1'],
        ['CD1', 'CE1', 'CZ'],
        ['CG', 'CD1', 'CE1', 'CZ'],
        ['CB', 'CG', 'CD2'],
        ['CA', 'CB', 'CG', 'CD2'],
        ['CG', 'CD2', 'CE2'],
        ['CB', 'CG', 'CD2', 'CE2'],
        # 'chi1': 1, 'chi2': 3
    ],
    'P': [
        ['CA', 'CB', 'CG'],
        ['N', 'CA', 'CB', 'CG', 'chi1'],   # chi1
        ['CB', 'CG', 'CD'],
        ['CA', 'CB', 'CG', 'CD', 'chi2'],   # chi2
        # 'chi1': 1, 'chi2': 3
    ],
    'S': [
        ['CA', 'CB', 'OG'],
        ['N', 'CA', 'CB', 'OG', 'chi1'],   # chi1
        # 'chi1': 1
    ],
    'T': [
        ['CA', 'CB', 'OG1'],
        ['N', 'CA', 'CB', 'OG1', 'chi1'],   # chi1
        ['CA', 'CB', 'CG2'],
        ['N', 'CA', 'CB', 'CG2'],
        # 'chi1': 1
    ],
    'C': [
        ['CA', 'CB', 'SG'],
        ['N', 'CA', 'CB', 'SG', 'chi1'],   # chi1
        # 'chi1': 1
    ],
    'N': [
        ['CA', 'CB', 'CG'],
        ['N', 'CA', 'CB', 'CG', 'chi1'],   # chi1
        ['CB', 'CG', 'OD1'],
        ['CA', 'CB', 'CG', 'OD1', 'chi2'],   # chi2
        ['CB', 'CG', 'ND2'],
        ['CA', 'CB', 'CG', 'ND2'],
        # 'chi1': 1, 'chi2': 1
    ],
    'Q': [
        ['CA', 'CB', 'CG'],
        ['N', 'CA', 'CB', 'CG', 'chi1'],   # chi1
        ['CB', 'CG', 'CD'],
        ['CA', 'CB', 'CG', 'CD', 'chi2'],   # chi2
        ['CG', 'CD', 'OE1'],
        ['CB', 'CG', 'CD', 'OE1', 'chi3'],   # chi3
        ['CG', 'CD', 'NE2'],
        ['CB', 'CG', 'CD', 'NE2'],
        # 'chi1': 1, 'chi2': 3, 'chi3': 5
    ],
    'Y': [
        ['CA', 'CB', 'CG'],
        ['N', 'CA', 'CB', 'CG', 'chi1'],   # chi1
        ['CB', 'CG', 'CD1'],
        ['CA', 'CB', 'CG', 'CD1', 'chi2'],   # chi2
        ['CG', 'CD1', 'CE1'],
        ['CB', 'CG', 'CD1', 'CE1'],
        ['CD1', 'CE1', 'CZ'],
        ['CG', 'CD1', 'CE1', 'CZ'],
        ['CE1', 'CZ', 'OH'],
        ['CD1', 'CE1', 'CZ', 'OH'],
        ['CB', 'CG', 'CD2'],
        ['CA', 'CB', 'CG', 'CD2'],
        ['CG', 'CD2', 'CE2'],
        ['CB', 'CG', 'CD2', 'CE2'],
        # 'chi1': 1, 'chi2': 3
    ],
    'W': [
        ['CA', 'CB', 'CG'],
        ['N', 'CA', 'CB', 'CG', 'chi1'],   # chi1
        ['CB', 'CG', 'CD1'],
        ['CA', 'CB', 'CG', 'CD1', 'chi2'],   # chi2
        ['CG', 'CD1', 'NE1'],
        ['CB', 'CG', 'CD1', 'NE1'],
        ['CB', 'CG', 'CD2'],
        ['CA', 'CB', 'CG', 'CD2'],
        ['CG', 'CD2', 'CE2'],
        ['CB', 'CG', 'CD2', 'CE2'],
        ['CD2', 'CE2', 'CZ2'],
        ['CG', 'CD2', 'CE2', 'CZ2'],
        ['CE2', 'CZ2', 'CH2'],
        ['CD2', 'CE2', 'CZ2', 'CH2'],
        ['CG', 'CD2', 'CE3'],
        ['CB', 'CG', 'CD2', 'CE3'],
        ['CD2', 'CE3', 'CZ3'],
        ['CG', 'CD2', 'CE3', 'CZ3'],
        # 'chi1': 1, 'chi2': 3
    ],
    'D': [
        ['CA', 'CB', 'CG'],
        ['N', 'CA', 'CB', 'CG', 'chi1'],   # chi1
        ['CB', 'CG', 'OD1'],
        ['CA', 'CB', 'CG', 'OD1', 'chi2'],   # chi2
        ['CB', 'CG', 'OD2'],
        ['CA', 'CB', 'CG', 'OD2'],
        # 'chi1': 1, 'chi2': 3
    ],
    'E': [
        ['CA', 'CB', 'CG'],
        ['N', 'CA', 'CB', 'CG', 'chi1'],   # chi1
        ['CB', 'CG', 'CD'],
        ['CA', 'CB', 'CG', 'CD', 'chi2'],   # chi2
        ['CG', 'CD', 'OE1'],
        ['CB', 'CG', 'CD', 'OE1', 'chi3'],   # chi3
        ['CG', 'CD', 'OE2'],
        ['CB', 'CG', 'CD', 'OE2'],
        # 'chi1': 1, 'chi2': 3, 'chi3': 5
    ],
    'H': [
        ['CA', 'CB', 'CG'],
        ['N', 'CA', 'CB', 'CG', 'chi1'],   # chi1
        ['CB', 'CG', 'ND1'],
        ['CA', 'CB', 'CG', 'ND1', 'chi2'],   # chi2
        ['CG', 'ND1', 'CE1'],
        ['CB', 'CG', 'ND1', 'CE1'],
        ['CB', 'CG', 'CD2'],
        ['CA', 'CB', 'CG', 'CD2'],
        ['CG', 'CD2', 'NE2'],
        ['CB', 'CG', 'CD2', 'NE2'],
        # 'chi1': 1, 'chi2': 3
    ],
    'K': [
        ['CA', 'CB', 'CG'],
        ['N', 'CA', 'CB', 'CG', 'chi1'],   # chi1
        ['CB', 'CG', 'CD'],
        ['CA', 'CB', 'CG', 'CD', 'chi2'],   # chi2
        ['CG', 'CD', 'CE'],
        ['CB', 'CG', 'CD', 'CE', 'chi3'],   # chi3
        ['CD', 'CE', 'NZ'],
        ['CG', 'CD', 'CE', 'NZ', 'chi4'],   # chi4
        # 'chi1': 1, 'chi2': 3, 'chi3': 5, 'chi4': 7
    ],
    'R': [
        ['CA', 'CB', 'CG'],
        ['N', 'CA', 'CB', 'CG', 'chi1'],   # chi1
        ['CB', 'CG', 'CD'],
        ['CA', 'CB', 'CG', 'CD', 'chi2'],   # chi2
        ['CG', 'CD', 'NE'],
        ['CB', 'CG', 'CD', 'NE', 'chi3'],   # chi3
        ['CD', 'NE', 'CZ'],
        ['CG', 'CD', 'NE', 'CZ', 'chi4'],   # chi4
        ['NE', 'CZ', 'NH1'],
        ['CD', 'NE', 'CZ', 'NH1', 'chi5'],   # chi5
        ['NE', 'CZ', 'NH2'],
        ['CD', 'NE', 'CZ', 'NH2'],
        # 'chi1': 1, 'chi2': 3, 'chi3': 5, 'chi4': 7, 'chi5': 9
    ]
}

# http://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/formuleAA/ 
# for naming of individual atoms
