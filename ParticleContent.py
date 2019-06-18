#!/usr/bin/env python
#Filename: Configuration.py


 ###################################################################################
#																					#
#								ParticleContent										#
#																					#
#	Purpose:	Specifies the particle content and all self-energies of the model.	#
#																					#
 ###################################################################################

# Specify the particles living in the 2HDM as allowed input values by N2HDMCalc
particles2HDM = {"H1", "H2", "H3", "G0", "A0", "Gp", "Hp", "Gm", "Hm", "A", "Z0", "Wp", "Wm", "NeuE", "NeuEBar", "NeuM", "NeuMBar", "NeuT", "NeuTBar", "El", "ElBar", "Mu", "MuBar", "Tau", "TauBar", "D", "DBar", "U", "UBar", "S", "SBar", "C", "CBar", "B", "BBar", "T", "TBar"}

# Specify the masses that correspond to the FeynCalc particles
particleMasses = {"Gam": "0", "Z0": "MZ", "Wm": "MW", "Wp": "MW", "H1": "MH1", "H2": "MH2", "H3": "MH3", "A0": "MA0", "G0": "0", "Hm": "MHp", "Hp": "MHp", "Gm": "0", "Gp": "0", "NeuE": "0", "NeuEBar": "0", "NeuM": "0", "NeuMBar": "0", "NeuT": "0", "NeuTBar": "0", "El": "ME", "ElBar": "ME", "Mu": "MM", "MuBar": "MM", "Tau": "ML", "TauBar": "ML", "D": "MD", "DBar": "MD", "U": "MU", "UBar": "MU", "S": "MS", "SBar": "MS", "C": "MC", "CBar": "MC", "B": "MB", "BBar": "MB", "T": "MT", "TBar": "MT"}

# Specify the necessary self-energies of the 2HDM.
feynartsScalars = {"H1", "H2", "H3", "G0", "A0", "Gp", "Hp"}      # NOTE: The order of particles has to be the same as in the FeynArts model file
feynartsVectors = {"A", "Z0", "Wp"}                         # NOTE: The order of particles has to be the same as in the FeynArts model file
feynartsFermions = {"NeuE", "NeuM", "NeuT", "El", "Mu", "Tau", "U", "C", "T", "D", "S", "B"}     # NOTE: The order of particles has to be the same as in the FeynArts model file
scalarSelfEnergies = [[1, 1], [1, 1], [1, 1], [1, 1], [1, 1], [1, 1], [2, 2], [2, 2], [2, 2], [3, 3], [3, 3], [3, 3]]       # Each doublet contains the FeynArts particle numbers of the two particles that form a self-energy
scalarSelfSubEnergies = [[1, 1], [2, 2], [3, 3], [1, 2], [1, 3], [2, 3], [1, 1], [2, 2], [1, 2], [1, 1], [2, 2], [1, 2]]
scalarSelfTypes = [["S", "S"], ["S", "S"], ["S", "S"], ["S", "S"], ["S", "S"], ["S", "S"], ["S", "S"], ["S", "S"], ["S", "S"], ["S", "S"], ["S", "S"], ["S", "S"]]
vectorSelfEnergies = [[1, 1], [2, 2], [3, 3], [1, 2]]       # Each doublet contains the FeynArts particle numbers of the two particles that form a self-energy
vectorSelfTypes = [["V", "V"], ["V", "V"], ["V", "V"], ["V", "V"]]
fermionSelfEnergies = [[[1, 1], [1, 1]], [[1, 2], [1, 2]], [[1, 3], [1, 3]], [[2, 1], [2, 1]], [[2, 2], [2, 2]], [[2, 3], [2, 3]], [[3, 1], [3, 1]], [[3, 2], [3, 2]], [[3, 3], [3, 3]], [[4, 1], [4, 1]], [[4, 2], [4, 2]], [[4, 3], [4, 3]], [[3,1], [3,2]], [[3,1], [3,3]], [[3,2], [3,3]], [[4,1], [4,2]], [[4,1], [4,3]], [[4,2], [4,3]]]       # Each doublet contains the FeynArts particle numbers of the two particles that form a self-energy (for fermions, we assume that three generations exist per fermion type)
fermionSelfTypes = [["F", "F"], ["F", "F"], ["F", "F"], ["F", "F"], ["F", "F"], ["F", "F"], ["F", "F"], ["F", "F"], ["F", "F"], ["F", "F"], ["F", "F"], ["F", "F"], ["F", "F"], ["F", "F"], ["F", "F"], ["F", "F"], ["F", "F"], ["F", "F"]]
scalarTadpoles = [[1, 0], [1, 0], [1, 0]]       # Each doublet contains the FeynArts particle numbers of the one particle that forms the tadpole
scalarSubTadpoles = [[1, 0], [2, 0], [3, 0]]
scalarTadpoleTypes = [["S", ""], ["S", ""], ["S", ""]]
