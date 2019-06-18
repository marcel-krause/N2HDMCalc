#!/usr/bin/env python
#Filename: ExportToEW2HDECAY.py


 ###################################################################################
#																					#
#							ExportToEW2HDECAY.py									#
#																					#
#	Purpose:	Mirrors all relevant files of N2HDMCalc to the ewN2HDECAY folder. 	#
#				WARNING: this UPDATES all files in the ewN2HDECAY folder. Changes 	#
#				to the relevant files should only be made in the N2HDMCalc folder!	#
#	Application: This file should be saved in the N2HDMCalc folder. The ewN2HDECAY	#
#				folder HAS TO BE a parallel folder of the N2HDMCalc folder, i.e. 	#
#				both the N2HDMCalc and ewN2HDECAY folders should be in the same 	#
#				directory.															#
#	Author: 	Marcel Krause (marcel.krause@kit.edu)								#
#																					#
 ###################################################################################


#------------------------------#
#		 Import Modules		   #
#------------------------------#
import sys
import os
from shutil import copyfile, rmtree


#-------------------------#
#		 Functions		  #
#-------------------------#

# def testFunc(testArg):
# 	'''
# 		Descr
# 	'''
	

#----------------------------#
#		 Main Program		 #
#----------------------------#

if __name__ == "__main__":
	# Exclude list for file extensions
	excludeExtensions = ['.o', '.exe', '.a', '.mod', '.stackdump', '.gch']

	# List of relevant folders
	rootFolders = ['BuildingBlocks', 'BuildingBlocks' + os.sep + 'SelfEnergies', 'BuildingBlocks' + os.sep + 'Processes', 'Parameters']
	relevantFolders = ['BuildingBlocks' + os.sep + 'ProcessDependentScheme', 'BuildingBlocks' + os.sep + 'SelfEnergies' + os.sep + 'Usual', 'BuildingBlocks' + os.sep + 'SelfEnergies' + os.sep + 'Alternative', 'BuildingBlocks' + os.sep + 'SelfEnergiesDerivatives', 'BuildingBlocks' + os.sep + 'Tadpoles']
	processFolder = 'BuildingBlocks' + os.sep + 'Processes'
	processFiles = ['Counterterm.F90', 'NLOTadWidthRed.F90', 'NLOWidthRed.F90', 'RealCorrections.F90', 'TreeLevelWidthRed.F90', 'processDescription.txt']
	additionalFiles = ['counterterms.F90']
	rootDirewN2HDECAY = '..' + os.sep + 'ewN2HDECAY Development'

	# Converter for HDECAY conventions
	convertList = [["H1 -> B,BBar", "H1 -> bb", True],["H1 -> Tau,TauBar", "H1 -> tautau", True],["H1 -> Mu,MuBar", "H1 -> mu mu", True],["H1 -> S,SBar", "H1 -> ss", True],["H1 -> C,CBar", "H1 -> cc", True],["H1 -> T,TBar", "H1 -> tt", True],["H1 -> Wm,Wp", "H1 -> WW", True],["H1 -> Z0,Z0", "H1 -> ZZ", True],["H1 -> A0,A0", "H1 -> AA", True],["H1 -> A0,Z0", "H1 -> ZA", True],["H1 -> Hp,Wm", "H1 -> H+ W-", True],["H1 -> Hm,Hp", "H1 -> H+ H-", True],["H2 -> B,BBar", "H2 -> bb", True],["H2 -> Tau,TauBar", "H2 -> tau tau", True],["H2 -> Mu,MuBar", "H2 -> mu mu", True],["H2 -> S,SBar", "H2 -> ss", True],["H2 -> C,CBar", "H2 -> cc", True],["H2 -> T,TBar", "H2 -> tt", True],["H2 -> Wm,Wp", "H2 -> WW", True],["H2 -> Z0,Z0", "H2 -> ZZ", True],["H2 -> H1,H1", "H2 -> H1H1", True],["H2 -> A0,A0", "H2 -> AA", True],["H2 -> A0,Z0", "H2 -> ZA", True],["H2 -> Hp,Wm", "H2 -> H+ W-", True],["H2 -> Hm,Hp", "H2 -> H+ H-", True],["H3 -> B,BBar", "H3 -> bb", True],["H3 -> Tau,TauBar", "H3 -> tau tau", True],["H3 -> Mu,MuBar", "H3 -> mu mu", True],["H3 -> S,SBar", "H3 -> ss", True],["H3 -> C,CBar", "H3 -> cc", True],["H3 -> T,TBar", "H3 -> tt", True],["H3 -> Wm,Wp", "H3 -> WW", True],["H3 -> Z0,Z0", "H3 -> ZZ", True],["H3 -> H1,H1", "H3 -> H1H1", True],["H3 -> H1,H2", "H3 -> H1H2", True],["H3 -> H2,H2", "H3 -> H2H2", True],["H3 -> A0,A0", "H3 -> AA", True],["H3 -> A0,Z0", "H3 -> ZA", True],["H3 -> Hp,Wm", "H3 -> H+ W-", True],["H3 -> Hm,Hp", "H3 -> H+ H-", True],["A0 -> B,BBar", "A -> bb", True],["A0 -> Tau,TauBar", "A -> tau,tau", True],["A0 -> Mu,MuBar", "A -> mu,mu", True],["A0 -> S,SBar", "A -> ss", True],["A0 -> C,CBar", "A -> cc", True],["A0 -> T,TBar", "A -> tt", True],["A0 -> H1,Z0", "A -> ZH1", True],["A0 -> H2,Z0", "A -> ZH2", True],["A0 -> H3,Z0", "A -> ZH3", True],["A0 -> Hm,Wp", "A -> H- W+", True],["Hp -> BBar,C", "H+ -> bbar c", True],["Hp -> NeuT,TauBar", "H+ -> taubar nu", True],["Hp -> MuBar,NeuM", "H+ -> mubar nu", True],["Hp -> SBar,U", "H+ -> sbar u", True],["Hp -> C,SBar", "H+ -> sbar c", True],["Hp -> BBar,T", "H+ -> bbar t", True],["Hp -> C,DBar", "H+ -> dbar c", True],["Hp -> BBar,U", "H+ -> bbar u", True],["Hp -> SBar,T", "H+ -> sbar t", True],["Hp -> DBar,T", "H+ -> dbar t", True],["Hp -> H1,Wp", "H+ -> W H1", True],["Hp -> H2,Wp", "H+ -> W H2", True],["Hp -> H3,Wp", "H+ -> W H3", True],["Hp -> A0,Wp", "H+ -> W A", True],["A0 -> D,DBar", "A0 -> D,DBar", False],["A0 -> El,ElBar", "A0 -> El,ElBar", False],["A0 -> U,UBar", "A0 -> U,UBar", False],["H1 -> D,DBar", "H1 -> D,DBar", False],["H1 -> El,ElBar", "H1 -> El,ElBar", False],["H1 -> U,UBar", "H1 -> U,UBar", False],["H2 -> D,DBar", "H2 -> D,DBar", False],["H2 -> El,ElBar", "H2 -> El,ElBar", False],["H2 -> U,UBar", "H2 -> U,UBar", False],["H3 -> D,DBar", "H3 -> D,DBar", False],["H3 -> El,ElBar", "H3 -> El,ElBar", False],["H3 -> U,UBar", "H3 -> U,UBar", False],["Hp -> DBar,U", "Hp -> DBar,U", False],["Hp -> ElBar,NeuE", "Hp -> ElBar,NeuE", False]]

	# Start message 
	print('Start of the copying process.\n')

	# Check if the ewN2HDECAY root dir exists and is in the right location
	print('Checking for ewN2HDECAY folder...')
	if os.path.isdir(rootDirewN2HDECAY):
		print('  ewN2HDECAY folder found.\n')
	else:
		print('  Critical error: ewN2HDECAY folder not found!\n')
		sys.exit()

	# Check if the necessary subfolders exists and if not, create them 
	print('Checking for main subfolders...')
	changedFolder = False
	for folderPath in rootFolders:
		tempPath = rootDirewN2HDECAY + os.sep + folderPath
		if not os.path.isdir(tempPath):
			os.makedirs(tempPath)
			print('  Folder "{}" created.'.format(tempPath))
			changedFolder = True
		# else:
		# 	print('  Folder "{}" already exists.'.format(tempPath))
	for folderPath in relevantFolders:
		tempPath = rootDirewN2HDECAY + os.sep + folderPath
		if not os.path.isdir(tempPath):
			os.makedirs(tempPath)
			print('  Folder "{}" created.'.format(tempPath))
			changedFolder = True
		# else:
		# 	print('  Folder "{}" already exists.'.format(tempPath))
	for folderPath in os.listdir(processFolder):
		tempPath = rootDirewN2HDECAY + os.sep + processFolder + os.sep + folderPath
		if not os.path.isdir(tempPath):
			os.makedirs(tempPath)
			print('  Folder "{}" created.'.format(tempPath))
			changedFolder = True
		# else:
		# 	print('  Folder "{}" already exists.'.format(tempPath))
	if not changedFolder:
		print('  All necessary folders already exist.')

	# Iterate through the relevant folders and update all files that have changed
	print('\nCopying relevant files...')
	for folderPath in relevantFolders:
		dirList = os.listdir(folderPath)
		targetList = []

		# Remove all files which have the file extensions as given in the extension list
		for fileCandidate in dirList:
			for excludeCandidate in excludeExtensions:
				if excludeCandidate in fileCandidate:
					dirList.remove(fileCandidate)
		dirListFinal = [folderPath + os.sep + s for s in dirList]

		# Iterate through the final file list; if the file does not exist, create it. Else, check if the two files differ (compare with latest modification timestamp) and if so, copy the file
		fileCreated = False
		fileUpdated = False
		fileCreatedCounter = 0
		fileUpdatedCounter = 0
		print('  Copying files from folder "{}".'.format(folderPath))
		for fileInN2HDMCalc in dirListFinal:
			fileInewN2HDECAY = rootDirewN2HDECAY + os.sep + fileInN2HDMCalc
			# Check if the file does not exist yet 
			if not os.path.isfile(fileInewN2HDECAY):
				copyfile(fileInN2HDMCalc, fileInewN2HDECAY)
				print('    File "{}" was created.'.format(fileInewN2HDECAY))
				fileCreated = True
				fileCreatedCounter += 1
			else:
				# Compare the timestamps of the two files to check if the file in the N2HDMCalc folder has changed w.r.t. the file in the ewN2HDECAY folder
				timestampN2HDMCalc = os.path.getmtime(fileInN2HDMCalc)
				timestampewN2HDECAY = os.path.getmtime(fileInewN2HDECAY)
				if timestampN2HDMCalc > timestampewN2HDECAY:
					copyfile(fileInN2HDMCalc, fileInewN2HDECAY)
					print('    File "{}" was updated.'.format(fileInewN2HDECAY))
					fileUpdated = True
					fileUpdatedCounter += 1
				# else:
				# 	print('    File "{}" is already up-to-date.'.format(fileInN2HDMCalc))
		if fileCreated:
			print('    {} files were created.'.format(fileCreatedCounter))
		if fileUpdated:
			print('    {} files were updated.'.format(fileUpdatedCounter))
		if (not fileCreated) and (not fileUpdated):
			print('    All files are already up-to-date. No files were copied.')
		
	# Copy the additional files
	print('  Copying additional files.')
	fileCreated = False
	fileUpdated = False
	fileCreatedCounter = 0
	fileUpdatedCounter = 0
	for fileInN2HDMCalc in additionalFiles:
		fileInewN2HDECAY = rootDirewN2HDECAY + os.sep + fileInN2HDMCalc
		# Check if the file does not exist yet 
		if not os.path.isfile(fileInewN2HDECAY):
			copyfile(fileInN2HDMCalc, fileInewN2HDECAY)
			print('    File "{}" was created.'.format(fileInewN2HDECAY))
			fileCreated = True
			fileCreatedCounter += 1
		else:
			# Compare the timestamps of the two files to check if the file in the N2HDMCalc folder has changed w.r.t. the file in the ewN2HDECAY folder
			timestampN2HDMCalc = os.path.getmtime(fileInN2HDMCalc)
			timestampewN2HDECAY = os.path.getmtime(fileInewN2HDECAY)
			if timestampN2HDMCalc > timestampewN2HDECAY:
				copyfile(fileInN2HDMCalc, fileInewN2HDECAY)
				print('    File "{}" was updated.'.format(fileInewN2HDECAY))
				fileUpdated = True
				fileUpdatedCounter += 1
	if fileCreated:
		print('    {} files were created.'.format(fileCreatedCounter))
	if fileUpdated:
		print('    {} files were updated.'.format(fileUpdatedCounter))
	if (not fileCreated) and (not fileUpdated):
		print('    All files are already up-to-date. No files were copied.')

	# Copy all process files 
	processList = []
	for item in convertList:
		processList.append(item[0])
	fileCreated = False
	fileUpdated = False
	fileCreatedCounter = 0
	fileUpdatedCounter = 0
	print('  Copying process files.')
	for folderPath in os.listdir(processFolder):
		processPathN2HDMCalc = processFolder + os.sep + folderPath
		processPathewN2HDECAY = rootDirewN2HDECAY + os.sep + processPathN2HDMCalc
		for processFile in processFiles:
			processFileN2HDMCalc = processPathN2HDMCalc + os.sep + processFile
			processFileewN2HDECAY = processPathewN2HDECAY + os.sep + processFile

			# Check if the file does not exist yet 
			if not os.path.isfile(processFileewN2HDECAY):
				if processFile == 'processDescription.txt':
					fileHandler = open(processFileN2HDMCalc, "r")
					convertedFile = ''
					lineCount = 1
					for line in fileHandler:
						if lineCount == 1:
							finalParticles = (line.split())[2].split(',')
							if finalParticles[0] == finalParticles[1]:
								symmetryFactor = '2'
							else:
								symmetryFactor = '1'
							processPosition = processList.index((line.split())[0] + ' ' + (line.split())[1] + ' ' + (line.split())[2])
							convertedFile += convertList[processPosition][1] + '\n'
						else:
							convertedFile += line
						lineCount += 1
					if convertList[processPosition][2]:
						convertedFile += '\n1\n' + str(processPosition + 1) + '\n' + symmetryFactor
					else:
						convertedFile += '\n0\n' + '0' + '\n' + symmetryFactor
					fileHandler.close()
					fileHandler = open(processFileewN2HDECAY, "w")
					fileHandler.write(convertedFile)
					fileHandler.close()
				else:
					copyfile(processFileN2HDMCalc, processFileewN2HDECAY)
				print('    File "{}" was created.'.format(processFileewN2HDECAY))
				fileCreated = True 
				fileCreatedCounter += 1
			else:
				# Compare the timestamps of the two files to check if the file in the N2HDMCalc folder has changed w.r.t. the file in the ewN2HDECAY folder
				timestampN2HDMCalc = os.path.getmtime(processFileN2HDMCalc)
				timestampewN2HDECAY = os.path.getmtime(processFileewN2HDECAY)
				if timestampN2HDMCalc > timestampewN2HDECAY:
					if processFile == 'processDescription.txt':
						fileHandler = open(processFileN2HDMCalc, "r")
						convertedFile = ''
						lineCount = 1
						for line in fileHandler:
							if lineCount == 1:
								finalParticles = (line.split())[2].split(',')
								if finalParticles[0] == finalParticles[1]:
									symmetryFactor = '2'
								else:
									symmetryFactor = '1'
								processPosition = processList.index((line.split())[0] + ' ' + (line.split())[1] + ' ' + (line.split())[2])
								convertedFile += convertList[processPosition][1] + '\n'
							else:
								convertedFile += line
							lineCount += 1
						if convertList[processPosition][2]:
							convertedFile += '\n1\n' + str(processPosition + 1) + '\n' + symmetryFactor
						else:
							convertedFile += '\n0\n' + '0' + '\n' + symmetryFactor
						fileHandler.close()
						fileHandler = open(processFileewN2HDECAY, "w")
						fileHandler.write(convertedFile)
						fileHandler.close()
					else:
						copyfile(processFileN2HDMCalc, processFileewN2HDECAY)
					print('    File "{}" was updated.'.format(processFileewN2HDECAY))
					fileUpdated = True 
					fileUpdatedCounter += 1
				# else:
				# 	print('    File "{}" is already up-to-date.'.format(processFileN2HDMCalc))	
	if fileCreated:
		print('    {} files were created.'.format(fileCreatedCounter))
	if fileUpdated:
		print('    {} files were updated.'.format(fileUpdatedCounter))
	if (not fileCreated) and (not fileUpdated):
		print('    All files are already up-to-date. No files were copied.')


	# End message 
	print('\nCopying process completed successfully.\n')