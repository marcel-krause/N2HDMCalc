(* ::Package:: *)

dir = DirectoryName[$InputFileName]; 
ImportPaths = FileNameJoin[{dir, "..", "Paths.m"}]; 
Get[ImportPaths]; 
$LoadFeynArts = $LoadPhi = False; 
Get[FeyncalcDirectory]; 
Get[FeynArtsFeynCalcToolsDirectory];
selfies = {"SelfH1H2Add", "SelfH1H3Add", "SelfH2H3Add", "SelfG0A0Add", "SelfGpHpAdd"};
selfiesExpressions = {
	-EL^2/(32*Pi^2*CW^2*SW^2)*(x-(MH1^2+MH2^2)/2)*((-CA2*(SB*CA1-CB*SA1))*(SA2*SA3*(SB*CA1-CB*SA1)+CA3*(CA1*CB+SA1*SB))*B0[x,MZ^2,MA0^2]+(CA2*(CA1*CB+SA1*SB))*(-SA2*SA3*(CA1*CB+SA1*SB)+CA3*(SB*CA1-CB*SA1))*B0[x,MZ^2,MZ^2]+2*CW^2*((-CA2*(SB*CA1-CB*SA1))*(SA2*SA3*(SB*CA1-CB*SA1)+CA3*(CA1*CB+SA1*SB))*B0[x,MW^2,MHp^2]+(CA2*(CA1*CB+SA1*SB))*(-SA2*SA3*(CA1*CB+SA1*SB)+CA3*(SB*CA1-CB*SA1))*B0[x,MW^2,MW^2])), 
	-EL^2/(32*Pi^2*CW^2*SW^2)*(x-(MH1^2+MH3^2)/2)*((-CA2*(SB*CA1-CB*SA1))*(CA3*SA2*(SB*CA1-CB*SA1)-SA3*(CA1*CB+SA1*SB))*B0[x,MZ^2,MA0^2]+(CA2*(CA1*CB+SA1*SB))*(-CA3*SA2*(CA1*CB+SA1*SB)-SA3*(SB*CA1-CB*SA1))*B0[x,MZ^2,MZ^2]+2*CW^2*((-CA2*(SB*CA1-CB*SA1))*(CA3*SA2*(SB*CA1-CB*SA1)-SA3*(CA1*CB+SA1*SB))*B0[x,MW^2,MHp^2]+(CA2*(CA1*CB+SA1*SB))*(-CA3*SA2*(CA1*CB+SA1*SB)-SA3*(SB*CA1-CB*SA1))*B0[x,MW^2,MW^2])), 
	-EL^2/(32*Pi^2*CW^2*SW^2)*(x-(MH2^2+MH3^2)/2)*((SA2*SA3*(SB*CA1-CB*SA1)+CA3*(CA1*CB+SA1*SB))*(CA3*SA2*(SB*CA1-CB*SA1)-SA3*(CA1*CB+SA1*SB))*B0[x,MZ^2,MA0^2]+(-SA2*SA3*(CA1*CB+SA1*SB)+CA3*(SB*CA1-CB*SA1))*(-CA3*SA2*(CA1*CB+SA1*SB)-SA3*(SB*CA1-CB*SA1))*B0[x,MZ^2,MZ^2]+2*CW^2*((SA2*SA3*(SB*CA1-CB*SA1)+CA3*(CA1*CB+SA1*SB))*(CA3*SA2*(SB*CA1-CB*SA1)-SA3*(CA1*CB+SA1*SB))*B0[x,MW^2,MHp^2]+(-SA2*SA3*(CA1*CB+SA1*SB)+CA3*(SB*CA1-CB*SA1))*(-CA3*SA2*(CA1*CB+SA1*SB)-SA3*(SB*CA1-CB*SA1))*B0[x,MW^2,MW^2])), 
	-EL^2/(32*Pi^2*CW^2*SW^2)*(x-MA0^2/2)*((CA2*(CA1*CB+SA1*SB))*(-CA2*(SB*CA1-CB*SA1))*B0[x,MZ^2,MH1^2]+(-SA2*SA3*(CA1*CB+SA1*SB)+CA3*(SB*CA1-CB*SA1))*(SA2*SA3*(SB*CA1-CB*SA1)+CA3*(CA1*CB+SA1*SB))*B0[x,MZ^2,MH2^2]+(-CA3*SA2*(CA1*CB+SA1*SB)-SA3*(SB*CA1-CB*SA1))*(CA3*SA2*(SB*CA1-CB*SA1)-SA3*(CA1*CB+SA1*SB))*B0[x,MZ^2,MH3^2]), 
	-EL^2/(16*Pi^2*SW^2)*(x-MHp^2/2)*((CA2*(CA1*CB+SA1*SB))*(-CA2*(SB*CA1-CB*SA1))*B0[x,MW^2,MH1^2]+(-SA2*SA3*(CA1*CB+SA1*SB)+CA3*(SB*CA1-CB*SA1))*(SA2*SA3*(SB*CA1-CB*SA1)+CA3*(CA1*CB+SA1*SB))*B0[x,MW^2,MH2^2]+(-CA3*SA2*(CA1*CB+SA1*SB)-SA3*(SB*CA1-CB*SA1))*(CA3*SA2*(SB*CA1-CB*SA1)-SA3*(CA1*CB+SA1*SB))*B0[x,MW^2,MH3^2])
};
For[counter=1,counter<=Length[selfies],counter++,
	Print["Calculating additional self-energy contribution "<> selfies[[counter]] <> " ..."];
	export = 0; 
	exportFortran = "";
	exportFortran = exportFortran <> "double complex function " <> selfies[[counter]] <> "(x)\n";
	exportFortran = exportFortran <> " use constants\n";
	exportFortran = exportFortran <> " implicit none\n";
	exportFortran = exportFortran <> "#include \"looptools.h\"\n";
	exportFortran = exportFortran <> " double precision, intent(in) :: x\n";
	exportFortran = exportFortran <> " double complex :: totalAmplitude\n\n";
	export = selfiesExpressions[[counter]];
	exportFortran=exportFortran<>feyncalcToFortran[export," totalAmplitude = "]<>"\n\n";
	exportFortran=exportFortran<>" " <> selfies[[counter]] <> " = totalAmplitude\n";
	exportFortran=exportFortran<>"end function " <> selfies[[counter]] <> "\n\n";
	filename = FileNameJoin[{dir, "..", "BuildingBlocks", "SelfEnergies", "Alternative", selfies[[counter]] <> ".txt"}];
	Export[filename, exportFortran]; 
	If[FileExistsQ[StringReplace[filename,".txt"->".F90"]],DeleteFile[StringReplace[filename,".txt"->".F90"]];];
	RenameFile[filename,StringReplace[filename,".txt"->".F90"]];
	Print["Additional self-energy contribution "<> selfies[[counter]] <> " finished."];
	Clear[FAamp, simpler, result, export, filename, exportFortran]; 
];
Exit[]; 
