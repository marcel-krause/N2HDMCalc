(* ::Package:: *)

dir = DirectoryName[$InputFileName]; 
ImportPaths = FileNameJoin[{dir, "..", "Paths.m"}]; 
Get[ImportPaths]; 
commandLine = $CommandLine; 
stringlist = StringSplit[commandLine]; 
ampFolder = StringTrim[stringlist[[-1]], "-"][[1]]; 
mOut2 = ToExpression[StringTrim[stringlist[[-2]], "-"][[1]]]; 
mOut1 = ToExpression[StringTrim[stringlist[[-3]], "-"][[1]]]; 
mIn1 = ToExpression[StringTrim[stringlist[[-4]], "-"][[1]]]; 
$LoadFeynArts = $LoadPhi = False; 
Get[FeyncalcDirectory]; 
Get[FeynArtsFeynCalcToolsDirectory];
ImportFileAmp = FileNameJoin[{dir, "..", "BuildingBlocks", "Processes", ampFolder, "TreeLevel", "tree.txt"}]; 
ImportFileAmpConj = FileNameJoin[{dir, "..", "BuildingBlocks", "Processes", ampFolder, "TreeLevel", "treeConj.txt"}]; 
\[Lambda][a_,b_,c_]:=Sqrt[a^2+b^2+c^2-2*a*b-2*a*c-2*b*c];
polarizationSums={EMPTY:>EMPTY};
polarizationSums2={EMPTY:>EMPTY};
ToExpression[Import[FileNameJoin[{dir, "..", "BuildingBlocks", "Processes", ampFolder, "polarizationReplace.txt"}]]];
momentaReplace = {p1:>mIn1, p3:>mOut1, p4:>mOut2};
amp = ToExpression[Import[ImportFileAmp]]; 
ampConj = ToExpression[Import[ImportFileAmpConj]]; 
ampTreeSquared = Simplify[amp*ampConj];
ampTreeSquaredPolSums = (ampTreeSquared /.polarizationSums/.polarizationSums2)/.momentaReplace;
export = (ampTreeSquaredPolSums/.{C0[0,x_,x_,0,0,b_]:>C0Mine[DBLE[0],DBLE[x],DBLE[x],DBLE[0],DBLE[0],DBLE[b]]}); 
exportFortran = "";
exportFortran = exportFortran <> "double precision function " <> ampFolder <> "Tree()\n";
exportFortran = exportFortran <> " use constants\n";
exportFortran = exportFortran <> " implicit none\n";
exportFortran = exportFortran <> "#include \"looptools.h\"\n";
exportFortran = exportFortran <> " double precision :: totalAmplitude\n\n";
exportFortran=exportFortran<>StringReplace[feyncalcToFortran[export," totalAmplitude = "],"Global`"->""]<>"\n\n";
exportFortran=exportFortran<>" " <> ampFolder <> "Tree = totalAmplitude\n";
exportFortran=exportFortran<>"end function " <> ampFolder <> "Tree";
filename = FileNameJoin[{dir, "..", "BuildingBlocks", "Processes", ampFolder, "TreeLevelWidthRed.txt"}];
Export[filename, exportFortran]; 
If[FileExistsQ[StringReplace[filename,".txt"->".F90"]],DeleteFile[StringReplace[filename,".txt"->".F90"]];];
RenameFile[filename,StringReplace[filename,".txt"->".F90"]];
Exit[];