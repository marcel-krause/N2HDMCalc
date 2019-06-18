(* ::Package:: *)

dir = DirectoryName[$InputFileName]; 
ImportPaths = FileNameJoin[{dir, "..", "Paths.m"}]; 
Get[ImportPaths]; 
$CKM = True; 
Get[FeynartsDirectory]; 
Get[FeynArtsFeynCalcToolsDirectory];
export = {}; 
topology = CreateTopologies[1, 1 -> 1, ExcludeTopologies -> {Internal}];
AA = InsertFields[topology, V[1] -> V[1], Model -> "N2HDMEWSB", LastSelections->{F,!F[3,{3}]}, InsertionLevel -> Particles];
(*amp = CreateFeynAmp[AA, GaugeRules -> {}]; *)
amp = CreateFeynAmp[AA]; 
amplitudes = PickLevel[Particles][amp]; 
amplituderange = Range[1, Length[amplitudes]]; 
FAamp = ""; 
ToExpression[Import[FeynartsparserDirectory]]; 
export = Append[export, StringJoin["SelfAA={", FAamp, "}"]]; 
ExportFile = FileNameJoin[{dir, "..", "BuildingBlocks", "SelfEnergies", "Usual", "SelfEnergyAALightFA.txt"}]; 
Export[ExportFile, export]; 
Exit[]; 
