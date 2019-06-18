(* ::Package:: *)

dir = DirectoryName[$InputFileName]; 
ImportPaths = FileNameJoin[{dir, "..", "Paths.m"}]; 
Get[ImportPaths]; 
commandLine = $CommandLine;
stringlist = StringSplit[commandLine];
process = StringTrim[stringlist[[-1]], "-"][[1]];
inputList = StringSplit[process, "to"];
incomingList = StringSplit[inputList[[1]], ","];
outgoingList = StringSplit[inputList[[2]], ","];
incomingAddCount=StringCount[inputList[[1]],"{"];
incomingCount=Length[incomingList]-incomingAddCount;
outgoingAddCount=StringCount[inputList[[2]],"{"];
outgoingCount=Length[outgoingList]-outgoingAddCount;
processString = ToString[incomingList] <> "->" <> ToString[outgoingList];
$CKM = True; 
Get[FeynartsDirectory]; 
topology = CreateTopologies[0, incomingCount -> outgoingCount, ExcludeTopologies -> {Internal}]; 
AA = InsertFields[topology, ToExpression[processString], Model -> "N2HDMEWSB", InsertionLevel -> Particles]; 
amp = CreateFeynAmp[AA, GaugeRules -> {}]; 
amplitudes = PickLevel[Particles][amp]; 
amplituderange = Range[1, Length[amplitudes]]; 
FAamp = ""; 
ToExpression[Import[FeynartsparserDirectory]]; 
ampString = StringJoin["{", FAamp, "}"]; 
ampList = ToExpression[ampString]; 
Export[FileNameJoin[{dir,"..","Temp","tree",StringJoin["amp.txt"]}], ampList[[1]]]; 
Exit[]; 
