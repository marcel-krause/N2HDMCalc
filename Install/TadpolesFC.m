(* ::Package:: *)

dir = DirectoryName[$InputFileName]; 
ImportPaths = FileNameJoin[{dir, "..", "Paths.m"}]; 
Get[ImportPaths]; 
$LoadFeynArts = $LoadPhi = False; 
Get[FeyncalcDirectory]; 
Get[FeynArtsFeynCalcToolsDirectory];
selfies := scalarTadpolesFAamp; 
saver := scalarTadpolesSaver; 
targetfile = "ScalarTadpolesFA.txt";
getParticleContent[ParticleContentFile,FileNameJoin[{dir, "..", "BuildingBlocks", "Tadpoles"}]];
ImportFile = FileNameJoin[{dir, "..", "BuildingBlocks", "Tadpoles", targetfile}]; 
ToExpression[StringReplace[Import[ImportFile],{"Sqrt[1/;1>1]*"->"","Sqrt[1/;2>1]*"->"","*Sqrt[1/;1>1]"->"","*Sqrt[1/;2>1]"->"","Sqrt[1 /; 1 > 1]*"->"","Sqrt[1 /; 2 > 1]*"->"","*Sqrt[1 /; 1 > 1]"->"","*Sqrt[1 /; 2 > 1]"->""}]]; 
p1 = p3 = k; 
For[counter=1,counter<=Length[selfies],counter++,
Print["Calculating tadpole "<> selfies[[counter]] <> " ..."];
export = 0; 
exportFortran = "";
exportFortran = exportFortran <> "double complex function " <> selfies[[counter]] <> "()\n";
exportFortran = exportFortran <> " use constants\n";
exportFortran = exportFortran <> " implicit none\n";
exportFortran = exportFortran <> "#include \"looptools.h\"\n";
exportFortran = exportFortran <> " integer :: j\n";
exportFortran = exportFortran <> " double complex :: totalAmplitude\n";
exportFortran = exportFortran <> " double complex :: amplitudes("<>ToString[Length[ToExpression[selfies[[counter]]]]]<>")\n\n";
RR={{CA1*CA2,CA2*SA1,SA2},{-CA3*SA1-CA1*SA2*SA3,CA1*CA3-SA1*SA2*SA3,CA2*SA3},{-CA1*CA3*SA2+SA1*SA3,-CA3*SA1*SA2-CA1*SA3,CA2*CA3}};
RA={{CB,SB},{-SB,CB}};
For[i = 1, i <= Length[ToExpression[selfies[[counter]]]], i++, 
 Print["Calculating sub-amplitude number "<>ToString[i]<>" ..."];
  (*fullAmp=ToExpression[selfies[[counter]]][[i]]/.{MassAh[1]:>Sqrt[GaugeXi[Z]]*MZ,MassHm[1]:>Sqrt[GaugeXi[W]]*MW};*)
	fullAmp=ToExpression[selfies[[counter]]][[i]]/.{MassAh[1]:>MZ,MassHm[1]:>MW};
	tracer = (((fullAmp)/.{PropagatorDenominator[0,0]:>PropagatorDenominator[0,\[Mu]]})/.{DiracTrace -> Tr}); If[tracer == 0, exportFortran=exportFortran<>" amplitudes("<>ToString[i]<>") = 0D0"<>"\n\n"; Continue[]; ]; If[FreeQ[tracer, l], exportFortran=exportFortran<>" amplitudes("<>ToString[i]<>") = 0D0"<>"\n\n"; Continue[]; ];
	simpler = OneLoopSimplify[tracer, l]; oneloop = OneLoop[l, simpler]; 
	(*result = Simplify[PaVeReduce[oneloop /. {Pair[Momentum[k], Momentum[k]] -> x}]]/.{ZEL[a_,a_]:>1,ZER[a_,a_]:>1}/.{ZEL[__]:>0,ZER[__]:>0}/.{ZUL[a_,a_]:>1,ZUR[a_,a_]:>1}/.{ZUL[__]:>0,ZUR[__]:>0}/.{Conjugate[Ye[1,1]]:>Sqrt[2]*ME/(2*MW*SW/EL),Conjugate[Ye[2,2]]:>Sqrt[2]*MM/(2*MW*SW/EL),Conjugate[Ye[3,3]]:>Sqrt[2]*ML/(2*MW*SW/EL),Conjugate[Yu[1,1]]:>Sqrt[2]*MU/(2*MW*SW/EL),Conjugate[Yu[2,2]]:>Sqrt[2]*MC/(2*MW*SW/EL),Conjugate[Yu[3,3]]:>Sqrt[2]*MT/(2*MW*SW/EL),Conjugate[Yd[1,1]]:>Sqrt[2]*MD/(2*MW*SW/EL),Conjugate[Yd[2,2]]:>Sqrt[2]*MS/(2*MW*SW/EL),Conjugate[Yd[3,3]]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{Ye[1,1]:>Sqrt[2]*ME/(2*MW*SW/EL),Ye[2,2]:>Sqrt[2]*MM/(2*MW*SW/EL),Ye[3,3]:>Sqrt[2]*ML/(2*MW*SW/EL),Yu[1,1]:>Sqrt[2]*MU/(2*MW*SW/EL),Yu[2,2]:>Sqrt[2]*MC/(2*MW*SW/EL),Yu[3,3]:>Sqrt[2]*MT/(2*MW*SW/EL),Yd[1,1]:>Sqrt[2]*MD/(2*MW*SW/EL),Yd[2,2]:>Sqrt[2]*MS/(2*MW*SW/EL),Yd[3,3]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{MassFe[1]:>ME,MassFe[2]:>MM,MassFe[3]:>ML,MassFv[_]:>0,MassFu[1]:>MU,MassFu[2]:>MC,MassFu[3]:>MT,MassFd[1]:>MD,MassFd[2]:>MS,MassFd[3]:>MB}/.{vd:>CB*(2*MW*SW/EL),vu:>SB*(2*MW*SW/EL)}/.{Lam1:>(-mu2*SB^2+MH1^2*RR[[1]][[1]]^2+MH2^2*RR[[2]][[1]]^2+MH3^2*RR[[3]][[1]]^2)/(CB^2*v^2),Lam2:>(-mu2*CB^2+MH1^2*RR[[1]][[2]]^2+MH2^2*RR[[2]][[2]]^2+MH3^2*RR[[3]][[2]]^2)/(SB^2*v^2),Lam3:>(-mu2+2*MHp^2+(MH1^2*RR[[1]][[1]]*RR[[1]][[2]]+MH2^2*RR[[2]][[1]]*RR[[2]][[2]]+MH3^2*RR[[3]][[1]]*RR[[3]][[2]])/SB/CB)/v^2,Lam4:>(mu2+MA0^2-2*MHp^2)/v^2,Lam5:>(mu2-MA0^2)/v^2,Lam6:>(MH1^2*RR[[1]][[3]]^2+MH2^2*RR[[2]][[3]]^2+MH3^2*RR[[3]][[3]]^2)/vS^2,Lam7:>(MH1^2*RR[[1]][[1]]*RR[[1]][[3]]+MH2^2*RR[[2]][[1]]*RR[[2]][[3]]+MH3^2*RR[[3]][[1]]*RR[[3]][[3]])/(v*vS*CB),Lam8:>(MH1^2*RR[[1]][[2]]*RR[[1]][[3]]+MH2^2*RR[[2]][[2]]*RR[[2]][[3]]+MH3^2*RR[[3]][[2]]*RR[[3]][[3]])/(v*vS*SB)}/.{Masshh[1]:>MH1,Masshh[2]:>MH2,Masshh[3]:>MH3,MassAh[2]:>MA0,MassHm[2]:>MHp}/.{ZH[aa_,bb_]:>RR[[aa]][[bb]]}/.{ZA[aa_,bb_]:>RA[[aa]][[bb]],ZP[aa_,bb_]:>RA[[aa]][[bb]]}/.{v:>(2*MW*SW/EL)}/.{mu2:>m12squared/SB/CB}/.{YukS1Lep[1]:>YukS1Lep1,YukS1Lep[2]:>YukS1Lep2,YukS1Lep[3]:>YukS1Lep3,YukS2Lep[1]:>YukS2Lep1,YukS2Lep[2]:>YukS2Lep2,YukS3Lep[1]:>YukS3Lep1,YukS3Lep[2]:>YukS3Lep2}/.{YukS1Quark[1]:>YukS1Quark1,YukS1Quark[2]:>YukS1Quark2,YukS1Quark[3]:>YukS1Quark3,YukS2Quark[1]:>YukS2Quark1,YukS2Quark[2]:>YukS2Quark2,YukS3Quark[1]:>YukS3Quark1,YukS3Quark[2]:>YukS3Quark2}/.{GaugeXi[U[2]]:>GaugeXiZ,GaugeXi[U[3]]:>GaugeXiW,GaugeXi[U[4]]:>GaugeXiW}; *)
	result = Simplify[PaVeReduce[oneloop /. {Pair[Momentum[k], Momentum[k]] -> x}]]/.{ZEL[a_,a_]:>1,ZER[a_,a_]:>1}/.{ZEL[__]:>0,ZER[__]:>0}/.{ZUL[a_,a_]:>1,ZUR[a_,a_]:>1}/.{ZUL[__]:>0,ZUR[__]:>0}/.{Conjugate[Ye[1,1]]:>Sqrt[2]*ME/(2*MW*SW/EL),Conjugate[Ye[2,2]]:>Sqrt[2]*MM/(2*MW*SW/EL),Conjugate[Ye[3,3]]:>Sqrt[2]*ML/(2*MW*SW/EL),Conjugate[Yu[1,1]]:>Sqrt[2]*MU/(2*MW*SW/EL),Conjugate[Yu[2,2]]:>Sqrt[2]*MC/(2*MW*SW/EL),Conjugate[Yu[3,3]]:>Sqrt[2]*MT/(2*MW*SW/EL),Conjugate[Yd[1,1]]:>Sqrt[2]*MD/(2*MW*SW/EL),Conjugate[Yd[2,2]]:>Sqrt[2]*MS/(2*MW*SW/EL),Conjugate[Yd[3,3]]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{Ye[1,1]:>Sqrt[2]*ME/(2*MW*SW/EL),Ye[2,2]:>Sqrt[2]*MM/(2*MW*SW/EL),Ye[3,3]:>Sqrt[2]*ML/(2*MW*SW/EL),Yu[1,1]:>Sqrt[2]*MU/(2*MW*SW/EL),Yu[2,2]:>Sqrt[2]*MC/(2*MW*SW/EL),Yu[3,3]:>Sqrt[2]*MT/(2*MW*SW/EL),Yd[1,1]:>Sqrt[2]*MD/(2*MW*SW/EL),Yd[2,2]:>Sqrt[2]*MS/(2*MW*SW/EL),Yd[3,3]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{MassFe[1]:>ME,MassFe[2]:>MM,MassFe[3]:>ML,MassFv[_]:>0,MassFu[1]:>MU,MassFu[2]:>MC,MassFu[3]:>MT,MassFd[1]:>MD,MassFd[2]:>MS,MassFd[3]:>MB}/.{vd:>CB*(2*MW*SW/EL),vu:>SB*(2*MW*SW/EL)}/.{Lam1:>(-mu2*SB^2+MH1^2*RR[[1]][[1]]^2+MH2^2*RR[[2]][[1]]^2+MH3^2*RR[[3]][[1]]^2)/(CB^2*v^2),Lam2:>(-mu2*CB^2+MH1^2*RR[[1]][[2]]^2+MH2^2*RR[[2]][[2]]^2+MH3^2*RR[[3]][[2]]^2)/(SB^2*v^2),Lam3:>(-mu2+2*MHp^2+(MH1^2*RR[[1]][[1]]*RR[[1]][[2]]+MH2^2*RR[[2]][[1]]*RR[[2]][[2]]+MH3^2*RR[[3]][[1]]*RR[[3]][[2]])/SB/CB)/v^2,Lam4:>(mu2+MA0^2-2*MHp^2)/v^2,Lam5:>(mu2-MA0^2)/v^2,Lam6:>(MH1^2*RR[[1]][[3]]^2+MH2^2*RR[[2]][[3]]^2+MH3^2*RR[[3]][[3]]^2)/vS^2,Lam7:>(MH1^2*RR[[1]][[1]]*RR[[1]][[3]]+MH2^2*RR[[2]][[1]]*RR[[2]][[3]]+MH3^2*RR[[3]][[1]]*RR[[3]][[3]])/(v*vS*CB),Lam8:>(MH1^2*RR[[1]][[2]]*RR[[1]][[3]]+MH2^2*RR[[2]][[2]]*RR[[2]][[3]]+MH3^2*RR[[3]][[2]]*RR[[3]][[3]])/(v*vS*SB)}/.{Masshh[1]:>MH1,Masshh[2]:>MH2,Masshh[3]:>MH3,MassAh[2]:>MA0,MassHm[2]:>MHp}/.{ZH[aa_,bb_]:>RR[[aa]][[bb]]}/.{ZA[aa_,bb_]:>RA[[aa]][[bb]],ZP[aa_,bb_]:>RA[[aa]][[bb]]}/.{v:>(2*MW*SW/EL)}/.{mu2:>m12squared/SB/CB}/.{YukS1Lep[1]:>YukS1Lep1,YukS1Lep[2]:>YukS1Lep2,YukS1Lep[3]:>YukS1Lep3,YukS2Lep[1]:>YukS2Lep1,YukS2Lep[2]:>YukS2Lep2,YukS3Lep[1]:>YukS3Lep1,YukS3Lep[2]:>YukS3Lep2}/.{YukS1Quark[1]:>YukS1Quark1,YukS1Quark[2]:>YukS1Quark2,YukS1Quark[3]:>YukS1Quark3,YukS2Quark[1]:>YukS2Quark1,YukS2Quark[2]:>YukS2Quark2,YukS3Quark[1]:>YukS3Quark1,YukS3Quark[2]:>YukS3Quark2}/.{GaugeXi[U[2]]:>1,GaugeXi[U[3]]:>1,GaugeXi[U[4]]:>1}/.{CS1S1S1[gt1_,gt2_,gt3_]:>ToExpression["CS1S1S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]],CS1S3S3[gt1_,gt2_,gt3_]:>ToExpression["CS1S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]],CS2S2S1[gt1_,gt2_,gt3_]:>ToExpression["CS2S2S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]],CS2S2S2S2[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS2S2S2S2f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS2S2S1S1[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS2S2S1S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS2S2S3S3[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS2S2S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS1S1S1S1[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS1S1S1S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS1S1S3S3[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS1S1S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS3S3S3S3[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS3S3S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]]}; 
	export = (result/.{C0[0,x_,x_,0,0,b_]:>C0Mine[DBLE[0],DBLE[x],DBLE[x],DBLE[0],DBLE[0],DBLE[b]]}); exportFortran=exportFortran<>feyncalcToFortran[export," amplitudes("<>ToString[i]<>") = "]<>"\n\n";
	Print["Sub-amplitude number "<>ToString[i]<>" finished."]; 
	ClearAll[tracer, fullAmp, simpler, oneloop, result, export]; ];
exportFortran=exportFortran<>"  totalAmplitude = (0D0,0D0)\n";
exportFortran=exportFortran<>" do j=1,"<>ToString[Length[ToExpression[selfies[[counter]]]]]<>"\n";
exportFortran=exportFortran<>"  totalAmplitude = totalAmplitude + amplitudes(j)\n";
exportFortran=exportFortran<>" end do\n";
exportFortran=exportFortran<>" " <> selfies[[counter]] <> " = totalAmplitude\n";
exportFortran=exportFortran<>"end function " <> selfies[[counter]] <> "\n\n";
filename = saver[[counter]]; 
Export[filename, exportFortran]; 
If[FileExistsQ[StringReplace[filename,".txt"->".F90"]],DeleteFile[StringReplace[filename,".txt"->".F90"]];];
RenameFile[filename,StringReplace[filename,".txt"->".F90"]];
Print["Tadpole "<> selfies[[counter]] <> " finished."];
Clear[FAamp, simpler, result, export, filename, exportFortran]; 
]
DeleteFile[ImportFile];
Exit[]; 
