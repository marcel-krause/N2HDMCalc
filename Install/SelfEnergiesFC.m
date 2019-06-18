(* ::Package:: *)

dir = DirectoryName[$InputFileName]; 
commandLine = StringSplit[$CommandLine[[4]],".m "][[2]];
argv = StringSplit[commandLine, " "];
ImportPaths = FileNameJoin[{dir, "..", "Paths.m"}]; 
Get[ImportPaths]; 
$LoadFeynArts = $LoadPhi = False; 
Get[FeyncalcDirectory]; 
Get[FeynArtsFeynCalcToolsDirectory];
SetOptions[B1,BReduce->False];
RR={{CA1*CA2,CA2*SA1,SA2},{-CA3*SA1-CA1*SA2*SA3,CA1*CA3-SA1*SA2*SA3,CA2*SA3},{-CA1*CA3*SA2+SA1*SA3,-CA3*SA1*SA2-CA1*SA3,CA2*CA3}};
RA={{CB,SB},{-SB,CB}};
Switch[argv[[1]],
"1", selfies := scalarSelfEnergiesFAamp; saver := scalarSelfEnergiesSaver; saverDeriv := scalarSelfEnergiesDerivSaver; targetfile = "SelfEnergiesScalarFA.txt"; fermionParts = {""};,
"2", selfies := vectorSelfEnergiesFAamp; saver := vectorSelfEnergiesSaver; saverDeriv := vectorSelfEnergiesDerivSaver; targetfile = "SelfEnergiesVectorFA.txt"; fermionParts = {""};,
"3", selfies := fermionSelfEnergiesFAamp; saver := fermionSelfEnergiesSaver; saverDeriv := fermionSelfEnergiesDerivSaver; targetfile = "SelfEnergiesFermionFA.txt"; fermionParts = {"Left", "Right", "Scalar"};,
_,Exit[];]
Switch[argv[[2]],
"1", targetDir = "Usual";,
"2", targetDir = "Alternative";,
_, Exit[];]
getParticleContent[ParticleContentFile,FileNameJoin[{dir, "..", "BuildingBlocks", "SelfEnergies", targetDir}]];
ampVersion = StringTake[targetDir,5];
ImportFile = FileNameJoin[{dir, "..", "BuildingBlocks", "SelfEnergies", targetDir, targetfile}]; 
(*ToExpression[StringReplace[StringReplace[Import[ImportFile],{"Sqrt[1/;1>1]*"->"","Sqrt[1/;2>1]*"->"","*Sqrt[1/;1>1]"->"","*Sqrt[1/;2>1]"->"","Sqrt[1 /; 1 > 1]*"->"","Sqrt[1 /; 2 > 1]*"->"","*Sqrt[1 /; 1 > 1]"->"","*Sqrt[1 /; 2 > 1]"->""}],{"MassFe[1]"->"ME","MassFe[2]"->"MM","MassFe[3]"->"ML","MassFv[_]"->"0","MassFu[1]"->"MU","MassFu[2]"->"MC","MassFu[3]"->"MT","MassFd[1]"->"MD","MassFd[2]"->"MS","MassFd[3]"->"MB","MassAh[1]"->"Sqrt[GaugeXi[Z]]*MZ","MassHm[1]"->"Sqrt[GaugeXi[W]]*MW"}]]; *)
ToExpression[StringReplace[StringReplace[Import[ImportFile],{"Sqrt[1/;1>1]*"->"","Sqrt[1/;2>1]*"->"","*Sqrt[1/;1>1]"->"","*Sqrt[1/;2>1]"->"","Sqrt[1 /; 1 > 1]*"->"","Sqrt[1 /; 2 > 1]*"->"","*Sqrt[1 /; 1 > 1]"->"","*Sqrt[1 /; 2 > 1]"->""}],{"MassFe[1]"->"ME","MassFe[2]"->"MM","MassFe[3]"->"ML","MassFv[_]"->"0","MassFu[1]"->"MU","MassFu[2]"->"MC","MassFu[3]"->"MT","MassFd[1]"->"MD","MassFd[2]"->"MS","MassFd[3]"->"MB","MassAh[1]"->"MZ","MassHm[1]"->"MW"}]]; 
p1 = p3 = k; 
If[argv[[1]]=="3",
	counterIsZeroMom = True;
,
	counterIsZeroMom = False;
];
For[counter=1,counter<=Length[ToExpression[selfies]],counter++,
	For[n=1,n<=Length[fermionParts],n++,
		If[counterIsZeroMom,
			If[argv[[1]]=="3",
				zeroMomDesc="";
				zeroMomReplace={Pair[Momentum[k], Momentum[k]] :> Pair[Momentum[k], Momentum[k]]};
			,
				zeroMomDesc="ZeroMom";
				zeroMomReplace={Pair[Momentum[k], Momentum[k]] -> 0};
			];
		,
			zeroMomReplace={Pair[Momentum[k], Momentum[k]] :> Pair[Momentum[k], Momentum[k]]};
			zeroMomDesc="";
		];
		Print["Calculating self-energy "<> selfies[[counter]] <> ampVersion <> zeroMomDesc <> " ..."];
		export = 0; 
		exportFortran = "";
		exportFortran = exportFortran <> "double complex function " <> selfies[[counter]] <> fermionParts[[n]] <> ampVersion <> zeroMomDesc <> "(x)\n";
		exportFortran = exportFortran <> " use constants\n";
		exportFortran = exportFortran <> " implicit none\n";
		exportFortran = exportFortran <> "#include \"looptools.h\"\n";
		exportFortran = exportFortran <> " double precision, intent(in) :: x\n";
		exportFortran = exportFortran <> " integer :: j\n";
		exportFortran = exportFortran <> " double complex :: totalAmplitude\n";
		exportFortran = exportFortran <> " double complex :: amplitudes("<>ToString[Length[ToExpression[selfies[[counter]]]]]<>")\n\n";
		exportFortranDeriv = "";
		exportFortranDeriv = exportFortranDeriv <> "double complex function D" <> selfies[[counter]] <> fermionParts[[n]] <> "(x)\n";
		exportFortranDeriv = exportFortranDeriv <> " use constants\n";
		exportFortranDeriv = exportFortranDeriv <> " implicit none\n";
		exportFortranDeriv = exportFortranDeriv <> "#include \"looptools.h\"\n";
		exportFortranDeriv = exportFortranDeriv <> " double precision, intent(in) :: x\n";
		exportFortranDeriv = exportFortranDeriv <> " integer :: j\n";
		exportFortranDeriv = exportFortranDeriv <> " double complex :: totalAmplitude\n";
		exportFortranDeriv = exportFortranDeriv <> " double complex :: amplitudes("<>ToString[Length[ToExpression[selfies[[counter]]]]]<>")\n\n";
		For[i = 1, i <= Length[ToExpression[selfies[[counter]]]], i++, 
			Print["Calculating sub-amplitude number "<>ToString[i]<>" ..."];
			If[argv[[1]]=="3",
				truncated={Spinor[Momentum[p_],m_,1]:>1,Spinor[Momentum[p_],m_,1]:>1};
				tracer = ((((ToExpression[selfies[[counter]]][[i]])/.{PropagatorDenominator[0,0]:>PropagatorDenominator[0,\[Mu]]})/.truncated)/.{DiracTrace -> Tr});
				If[tracer == 0, 
					exportFortran=exportFortran<>" amplitudes("<>ToString[i]<>") = 0D0"<>"\n\n"; 
					Print["Sub-amplitude number "<>ToString[i]<>" finished."];
					Continue[]; 
				]; 
			,
				tracer = (((ToExpression[selfies[[counter]]][[i]])/.{PropagatorDenominator[0,0]:>PropagatorDenominator[0,\[Mu]]})/.{DiracTrace -> Tr}); 
			];
		If[tracer == 0, 
			exportFortran=exportFortran<>" amplitudes("<>ToString[i]<>") = 0D0"<>"\n\n"; 
			If[argv[[2]]=="1",
				exportFortranDeriv=exportFortranDeriv<>" amplitudes("<>ToString[i]<>") = 0D0"<>"\n\n"; 
			,
				False;
			];
			Print["Sub-amplitude number "<>ToString[i]<>" finished."];
			Continue[]; 
		]; 
		If[FreeQ[tracer, l], 
			exportFortran=exportFortran<>" amplitudes("<>ToString[i]<>") = 0D0"<>"\n\n"; Continue[]; 
		];
		If[counterIsZeroMom,
			simpler=tracer;
		,
			simpler = OneLoopSimplify[tracer, l];
		]; 
		If[argv[[1]]=="2",
			If[counter==1,
				simpler = tracer;
			,
				False;
			];
		,
			False;
		];
		oneloop = OneLoop[l, simpler]/.zeroMomReplace; 
		If[argv[[1]]=="2",
			If[counter==1,
				result = TrigExpand[PaVeReduce[oneloop /. {Pair[Momentum[k], Momentum[k]] -> x}]];,
				result = Simplify[PaVeReduce[oneloop /. {Pair[Momentum[k], Momentum[k]] -> x}]];
			];,
			result = Simplify[PaVeReduce[oneloop /. {Pair[Momentum[k], Momentum[k]] -> x}]];
		];
		If[argv[[1]]=="3",
			str=ToString[ToExpression[selfies[[counter]]],InputForm];
			begin=StringPosition[str,"Spinor[Momentum[k],"][[1]][[2]]+2;
			end=StringPosition[str,"Spinor[Momentum[k],"][[1]][[2]]+7;
			mass=ToExpression[StringSplit[StringTake[str,{begin,end}],","][[1]]];
			beginSecond=StringPosition[str,"Spinor[Momentum[k],"][[2]][[2]]+2;
			endSecond=StringPosition[str,"Spinor[Momentum[k],"][[2]][[2]]+7;
			massSecond=ToExpression[StringSplit[StringTake[str,{beginSecond,endSecond}],","][[1]]];
			If[mass==0,
				mass=1;
			,
				False;
			];
			(*resultTemp = result/.{IndexDelta[__]:>1,SumOver[__]:>1}/.{ZEL[a_,a_]:>1,ZER[a_,a_]:>1}/.{ZEL[__]:>0,ZER[__]:>0}/.{ZUL[a_,a_]:>1,ZUR[a_,a_]:>1}/.{ZUL[__]:>0,ZUR[__]:>0}/.{Conjugate[Ye[1,1]]:>Sqrt[2]*ME/(2*MW*SW/EL),Conjugate[Ye[2,2]]:>Sqrt[2]*MM/(2*MW*SW/EL),Conjugate[Ye[3,3]]:>Sqrt[2]*ML/(2*MW*SW/EL),Conjugate[Yu[1,1]]:>Sqrt[2]*MU/(2*MW*SW/EL),Conjugate[Yu[2,2]]:>Sqrt[2]*MC/(2*MW*SW/EL),Conjugate[Yu[3,3]]:>Sqrt[2]*MT/(2*MW*SW/EL),Conjugate[Yd[1,1]]:>Sqrt[2]*MD/(2*MW*SW/EL),Conjugate[Yd[2,2]]:>Sqrt[2]*MS/(2*MW*SW/EL),Conjugate[Yd[3,3]]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{Ye[1,1]:>Sqrt[2]*ME/(2*MW*SW/EL),Ye[2,2]:>Sqrt[2]*MM/(2*MW*SW/EL),Ye[3,3]:>Sqrt[2]*ML/(2*MW*SW/EL),Yu[1,1]:>Sqrt[2]*MU/(2*MW*SW/EL),Yu[2,2]:>Sqrt[2]*MC/(2*MW*SW/EL),Yu[3,3]:>Sqrt[2]*MT/(2*MW*SW/EL),Yd[1,1]:>Sqrt[2]*MD/(2*MW*SW/EL),Yd[2,2]:>Sqrt[2]*MS/(2*MW*SW/EL),Yd[3,3]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{MassFe[1]:>ME,MassFe[2]:>MM,MassFe[3]:>ML,MassFv[_]:>0,MassFu[1]:>MU,MassFu[2]:>MC,MassFu[3]:>MT,MassFd[1]:>MD,MassFd[2]:>MS,MassFd[3]:>MB}/.{vd:>CB*(2*MW*SW/EL),vu:>SB*(2*MW*SW/EL)}/.{Lam1:>(-mu2*SB^2+MH1^2*RR[[1]][[1]]^2+MH2^2*RR[[2]][[1]]^2+MH3^2*RR[[3]][[1]]^2)/(CB^2*v^2),Lam2:>(-mu2*CB^2+MH1^2*RR[[1]][[2]]^2+MH2^2*RR[[2]][[2]]^2+MH3^2*RR[[3]][[2]]^2)/(SB^2*v^2),Lam3:>(-mu2+2*MHp^2+(MH1^2*RR[[1]][[1]]*RR[[1]][[2]]+MH2^2*RR[[2]][[1]]*RR[[2]][[2]]+MH3^2*RR[[3]][[1]]*RR[[3]][[2]])/SB/CB)/v^2,Lam4:>(mu2+MA0^2-2*MHp^2)/v^2,Lam5:>(mu2-MA0^2)/v^2,Lam6:>(MH1^2*RR[[1]][[3]]^2+MH2^2*RR[[2]][[3]]^2+MH3^2*RR[[3]][[3]]^2)/vS^2,Lam7:>(MH1^2*RR[[1]][[1]]*RR[[1]][[3]]+MH2^2*RR[[2]][[1]]*RR[[2]][[3]]+MH3^2*RR[[3]][[1]]*RR[[3]][[3]])/(v*vS*CB),Lam8:>(MH1^2*RR[[1]][[2]]*RR[[1]][[3]]+MH2^2*RR[[2]][[2]]*RR[[2]][[3]]+MH3^2*RR[[3]][[2]]*RR[[3]][[3]])/(v*vS*SB)}/.{Masshh[1]:>MH1,Masshh[2]:>MH2,Masshh[3]:>MH3,MassAh[2]:>MA0,MassHm[2]:>MHp}/.{ZH[aa_,bb_]:>RR[[aa]][[bb]]}/.{ZA[aa_,bb_]:>RA[[aa]][[bb]],ZP[aa_,bb_]:>RA[[aa]][[bb]]}/.{v:>(2*MW*SW/EL)}/.{mu2:>m12squared/SB/CB}/.{YukS1Lep[1]:>YukS1Lep1,YukS1Lep[2]:>YukS1Lep2,YukS1Lep[3]:>YukS1Lep3,YukS2Lep[1]:>YukS2Lep1,YukS2Lep[2]:>YukS2Lep2,YukS3Lep[1]:>YukS3Lep1,YukS3Lep[2]:>YukS3Lep2}/.{YukS1Quark[1]:>YukS1Quark1,YukS1Quark[2]:>YukS1Quark2,YukS1Quark[3]:>YukS1Quark3,YukS2Quark[1]:>YukS2Quark1,YukS2Quark[2]:>YukS2Quark2,YukS3Quark[1]:>YukS3Quark1,YukS3Quark[2]:>YukS3Quark2}/.{GaugeXi[U[2]]:>GaugeXiZ,GaugeXi[U[3]]:>GaugeXiW,GaugeXi[U[4]]:>GaugeXiW};*)
			resultTemp = result/.{IndexDelta[__]:>1,SumOver[__]:>1}/.{ZEL[a_,a_]:>1,ZER[a_,a_]:>1}/.{ZEL[__]:>0,ZER[__]:>0}/.{ZUL[a_,a_]:>1,ZUR[a_,a_]:>1}/.{ZUL[__]:>0,ZUR[__]:>0}/.{Conjugate[Ye[1,1]]:>Sqrt[2]*ME/(2*MW*SW/EL),Conjugate[Ye[2,2]]:>Sqrt[2]*MM/(2*MW*SW/EL),Conjugate[Ye[3,3]]:>Sqrt[2]*ML/(2*MW*SW/EL),Conjugate[Yu[1,1]]:>Sqrt[2]*MU/(2*MW*SW/EL),Conjugate[Yu[2,2]]:>Sqrt[2]*MC/(2*MW*SW/EL),Conjugate[Yu[3,3]]:>Sqrt[2]*MT/(2*MW*SW/EL),Conjugate[Yd[1,1]]:>Sqrt[2]*MD/(2*MW*SW/EL),Conjugate[Yd[2,2]]:>Sqrt[2]*MS/(2*MW*SW/EL),Conjugate[Yd[3,3]]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{Ye[1,1]:>Sqrt[2]*ME/(2*MW*SW/EL),Ye[2,2]:>Sqrt[2]*MM/(2*MW*SW/EL),Ye[3,3]:>Sqrt[2]*ML/(2*MW*SW/EL),Yu[1,1]:>Sqrt[2]*MU/(2*MW*SW/EL),Yu[2,2]:>Sqrt[2]*MC/(2*MW*SW/EL),Yu[3,3]:>Sqrt[2]*MT/(2*MW*SW/EL),Yd[1,1]:>Sqrt[2]*MD/(2*MW*SW/EL),Yd[2,2]:>Sqrt[2]*MS/(2*MW*SW/EL),Yd[3,3]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{MassFe[1]:>ME,MassFe[2]:>MM,MassFe[3]:>ML,MassFv[_]:>0,MassFu[1]:>MU,MassFu[2]:>MC,MassFu[3]:>MT,MassFd[1]:>MD,MassFd[2]:>MS,MassFd[3]:>MB}/.{vd:>CB*(2*MW*SW/EL),vu:>SB*(2*MW*SW/EL)}/.{Lam1:>(-mu2*SB^2+MH1^2*RR[[1]][[1]]^2+MH2^2*RR[[2]][[1]]^2+MH3^2*RR[[3]][[1]]^2)/(CB^2*v^2),Lam2:>(-mu2*CB^2+MH1^2*RR[[1]][[2]]^2+MH2^2*RR[[2]][[2]]^2+MH3^2*RR[[3]][[2]]^2)/(SB^2*v^2),Lam3:>(-mu2+2*MHp^2+(MH1^2*RR[[1]][[1]]*RR[[1]][[2]]+MH2^2*RR[[2]][[1]]*RR[[2]][[2]]+MH3^2*RR[[3]][[1]]*RR[[3]][[2]])/SB/CB)/v^2,Lam4:>(mu2+MA0^2-2*MHp^2)/v^2,Lam5:>(mu2-MA0^2)/v^2,Lam6:>(MH1^2*RR[[1]][[3]]^2+MH2^2*RR[[2]][[3]]^2+MH3^2*RR[[3]][[3]]^2)/vS^2,Lam7:>(MH1^2*RR[[1]][[1]]*RR[[1]][[3]]+MH2^2*RR[[2]][[1]]*RR[[2]][[3]]+MH3^2*RR[[3]][[1]]*RR[[3]][[3]])/(v*vS*CB),Lam8:>(MH1^2*RR[[1]][[2]]*RR[[1]][[3]]+MH2^2*RR[[2]][[2]]*RR[[2]][[3]]+MH3^2*RR[[3]][[2]]*RR[[3]][[3]])/(v*vS*SB)}/.{Masshh[1]:>MH1,Masshh[2]:>MH2,Masshh[3]:>MH3,MassAh[2]:>MA0,MassHm[2]:>MHp}/.{ZH[aa_,bb_]:>RR[[aa]][[bb]]}/.{ZA[aa_,bb_]:>RA[[aa]][[bb]],ZP[aa_,bb_]:>RA[[aa]][[bb]]}/.{v:>(2*MW*SW/EL)}/.{mu2:>m12squared/SB/CB}/.{YukS1Lep[1]:>YukS1Lep1,YukS1Lep[2]:>YukS1Lep2,YukS1Lep[3]:>YukS1Lep3,YukS2Lep[1]:>YukS2Lep1,YukS2Lep[2]:>YukS2Lep2,YukS3Lep[1]:>YukS3Lep1,YukS3Lep[2]:>YukS3Lep2}/.{YukS1Quark[1]:>YukS1Quark1,YukS1Quark[2]:>YukS1Quark2,YukS1Quark[3]:>YukS1Quark3,YukS2Quark[1]:>YukS2Quark1,YukS2Quark[2]:>YukS2Quark2,YukS3Quark[1]:>YukS3Quark1,YukS3Quark[2]:>YukS3Quark2}/.{GaugeXi[U[2]]:>1,GaugeXi[U[3]]:>1,GaugeXi[U[4]]:>1}/.{CS1S1S1[gt1_,gt2_,gt3_]:>ToExpression["CS1S1S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]],CS1S3S3[gt1_,gt2_,gt3_]:>ToExpression["CS1S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]],CS2S2S1[gt1_,gt2_,gt3_]:>ToExpression["CS2S2S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]],CS2S2S2S2[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS2S2S2S2f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS2S2S1S1[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS2S2S1S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS2S2S3S3[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS2S2S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS1S1S1S1[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS1S1S1S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS1S1S3S3[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS1S1S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS3S3S3S3[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS3S3S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]]};
			resultLeftHanded = Coefficient[resultTemp,\!\(TraditionalForm\`Dot[DiracGamma[Momentum[p]], DiracGamma[7]]\),1];
			resultTemp = Coefficient[resultTemp,\!\(TraditionalForm\`Dot[DiracGamma[Momentum[p]], DiracGamma[7]]\),0];
			resultLeftHanded += Coefficient[resultTemp,\!\(TraditionalForm\`DiracGamma[Momentum[k]] . DiracGamma[7]\),1];
			resultTemp = Coefficient[resultTemp,\!\(TraditionalForm\`DiracGamma[Momentum[k]] . DiracGamma[7]\),0];
			resultRightHanded = Coefficient[resultTemp,Dot[DiracGamma[Momentum[p]],DiracGamma[6]],1];
			resultTemp = Coefficient[resultTemp,\!\(TraditionalForm\`Dot[DiracGamma[Momentum[p]], DiracGamma[6]]\),0];
			resultRightHanded += Coefficient[resultTemp,DiracGamma[Momentum[k]].DiracGamma[6],1];
			resultTemp = Coefficient[resultTemp,DiracGamma[Momentum[k]].DiracGamma[6],0];
			If[ToString[mass]==ToString[massSecond],
				resultScalar = (Simplify[resultTemp]/.{(DiracGamma[6]+DiracGamma[7]):>1})/mass;
			,
				resultScalar = (Simplify[resultTemp]/.{(m1_*DiracGamma[6]+m2_*DiracGamma[7]):>1});
			];
			Switch[n,
				1, exportFinal = resultLeftHanded/.{C0[0,x_,x_,0,0,b_]:>C0Mine[DBLE[0],DBLE[x],DBLE[x],DBLE[0],DBLE[0],DBLE[b]]}; exportToDeriv = resultLeftHanded;,
				2, exportFinal = resultRightHanded/.{C0[0,x_,x_,0,0,b_]:>C0Mine[DBLE[0],DBLE[x],DBLE[x],DBLE[0],DBLE[0],DBLE[b]]}; exportToDeriv = resultRightHanded;,
				3, exportFinal = resultScalar/.{C0[0,x_,x_,0,0,b_]:>C0Mine[DBLE[0],DBLE[x],DBLE[x],DBLE[0],DBLE[0],DBLE[b]]}; exportToDeriv = resultScalar;
			];
			(*exportFinal=exportFinal/.{ZEL[a_,a_]:>1,ZER[a_,a_]:>1}/.{ZEL[__]:>0,ZER[__]:>0}/.{ZUL[a_,a_]:>1,ZUR[a_,a_]:>1}/.{ZUL[__]:>0,ZUR[__]:>0}/.{Conjugate[Ye[1,1]]:>Sqrt[2]*ME/(2*MW*SW/EL),Conjugate[Ye[2,2]]:>Sqrt[2]*MM/(2*MW*SW/EL),Conjugate[Ye[3,3]]:>Sqrt[2]*ML/(2*MW*SW/EL),Conjugate[Yu[1,1]]:>Sqrt[2]*MU/(2*MW*SW/EL),Conjugate[Yu[2,2]]:>Sqrt[2]*MC/(2*MW*SW/EL),Conjugate[Yu[3,3]]:>Sqrt[2]*MT/(2*MW*SW/EL),Conjugate[Yd[1,1]]:>Sqrt[2]*MD/(2*MW*SW/EL),Conjugate[Yd[2,2]]:>Sqrt[2]*MS/(2*MW*SW/EL),Conjugate[Yd[3,3]]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{Ye[1,1]:>Sqrt[2]*ME/(2*MW*SW/EL),Ye[2,2]:>Sqrt[2]*MM/(2*MW*SW/EL),Ye[3,3]:>Sqrt[2]*ML/(2*MW*SW/EL),Yu[1,1]:>Sqrt[2]*MU/(2*MW*SW/EL),Yu[2,2]:>Sqrt[2]*MC/(2*MW*SW/EL),Yu[3,3]:>Sqrt[2]*MT/(2*MW*SW/EL),Yd[1,1]:>Sqrt[2]*MD/(2*MW*SW/EL),Yd[2,2]:>Sqrt[2]*MS/(2*MW*SW/EL),Yd[3,3]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{MassFe[1]:>ME,MassFe[2]:>MM,MassFe[3]:>ML,MassFv[_]:>0,MassFu[1]:>MU,MassFu[2]:>MC,MassFu[3]:>MT,MassFd[1]:>MD,MassFd[2]:>MS,MassFd[3]:>MB}/.{vd:>CB*(2*MW*SW/EL),vu:>SB*(2*MW*SW/EL)}/.{Lam1:>(-mu2*SB^2+MH1^2*RR[[1]][[1]]^2+MH2^2*RR[[2]][[1]]^2+MH3^2*RR[[3]][[1]]^2)/(CB^2*v^2),Lam2:>(-mu2*CB^2+MH1^2*RR[[1]][[2]]^2+MH2^2*RR[[2]][[2]]^2+MH3^2*RR[[3]][[2]]^2)/(SB^2*v^2),Lam3:>(-mu2+2*MHp^2+(MH1^2*RR[[1]][[1]]*RR[[1]][[2]]+MH2^2*RR[[2]][[1]]*RR[[2]][[2]]+MH3^2*RR[[3]][[1]]*RR[[3]][[2]])/SB/CB)/v^2,Lam4:>(mu2+MA0^2-2*MHp^2)/v^2,Lam5:>(mu2-MA0^2)/v^2,Lam6:>(MH1^2*RR[[1]][[3]]^2+MH2^2*RR[[2]][[3]]^2+MH3^2*RR[[3]][[3]]^2)/vS^2,Lam7:>(MH1^2*RR[[1]][[1]]*RR[[1]][[3]]+MH2^2*RR[[2]][[1]]*RR[[2]][[3]]+MH3^2*RR[[3]][[1]]*RR[[3]][[3]])/(v*vS*CB),Lam8:>(MH1^2*RR[[1]][[2]]*RR[[1]][[3]]+MH2^2*RR[[2]][[2]]*RR[[2]][[3]]+MH3^2*RR[[3]][[2]]*RR[[3]][[3]])/(v*vS*SB)}/.{Masshh[1]:>MH1,Masshh[2]:>MH2,Masshh[3]:>MH3,MassAh[2]:>MA0,MassHm[2]:>MHp}/.{ZH[aa_,bb_]:>RR[[aa]][[bb]]}/.{ZA[aa_,bb_]:>RA[[aa]][[bb]],ZP[aa_,bb_]:>RA[[aa]][[bb]]}/.{v:>(2*MW*SW/EL)}/.{mu2:>m12squared/SB/CB}/.{YukS1Lep[1]:>YukS1Lep1,YukS1Lep[2]:>YukS1Lep2,YukS1Lep[3]:>YukS1Lep3,YukS2Lep[1]:>YukS2Lep1,YukS2Lep[2]:>YukS2Lep2,YukS3Lep[1]:>YukS3Lep1,YukS3Lep[2]:>YukS3Lep2}/.{YukS1Quark[1]:>YukS1Quark1,YukS1Quark[2]:>YukS1Quark2,YukS1Quark[3]:>YukS1Quark3,YukS2Quark[1]:>YukS2Quark1,YukS2Quark[2]:>YukS2Quark2,YukS3Quark[1]:>YukS3Quark1,YukS3Quark[2]:>YukS3Quark2}/.{GaugeXi[U[2]]:>GaugeXiZ,GaugeXi[U[3]]:>GaugeXiW,GaugeXi[U[4]]:>GaugeXiW};*)
			exportFinal=exportFinal/.{ZEL[a_,a_]:>1,ZER[a_,a_]:>1}/.{ZEL[__]:>0,ZER[__]:>0}/.{ZUL[a_,a_]:>1,ZUR[a_,a_]:>1}/.{ZUL[__]:>0,ZUR[__]:>0}/.{Conjugate[Ye[1,1]]:>Sqrt[2]*ME/(2*MW*SW/EL),Conjugate[Ye[2,2]]:>Sqrt[2]*MM/(2*MW*SW/EL),Conjugate[Ye[3,3]]:>Sqrt[2]*ML/(2*MW*SW/EL),Conjugate[Yu[1,1]]:>Sqrt[2]*MU/(2*MW*SW/EL),Conjugate[Yu[2,2]]:>Sqrt[2]*MC/(2*MW*SW/EL),Conjugate[Yu[3,3]]:>Sqrt[2]*MT/(2*MW*SW/EL),Conjugate[Yd[1,1]]:>Sqrt[2]*MD/(2*MW*SW/EL),Conjugate[Yd[2,2]]:>Sqrt[2]*MS/(2*MW*SW/EL),Conjugate[Yd[3,3]]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{Ye[1,1]:>Sqrt[2]*ME/(2*MW*SW/EL),Ye[2,2]:>Sqrt[2]*MM/(2*MW*SW/EL),Ye[3,3]:>Sqrt[2]*ML/(2*MW*SW/EL),Yu[1,1]:>Sqrt[2]*MU/(2*MW*SW/EL),Yu[2,2]:>Sqrt[2]*MC/(2*MW*SW/EL),Yu[3,3]:>Sqrt[2]*MT/(2*MW*SW/EL),Yd[1,1]:>Sqrt[2]*MD/(2*MW*SW/EL),Yd[2,2]:>Sqrt[2]*MS/(2*MW*SW/EL),Yd[3,3]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{MassFe[1]:>ME,MassFe[2]:>MM,MassFe[3]:>ML,MassFv[_]:>0,MassFu[1]:>MU,MassFu[2]:>MC,MassFu[3]:>MT,MassFd[1]:>MD,MassFd[2]:>MS,MassFd[3]:>MB}/.{vd:>CB*(2*MW*SW/EL),vu:>SB*(2*MW*SW/EL)}/.{Lam1:>(-mu2*SB^2+MH1^2*RR[[1]][[1]]^2+MH2^2*RR[[2]][[1]]^2+MH3^2*RR[[3]][[1]]^2)/(CB^2*v^2),Lam2:>(-mu2*CB^2+MH1^2*RR[[1]][[2]]^2+MH2^2*RR[[2]][[2]]^2+MH3^2*RR[[3]][[2]]^2)/(SB^2*v^2),Lam3:>(-mu2+2*MHp^2+(MH1^2*RR[[1]][[1]]*RR[[1]][[2]]+MH2^2*RR[[2]][[1]]*RR[[2]][[2]]+MH3^2*RR[[3]][[1]]*RR[[3]][[2]])/SB/CB)/v^2,Lam4:>(mu2+MA0^2-2*MHp^2)/v^2,Lam5:>(mu2-MA0^2)/v^2,Lam6:>(MH1^2*RR[[1]][[3]]^2+MH2^2*RR[[2]][[3]]^2+MH3^2*RR[[3]][[3]]^2)/vS^2,Lam7:>(MH1^2*RR[[1]][[1]]*RR[[1]][[3]]+MH2^2*RR[[2]][[1]]*RR[[2]][[3]]+MH3^2*RR[[3]][[1]]*RR[[3]][[3]])/(v*vS*CB),Lam8:>(MH1^2*RR[[1]][[2]]*RR[[1]][[3]]+MH2^2*RR[[2]][[2]]*RR[[2]][[3]]+MH3^2*RR[[3]][[2]]*RR[[3]][[3]])/(v*vS*SB)}/.{Masshh[1]:>MH1,Masshh[2]:>MH2,Masshh[3]:>MH3,MassAh[2]:>MA0,MassHm[2]:>MHp}/.{ZH[aa_,bb_]:>RR[[aa]][[bb]]}/.{ZA[aa_,bb_]:>RA[[aa]][[bb]],ZP[aa_,bb_]:>RA[[aa]][[bb]]}/.{v:>(2*MW*SW/EL)}/.{mu2:>m12squared/SB/CB}/.{YukS1Lep[1]:>YukS1Lep1,YukS1Lep[2]:>YukS1Lep2,YukS1Lep[3]:>YukS1Lep3,YukS2Lep[1]:>YukS2Lep1,YukS2Lep[2]:>YukS2Lep2,YukS3Lep[1]:>YukS3Lep1,YukS3Lep[2]:>YukS3Lep2}/.{YukS1Quark[1]:>YukS1Quark1,YukS1Quark[2]:>YukS1Quark2,YukS1Quark[3]:>YukS1Quark3,YukS2Quark[1]:>YukS2Quark1,YukS2Quark[2]:>YukS2Quark2,YukS3Quark[1]:>YukS3Quark1,YukS3Quark[2]:>YukS3Quark2}/.{GaugeXi[U[2]]:>1,GaugeXi[U[3]]:>1,GaugeXi[U[4]]:>1}/.{CS1S1S1[gt1_,gt2_,gt3_]:>ToExpression["CS1S1S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]],CS1S3S3[gt1_,gt2_,gt3_]:>ToExpression["CS1S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]],CS2S2S1[gt1_,gt2_,gt3_]:>ToExpression["CS2S2S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]],CS2S2S2S2[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS2S2S2S2f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS2S2S1S1[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS2S2S1S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS2S2S3S3[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS2S2S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS1S1S1S1[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS1S1S1S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS1S1S3S3[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS1S1S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS3S3S3S3[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS3S3S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]]};
			exportFortran=exportFortran<>feyncalcToFortran[exportFinal," amplitudes("<>ToString[i]<>") = "]<>"\n\n";
			resultDeriv = D[exportToDeriv,x];
			exportDeriv = (resultDeriv/.{C0[0,x_,x_,0,0,b_]:>C0Mine[DBLE[0],DBLE[x],DBLE[x],DBLE[0],DBLE[0],DBLE[b]]}/.{Derivative[1,0,0][B00][x_,m1_,m2_]:>DB00[DBLE[x],DBLE[m1],DBLE[m2]],Derivative[1,0,0][B11][x_,m1_,m2_]:>DB11[DBLE[x],DBLE[m1],DBLE[m2]],Derivative[0,1,0,0,0,0][C0][0,x_,x_,0,0,m_]:>DC01Mine[DBLE[0],DBLE[x],DBLE[x],DBLE[0],DBLE[0],DBLE[m]],Derivative[0,0,1,0,0,0][C0][0,x_,x_,0,0,m_]:>DC02Mine[DBLE[0],DBLE[x],DBLE[x],DBLE[0],DBLE[0],DBLE[m]]});
			(*exportDeriv=exportDeriv/.{ZEL[a_,a_]:>1,ZER[a_,a_]:>1}/.{ZEL[__]:>0,ZER[__]:>0}/.{ZUL[a_,a_]:>1,ZUR[a_,a_]:>1}/.{ZUL[__]:>0,ZUR[__]:>0}/.{Conjugate[Ye[1,1]]:>Sqrt[2]*ME/(2*MW*SW/EL),Conjugate[Ye[2,2]]:>Sqrt[2]*MM/(2*MW*SW/EL),Conjugate[Ye[3,3]]:>Sqrt[2]*ML/(2*MW*SW/EL),Conjugate[Yu[1,1]]:>Sqrt[2]*MU/(2*MW*SW/EL),Conjugate[Yu[2,2]]:>Sqrt[2]*MC/(2*MW*SW/EL),Conjugate[Yu[3,3]]:>Sqrt[2]*MT/(2*MW*SW/EL),Conjugate[Yd[1,1]]:>Sqrt[2]*MD/(2*MW*SW/EL),Conjugate[Yd[2,2]]:>Sqrt[2]*MS/(2*MW*SW/EL),Conjugate[Yd[3,3]]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{Ye[1,1]:>Sqrt[2]*ME/(2*MW*SW/EL),Ye[2,2]:>Sqrt[2]*MM/(2*MW*SW/EL),Ye[3,3]:>Sqrt[2]*ML/(2*MW*SW/EL),Yu[1,1]:>Sqrt[2]*MU/(2*MW*SW/EL),Yu[2,2]:>Sqrt[2]*MC/(2*MW*SW/EL),Yu[3,3]:>Sqrt[2]*MT/(2*MW*SW/EL),Yd[1,1]:>Sqrt[2]*MD/(2*MW*SW/EL),Yd[2,2]:>Sqrt[2]*MS/(2*MW*SW/EL),Yd[3,3]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{MassFe[1]:>ME,MassFe[2]:>MM,MassFe[3]:>ML,MassFv[_]:>0,MassFu[1]:>MU,MassFu[2]:>MC,MassFu[3]:>MT,MassFd[1]:>MD,MassFd[2]:>MS,MassFd[3]:>MB}/.{vd:>CB*(2*MW*SW/EL),vu:>SB*(2*MW*SW/EL)}/.{Lam1:>(-mu2*SB^2+MH1^2*RR[[1]][[1]]^2+MH2^2*RR[[2]][[1]]^2+MH3^2*RR[[3]][[1]]^2)/(CB^2*v^2),Lam2:>(-mu2*CB^2+MH1^2*RR[[1]][[2]]^2+MH2^2*RR[[2]][[2]]^2+MH3^2*RR[[3]][[2]]^2)/(SB^2*v^2),Lam3:>(-mu2+2*MHp^2+(MH1^2*RR[[1]][[1]]*RR[[1]][[2]]+MH2^2*RR[[2]][[1]]*RR[[2]][[2]]+MH3^2*RR[[3]][[1]]*RR[[3]][[2]])/SB/CB)/v^2,Lam4:>(mu2+MA0^2-2*MHp^2)/v^2,Lam5:>(mu2-MA0^2)/v^2,Lam6:>(MH1^2*RR[[1]][[3]]^2+MH2^2*RR[[2]][[3]]^2+MH3^2*RR[[3]][[3]]^2)/vS^2,Lam7:>(MH1^2*RR[[1]][[1]]*RR[[1]][[3]]+MH2^2*RR[[2]][[1]]*RR[[2]][[3]]+MH3^2*RR[[3]][[1]]*RR[[3]][[3]])/(v*vS*CB),Lam8:>(MH1^2*RR[[1]][[2]]*RR[[1]][[3]]+MH2^2*RR[[2]][[2]]*RR[[2]][[3]]+MH3^2*RR[[3]][[2]]*RR[[3]][[3]])/(v*vS*SB)}/.{Masshh[1]:>MH1,Masshh[2]:>MH2,Masshh[3]:>MH3,MassAh[2]:>MA0,MassHm[2]:>MHp}/.{ZH[aa_,bb_]:>RR[[aa]][[bb]]}/.{ZA[aa_,bb_]:>RA[[aa]][[bb]],ZP[aa_,bb_]:>RA[[aa]][[bb]]}/.{v:>(2*MW*SW/EL)}/.{mu2:>m12squared/SB/CB}/.{YukS1Lep[1]:>YukS1Lep1,YukS1Lep[2]:>YukS1Lep2,YukS1Lep[3]:>YukS1Lep3,YukS2Lep[1]:>YukS2Lep1,YukS2Lep[2]:>YukS2Lep2,YukS3Lep[1]:>YukS3Lep1,YukS3Lep[2]:>YukS3Lep2}/.{YukS1Quark[1]:>YukS1Quark1,YukS1Quark[2]:>YukS1Quark2,YukS1Quark[3]:>YukS1Quark3,YukS2Quark[1]:>YukS2Quark1,YukS2Quark[2]:>YukS2Quark2,YukS3Quark[1]:>YukS3Quark1,YukS3Quark[2]:>YukS3Quark2}/.{GaugeXi[U[2]]:>GaugeXiZ,GaugeXi[U[3]]:>GaugeXiW,GaugeXi[U[4]]:>GaugeXiW};*)
			exportDeriv=exportDeriv/.{ZEL[a_,a_]:>1,ZER[a_,a_]:>1}/.{ZEL[__]:>0,ZER[__]:>0}/.{ZUL[a_,a_]:>1,ZUR[a_,a_]:>1}/.{ZUL[__]:>0,ZUR[__]:>0}/.{Conjugate[Ye[1,1]]:>Sqrt[2]*ME/(2*MW*SW/EL),Conjugate[Ye[2,2]]:>Sqrt[2]*MM/(2*MW*SW/EL),Conjugate[Ye[3,3]]:>Sqrt[2]*ML/(2*MW*SW/EL),Conjugate[Yu[1,1]]:>Sqrt[2]*MU/(2*MW*SW/EL),Conjugate[Yu[2,2]]:>Sqrt[2]*MC/(2*MW*SW/EL),Conjugate[Yu[3,3]]:>Sqrt[2]*MT/(2*MW*SW/EL),Conjugate[Yd[1,1]]:>Sqrt[2]*MD/(2*MW*SW/EL),Conjugate[Yd[2,2]]:>Sqrt[2]*MS/(2*MW*SW/EL),Conjugate[Yd[3,3]]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{Ye[1,1]:>Sqrt[2]*ME/(2*MW*SW/EL),Ye[2,2]:>Sqrt[2]*MM/(2*MW*SW/EL),Ye[3,3]:>Sqrt[2]*ML/(2*MW*SW/EL),Yu[1,1]:>Sqrt[2]*MU/(2*MW*SW/EL),Yu[2,2]:>Sqrt[2]*MC/(2*MW*SW/EL),Yu[3,3]:>Sqrt[2]*MT/(2*MW*SW/EL),Yd[1,1]:>Sqrt[2]*MD/(2*MW*SW/EL),Yd[2,2]:>Sqrt[2]*MS/(2*MW*SW/EL),Yd[3,3]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{MassFe[1]:>ME,MassFe[2]:>MM,MassFe[3]:>ML,MassFv[_]:>0,MassFu[1]:>MU,MassFu[2]:>MC,MassFu[3]:>MT,MassFd[1]:>MD,MassFd[2]:>MS,MassFd[3]:>MB}/.{vd:>CB*(2*MW*SW/EL),vu:>SB*(2*MW*SW/EL)}/.{Lam1:>(-mu2*SB^2+MH1^2*RR[[1]][[1]]^2+MH2^2*RR[[2]][[1]]^2+MH3^2*RR[[3]][[1]]^2)/(CB^2*v^2),Lam2:>(-mu2*CB^2+MH1^2*RR[[1]][[2]]^2+MH2^2*RR[[2]][[2]]^2+MH3^2*RR[[3]][[2]]^2)/(SB^2*v^2),Lam3:>(-mu2+2*MHp^2+(MH1^2*RR[[1]][[1]]*RR[[1]][[2]]+MH2^2*RR[[2]][[1]]*RR[[2]][[2]]+MH3^2*RR[[3]][[1]]*RR[[3]][[2]])/SB/CB)/v^2,Lam4:>(mu2+MA0^2-2*MHp^2)/v^2,Lam5:>(mu2-MA0^2)/v^2,Lam6:>(MH1^2*RR[[1]][[3]]^2+MH2^2*RR[[2]][[3]]^2+MH3^2*RR[[3]][[3]]^2)/vS^2,Lam7:>(MH1^2*RR[[1]][[1]]*RR[[1]][[3]]+MH2^2*RR[[2]][[1]]*RR[[2]][[3]]+MH3^2*RR[[3]][[1]]*RR[[3]][[3]])/(v*vS*CB),Lam8:>(MH1^2*RR[[1]][[2]]*RR[[1]][[3]]+MH2^2*RR[[2]][[2]]*RR[[2]][[3]]+MH3^2*RR[[3]][[2]]*RR[[3]][[3]])/(v*vS*SB)}/.{Masshh[1]:>MH1,Masshh[2]:>MH2,Masshh[3]:>MH3,MassAh[2]:>MA0,MassHm[2]:>MHp}/.{ZH[aa_,bb_]:>RR[[aa]][[bb]]}/.{ZA[aa_,bb_]:>RA[[aa]][[bb]],ZP[aa_,bb_]:>RA[[aa]][[bb]]}/.{v:>(2*MW*SW/EL)}/.{mu2:>m12squared/SB/CB}/.{YukS1Lep[1]:>YukS1Lep1,YukS1Lep[2]:>YukS1Lep2,YukS1Lep[3]:>YukS1Lep3,YukS2Lep[1]:>YukS2Lep1,YukS2Lep[2]:>YukS2Lep2,YukS3Lep[1]:>YukS3Lep1,YukS3Lep[2]:>YukS3Lep2}/.{YukS1Quark[1]:>YukS1Quark1,YukS1Quark[2]:>YukS1Quark2,YukS1Quark[3]:>YukS1Quark3,YukS2Quark[1]:>YukS2Quark1,YukS2Quark[2]:>YukS2Quark2,YukS3Quark[1]:>YukS3Quark1,YukS3Quark[2]:>YukS3Quark2}/.{GaugeXi[U[2]]:>1,GaugeXi[U[3]]:>1,GaugeXi[U[4]]:>1}/.{CS1S1S1[gt1_,gt2_,gt3_]:>ToExpression["CS1S1S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]],CS1S3S3[gt1_,gt2_,gt3_]:>ToExpression["CS1S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]],CS2S2S1[gt1_,gt2_,gt3_]:>ToExpression["CS2S2S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]],CS2S2S2S2[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS2S2S2S2f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS2S2S1S1[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS2S2S1S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS2S2S3S3[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS2S2S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS1S1S1S1[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS1S1S1S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS1S1S3S3[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS1S1S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS3S3S3S3[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS3S3S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]]};
			exportFortranDeriv=exportFortranDeriv<>feyncalcToFortran[exportDeriv," amplitudes("<>ToString[i]<>") = "]<>"\n\n";
			Print["Sub-amplitude number "<>ToString[ampCounter]<>" finished."];
			ClearAll[tracer, simpler, oneloop, result, resultDeriv, export, exportDeriv];
		,
			export = (result/.{C0[0,x_,x_,0,0,b_]:>C0Mine[DBLE[0],DBLE[x],DBLE[x],DBLE[0],DBLE[0],DBLE[b]]}); 
			(*export=export/.{ZEL[a_,a_]:>1,ZER[a_,a_]:>1}/.{ZEL[__]:>0,ZER[__]:>0}/.{ZUL[a_,a_]:>1,ZUR[a_,a_]:>1}/.{ZUL[__]:>0,ZUR[__]:>0}/.{Conjugate[Ye[1,1]]:>Sqrt[2]*ME/(2*MW*SW/EL),Conjugate[Ye[2,2]]:>Sqrt[2]*MM/(2*MW*SW/EL),Conjugate[Ye[3,3]]:>Sqrt[2]*ML/(2*MW*SW/EL),Conjugate[Yu[1,1]]:>Sqrt[2]*MU/(2*MW*SW/EL),Conjugate[Yu[2,2]]:>Sqrt[2]*MC/(2*MW*SW/EL),Conjugate[Yu[3,3]]:>Sqrt[2]*MT/(2*MW*SW/EL),Conjugate[Yd[1,1]]:>Sqrt[2]*MD/(2*MW*SW/EL),Conjugate[Yd[2,2]]:>Sqrt[2]*MS/(2*MW*SW/EL),Conjugate[Yd[3,3]]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{Ye[1,1]:>Sqrt[2]*ME/(2*MW*SW/EL),Ye[2,2]:>Sqrt[2]*MM/(2*MW*SW/EL),Ye[3,3]:>Sqrt[2]*ML/(2*MW*SW/EL),Yu[1,1]:>Sqrt[2]*MU/(2*MW*SW/EL),Yu[2,2]:>Sqrt[2]*MC/(2*MW*SW/EL),Yu[3,3]:>Sqrt[2]*MT/(2*MW*SW/EL),Yd[1,1]:>Sqrt[2]*MD/(2*MW*SW/EL),Yd[2,2]:>Sqrt[2]*MS/(2*MW*SW/EL),Yd[3,3]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{MassFe[1]:>ME,MassFe[2]:>MM,MassFe[3]:>ML,MassFv[_]:>0,MassFu[1]:>MU,MassFu[2]:>MC,MassFu[3]:>MT,MassFd[1]:>MD,MassFd[2]:>MS,MassFd[3]:>MB}/.{vd:>CB*(2*MW*SW/EL),vu:>SB*(2*MW*SW/EL)}/.{Lam1:>(-mu2*SB^2+MH1^2*RR[[1]][[1]]^2+MH2^2*RR[[2]][[1]]^2+MH3^2*RR[[3]][[1]]^2)/(CB^2*v^2),Lam2:>(-mu2*CB^2+MH1^2*RR[[1]][[2]]^2+MH2^2*RR[[2]][[2]]^2+MH3^2*RR[[3]][[2]]^2)/(SB^2*v^2),Lam3:>(-mu2+2*MHp^2+(MH1^2*RR[[1]][[1]]*RR[[1]][[2]]+MH2^2*RR[[2]][[1]]*RR[[2]][[2]]+MH3^2*RR[[3]][[1]]*RR[[3]][[2]])/SB/CB)/v^2,Lam4:>(mu2+MA0^2-2*MHp^2)/v^2,Lam5:>(mu2-MA0^2)/v^2,Lam6:>(MH1^2*RR[[1]][[3]]^2+MH2^2*RR[[2]][[3]]^2+MH3^2*RR[[3]][[3]]^2)/vS^2,Lam7:>(MH1^2*RR[[1]][[1]]*RR[[1]][[3]]+MH2^2*RR[[2]][[1]]*RR[[2]][[3]]+MH3^2*RR[[3]][[1]]*RR[[3]][[3]])/(v*vS*CB),Lam8:>(MH1^2*RR[[1]][[2]]*RR[[1]][[3]]+MH2^2*RR[[2]][[2]]*RR[[2]][[3]]+MH3^2*RR[[3]][[2]]*RR[[3]][[3]])/(v*vS*SB)}/.{Masshh[1]:>MH1,Masshh[2]:>MH2,Masshh[3]:>MH3,MassAh[2]:>MA0,MassHm[2]:>MHp}/.{ZH[aa_,bb_]:>RR[[aa]][[bb]]}/.{ZA[aa_,bb_]:>RA[[aa]][[bb]],ZP[aa_,bb_]:>RA[[aa]][[bb]]}/.{v:>(2*MW*SW/EL)}/.{mu2:>m12squared/SB/CB}/.{YukS1Lep[1]:>YukS1Lep1,YukS1Lep[2]:>YukS1Lep2,YukS1Lep[3]:>YukS1Lep3,YukS2Lep[1]:>YukS2Lep1,YukS2Lep[2]:>YukS2Lep2,YukS3Lep[1]:>YukS3Lep1,YukS3Lep[2]:>YukS3Lep2}/.{YukS1Quark[1]:>YukS1Quark1,YukS1Quark[2]:>YukS1Quark2,YukS1Quark[3]:>YukS1Quark3,YukS2Quark[1]:>YukS2Quark1,YukS2Quark[2]:>YukS2Quark2,YukS3Quark[1]:>YukS3Quark1,YukS3Quark[2]:>YukS3Quark2}/.{GaugeXi[U[2]]:>GaugeXiZ,GaugeXi[U[3]]:>GaugeXiW,GaugeXi[U[4]]:>GaugeXiW};*)
			export=export/.{ZEL[a_,a_]:>1,ZER[a_,a_]:>1}/.{ZEL[__]:>0,ZER[__]:>0}/.{ZUL[a_,a_]:>1,ZUR[a_,a_]:>1}/.{ZUL[__]:>0,ZUR[__]:>0}/.{Conjugate[Ye[1,1]]:>Sqrt[2]*ME/(2*MW*SW/EL),Conjugate[Ye[2,2]]:>Sqrt[2]*MM/(2*MW*SW/EL),Conjugate[Ye[3,3]]:>Sqrt[2]*ML/(2*MW*SW/EL),Conjugate[Yu[1,1]]:>Sqrt[2]*MU/(2*MW*SW/EL),Conjugate[Yu[2,2]]:>Sqrt[2]*MC/(2*MW*SW/EL),Conjugate[Yu[3,3]]:>Sqrt[2]*MT/(2*MW*SW/EL),Conjugate[Yd[1,1]]:>Sqrt[2]*MD/(2*MW*SW/EL),Conjugate[Yd[2,2]]:>Sqrt[2]*MS/(2*MW*SW/EL),Conjugate[Yd[3,3]]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{Ye[1,1]:>Sqrt[2]*ME/(2*MW*SW/EL),Ye[2,2]:>Sqrt[2]*MM/(2*MW*SW/EL),Ye[3,3]:>Sqrt[2]*ML/(2*MW*SW/EL),Yu[1,1]:>Sqrt[2]*MU/(2*MW*SW/EL),Yu[2,2]:>Sqrt[2]*MC/(2*MW*SW/EL),Yu[3,3]:>Sqrt[2]*MT/(2*MW*SW/EL),Yd[1,1]:>Sqrt[2]*MD/(2*MW*SW/EL),Yd[2,2]:>Sqrt[2]*MS/(2*MW*SW/EL),Yd[3,3]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{MassFe[1]:>ME,MassFe[2]:>MM,MassFe[3]:>ML,MassFv[_]:>0,MassFu[1]:>MU,MassFu[2]:>MC,MassFu[3]:>MT,MassFd[1]:>MD,MassFd[2]:>MS,MassFd[3]:>MB}/.{vd:>CB*(2*MW*SW/EL),vu:>SB*(2*MW*SW/EL)}/.{Lam1:>(-mu2*SB^2+MH1^2*RR[[1]][[1]]^2+MH2^2*RR[[2]][[1]]^2+MH3^2*RR[[3]][[1]]^2)/(CB^2*v^2),Lam2:>(-mu2*CB^2+MH1^2*RR[[1]][[2]]^2+MH2^2*RR[[2]][[2]]^2+MH3^2*RR[[3]][[2]]^2)/(SB^2*v^2),Lam3:>(-mu2+2*MHp^2+(MH1^2*RR[[1]][[1]]*RR[[1]][[2]]+MH2^2*RR[[2]][[1]]*RR[[2]][[2]]+MH3^2*RR[[3]][[1]]*RR[[3]][[2]])/SB/CB)/v^2,Lam4:>(mu2+MA0^2-2*MHp^2)/v^2,Lam5:>(mu2-MA0^2)/v^2,Lam6:>(MH1^2*RR[[1]][[3]]^2+MH2^2*RR[[2]][[3]]^2+MH3^2*RR[[3]][[3]]^2)/vS^2,Lam7:>(MH1^2*RR[[1]][[1]]*RR[[1]][[3]]+MH2^2*RR[[2]][[1]]*RR[[2]][[3]]+MH3^2*RR[[3]][[1]]*RR[[3]][[3]])/(v*vS*CB),Lam8:>(MH1^2*RR[[1]][[2]]*RR[[1]][[3]]+MH2^2*RR[[2]][[2]]*RR[[2]][[3]]+MH3^2*RR[[3]][[2]]*RR[[3]][[3]])/(v*vS*SB)}/.{Masshh[1]:>MH1,Masshh[2]:>MH2,Masshh[3]:>MH3,MassAh[2]:>MA0,MassHm[2]:>MHp}/.{ZH[aa_,bb_]:>RR[[aa]][[bb]]}/.{ZA[aa_,bb_]:>RA[[aa]][[bb]],ZP[aa_,bb_]:>RA[[aa]][[bb]]}/.{v:>(2*MW*SW/EL)}/.{mu2:>m12squared/SB/CB}/.{YukS1Lep[1]:>YukS1Lep1,YukS1Lep[2]:>YukS1Lep2,YukS1Lep[3]:>YukS1Lep3,YukS2Lep[1]:>YukS2Lep1,YukS2Lep[2]:>YukS2Lep2,YukS3Lep[1]:>YukS3Lep1,YukS3Lep[2]:>YukS3Lep2}/.{YukS1Quark[1]:>YukS1Quark1,YukS1Quark[2]:>YukS1Quark2,YukS1Quark[3]:>YukS1Quark3,YukS2Quark[1]:>YukS2Quark1,YukS2Quark[2]:>YukS2Quark2,YukS3Quark[1]:>YukS3Quark1,YukS3Quark[2]:>YukS3Quark2}/.{GaugeXi[U[2]]:>1,GaugeXi[U[3]]:>1,GaugeXi[U[4]]:>1}/.{CS1S1S1[gt1_,gt2_,gt3_]:>ToExpression["CS1S1S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]],CS1S3S3[gt1_,gt2_,gt3_]:>ToExpression["CS1S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]],CS2S2S1[gt1_,gt2_,gt3_]:>ToExpression["CS2S2S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]],CS2S2S2S2[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS2S2S2S2f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS2S2S1S1[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS2S2S1S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS2S2S3S3[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS2S2S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS1S1S1S1[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS1S1S1S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS1S1S3S3[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS1S1S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS3S3S3S3[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS3S3S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]]};
			exportFortran=exportFortran<>feyncalcToFortran[export," amplitudes("<>ToString[i]<>") = "]<>"\n\n";
			If[argv[[2]]=="1",
				resultDeriv = D[result,x];
				exportDeriv = (resultDeriv/.{C0[0,x_,x_,0,0,b_]:>C0Mine[DBLE[0],DBLE[x],DBLE[x],DBLE[0],DBLE[0],DBLE[b]]}/.{Derivative[1,0,0][B00][x_,m1_,m2_]:>DB00[DBLE[x],DBLE[m1],DBLE[m2]],Derivative[1,0,0][B11][x_,m1_,m2_]:>DB11[DBLE[x],DBLE[m1],DBLE[m2]],Derivative[0,1,0,0,0,0][C0][0,x_,x_,0,0,m_]:>DC01Mine[DBLE[0],DBLE[x],DBLE[x],DBLE[0],DBLE[0],DBLE[m]],Derivative[0,0,1,0,0,0][C0][0,x_,x_,0,0,m_]:>DC02Mine[DBLE[0],DBLE[x],DBLE[x],DBLE[0],DBLE[0],DBLE[m]]});
				(*exportDeriv=exportDeriv/.{ZEL[a_,a_]:>1,ZER[a_,a_]:>1}/.{ZEL[__]:>0,ZER[__]:>0}/.{ZUL[a_,a_]:>1,ZUR[a_,a_]:>1}/.{ZUL[__]:>0,ZUR[__]:>0}/.{Conjugate[Ye[1,1]]:>Sqrt[2]*ME/(2*MW*SW/EL),Conjugate[Ye[2,2]]:>Sqrt[2]*MM/(2*MW*SW/EL),Conjugate[Ye[3,3]]:>Sqrt[2]*ML/(2*MW*SW/EL),Conjugate[Yu[1,1]]:>Sqrt[2]*MU/(2*MW*SW/EL),Conjugate[Yu[2,2]]:>Sqrt[2]*MC/(2*MW*SW/EL),Conjugate[Yu[3,3]]:>Sqrt[2]*MT/(2*MW*SW/EL),Conjugate[Yd[1,1]]:>Sqrt[2]*MD/(2*MW*SW/EL),Conjugate[Yd[2,2]]:>Sqrt[2]*MS/(2*MW*SW/EL),Conjugate[Yd[3,3]]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{Ye[1,1]:>Sqrt[2]*ME/(2*MW*SW/EL),Ye[2,2]:>Sqrt[2]*MM/(2*MW*SW/EL),Ye[3,3]:>Sqrt[2]*ML/(2*MW*SW/EL),Yu[1,1]:>Sqrt[2]*MU/(2*MW*SW/EL),Yu[2,2]:>Sqrt[2]*MC/(2*MW*SW/EL),Yu[3,3]:>Sqrt[2]*MT/(2*MW*SW/EL),Yd[1,1]:>Sqrt[2]*MD/(2*MW*SW/EL),Yd[2,2]:>Sqrt[2]*MS/(2*MW*SW/EL),Yd[3,3]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{MassFe[1]:>ME,MassFe[2]:>MM,MassFe[3]:>ML,MassFv[_]:>0,MassFu[1]:>MU,MassFu[2]:>MC,MassFu[3]:>MT,MassFd[1]:>MD,MassFd[2]:>MS,MassFd[3]:>MB}/.{vd:>CB*(2*MW*SW/EL),vu:>SB*(2*MW*SW/EL)}/.{Lam1:>(-mu2*SB^2+MH1^2*RR[[1]][[1]]^2+MH2^2*RR[[2]][[1]]^2+MH3^2*RR[[3]][[1]]^2)/(CB^2*v^2),Lam2:>(-mu2*CB^2+MH1^2*RR[[1]][[2]]^2+MH2^2*RR[[2]][[2]]^2+MH3^2*RR[[3]][[2]]^2)/(SB^2*v^2),Lam3:>(-mu2+2*MHp^2+(MH1^2*RR[[1]][[1]]*RR[[1]][[2]]+MH2^2*RR[[2]][[1]]*RR[[2]][[2]]+MH3^2*RR[[3]][[1]]*RR[[3]][[2]])/SB/CB)/v^2,Lam4:>(mu2+MA0^2-2*MHp^2)/v^2,Lam5:>(mu2-MA0^2)/v^2,Lam6:>(MH1^2*RR[[1]][[3]]^2+MH2^2*RR[[2]][[3]]^2+MH3^2*RR[[3]][[3]]^2)/vS^2,Lam7:>(MH1^2*RR[[1]][[1]]*RR[[1]][[3]]+MH2^2*RR[[2]][[1]]*RR[[2]][[3]]+MH3^2*RR[[3]][[1]]*RR[[3]][[3]])/(v*vS*CB),Lam8:>(MH1^2*RR[[1]][[2]]*RR[[1]][[3]]+MH2^2*RR[[2]][[2]]*RR[[2]][[3]]+MH3^2*RR[[3]][[2]]*RR[[3]][[3]])/(v*vS*SB)}/.{Masshh[1]:>MH1,Masshh[2]:>MH2,Masshh[3]:>MH3,MassAh[2]:>MA0,MassHm[2]:>MHp}/.{ZH[aa_,bb_]:>RR[[aa]][[bb]]}/.{ZA[aa_,bb_]:>RA[[aa]][[bb]],ZP[aa_,bb_]:>RA[[aa]][[bb]]}/.{v:>(2*MW*SW/EL)}/.{mu2:>m12squared/SB/CB}/.{YukS1Lep[1]:>YukS1Lep1,YukS1Lep[2]:>YukS1Lep2,YukS1Lep[3]:>YukS1Lep3,YukS2Lep[1]:>YukS2Lep1,YukS2Lep[2]:>YukS2Lep2,YukS3Lep[1]:>YukS3Lep1,YukS3Lep[2]:>YukS3Lep2}/.{YukS1Quark[1]:>YukS1Quark1,YukS1Quark[2]:>YukS1Quark2,YukS1Quark[3]:>YukS1Quark3,YukS2Quark[1]:>YukS2Quark1,YukS2Quark[2]:>YukS2Quark2,YukS3Quark[1]:>YukS3Quark1,YukS3Quark[2]:>YukS3Quark2}/.{GaugeXi[U[2]]:>GaugeXiZ,GaugeXi[U[3]]:>GaugeXiW,GaugeXi[U[4]]:>GaugeXiW};*)
				exportDeriv=exportDeriv/.{ZEL[a_,a_]:>1,ZER[a_,a_]:>1}/.{ZEL[__]:>0,ZER[__]:>0}/.{ZUL[a_,a_]:>1,ZUR[a_,a_]:>1}/.{ZUL[__]:>0,ZUR[__]:>0}/.{Conjugate[Ye[1,1]]:>Sqrt[2]*ME/(2*MW*SW/EL),Conjugate[Ye[2,2]]:>Sqrt[2]*MM/(2*MW*SW/EL),Conjugate[Ye[3,3]]:>Sqrt[2]*ML/(2*MW*SW/EL),Conjugate[Yu[1,1]]:>Sqrt[2]*MU/(2*MW*SW/EL),Conjugate[Yu[2,2]]:>Sqrt[2]*MC/(2*MW*SW/EL),Conjugate[Yu[3,3]]:>Sqrt[2]*MT/(2*MW*SW/EL),Conjugate[Yd[1,1]]:>Sqrt[2]*MD/(2*MW*SW/EL),Conjugate[Yd[2,2]]:>Sqrt[2]*MS/(2*MW*SW/EL),Conjugate[Yd[3,3]]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{Ye[1,1]:>Sqrt[2]*ME/(2*MW*SW/EL),Ye[2,2]:>Sqrt[2]*MM/(2*MW*SW/EL),Ye[3,3]:>Sqrt[2]*ML/(2*MW*SW/EL),Yu[1,1]:>Sqrt[2]*MU/(2*MW*SW/EL),Yu[2,2]:>Sqrt[2]*MC/(2*MW*SW/EL),Yu[3,3]:>Sqrt[2]*MT/(2*MW*SW/EL),Yd[1,1]:>Sqrt[2]*MD/(2*MW*SW/EL),Yd[2,2]:>Sqrt[2]*MS/(2*MW*SW/EL),Yd[3,3]:>Sqrt[2]*MB/(2*MW*SW/EL)}/.{MassFe[1]:>ME,MassFe[2]:>MM,MassFe[3]:>ML,MassFv[_]:>0,MassFu[1]:>MU,MassFu[2]:>MC,MassFu[3]:>MT,MassFd[1]:>MD,MassFd[2]:>MS,MassFd[3]:>MB}/.{vd:>CB*(2*MW*SW/EL),vu:>SB*(2*MW*SW/EL)}/.{Lam1:>(-mu2*SB^2+MH1^2*RR[[1]][[1]]^2+MH2^2*RR[[2]][[1]]^2+MH3^2*RR[[3]][[1]]^2)/(CB^2*v^2),Lam2:>(-mu2*CB^2+MH1^2*RR[[1]][[2]]^2+MH2^2*RR[[2]][[2]]^2+MH3^2*RR[[3]][[2]]^2)/(SB^2*v^2),Lam3:>(-mu2+2*MHp^2+(MH1^2*RR[[1]][[1]]*RR[[1]][[2]]+MH2^2*RR[[2]][[1]]*RR[[2]][[2]]+MH3^2*RR[[3]][[1]]*RR[[3]][[2]])/SB/CB)/v^2,Lam4:>(mu2+MA0^2-2*MHp^2)/v^2,Lam5:>(mu2-MA0^2)/v^2,Lam6:>(MH1^2*RR[[1]][[3]]^2+MH2^2*RR[[2]][[3]]^2+MH3^2*RR[[3]][[3]]^2)/vS^2,Lam7:>(MH1^2*RR[[1]][[1]]*RR[[1]][[3]]+MH2^2*RR[[2]][[1]]*RR[[2]][[3]]+MH3^2*RR[[3]][[1]]*RR[[3]][[3]])/(v*vS*CB),Lam8:>(MH1^2*RR[[1]][[2]]*RR[[1]][[3]]+MH2^2*RR[[2]][[2]]*RR[[2]][[3]]+MH3^2*RR[[3]][[2]]*RR[[3]][[3]])/(v*vS*SB)}/.{Masshh[1]:>MH1,Masshh[2]:>MH2,Masshh[3]:>MH3,MassAh[2]:>MA0,MassHm[2]:>MHp}/.{ZH[aa_,bb_]:>RR[[aa]][[bb]]}/.{ZA[aa_,bb_]:>RA[[aa]][[bb]],ZP[aa_,bb_]:>RA[[aa]][[bb]]}/.{v:>(2*MW*SW/EL)}/.{mu2:>m12squared/SB/CB}/.{YukS1Lep[1]:>YukS1Lep1,YukS1Lep[2]:>YukS1Lep2,YukS1Lep[3]:>YukS1Lep3,YukS2Lep[1]:>YukS2Lep1,YukS2Lep[2]:>YukS2Lep2,YukS3Lep[1]:>YukS3Lep1,YukS3Lep[2]:>YukS3Lep2}/.{YukS1Quark[1]:>YukS1Quark1,YukS1Quark[2]:>YukS1Quark2,YukS1Quark[3]:>YukS1Quark3,YukS2Quark[1]:>YukS2Quark1,YukS2Quark[2]:>YukS2Quark2,YukS3Quark[1]:>YukS3Quark1,YukS3Quark[2]:>YukS3Quark2}/.{GaugeXi[U[2]]:>1,GaugeXi[U[3]]:>1,GaugeXi[U[4]]:>1}/.{CS1S1S1[gt1_,gt2_,gt3_]:>ToExpression["CS1S1S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]],CS1S3S3[gt1_,gt2_,gt3_]:>ToExpression["CS1S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]],CS2S2S1[gt1_,gt2_,gt3_]:>ToExpression["CS2S2S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]],CS2S2S2S2[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS2S2S2S2f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS2S2S1S1[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS2S2S1S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS2S2S3S3[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS2S2S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS1S1S1S1[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS1S1S1S1f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS1S1S3S3[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS1S1S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]],CS3S3S3S3[gt1_,gt2_,gt3_,gt4_]:>ToExpression["CS3S3S3S3f"<>ToString[gt1]<>ToString[gt2]<>ToString[gt3]<>ToString[gt4]]};
				exportFortranDeriv=exportFortranDeriv<>feyncalcToFortran[exportDeriv," amplitudes("<>ToString[i]<>") = "]<>"\n\n";
				,
				False;
			];
			Print["Sub-amplitude number "<>ToString[i]<>" finished."];   
			ClearAll[tracer, simpler, oneloop, result, resultDeriv, export, exportDeriv]; ];
		];
		exportFortran=exportFortran<>"  totalAmplitude = (0D0,0D0)\n";
		exportFortran=exportFortran<>" do j=1,"<>ToString[Length[ToExpression[selfies[[counter]]]]]<>"\n";
		exportFortran=exportFortran<>"  totalAmplitude = totalAmplitude + amplitudes(j)\n";
		exportFortran=exportFortran<>" end do\n";
		exportFortran=exportFortran<>" " <> selfies[[counter]] <> fermionParts[[n]] <> ampVersion <> zeroMomDesc <> " = totalAmplitude\n";
		exportFortran=exportFortran<>"end function " <> selfies[[counter]] <> fermionParts[[n]] <> ampVersion <> zeroMomDesc <> "\n\n";
		exportFortranDeriv=exportFortranDeriv<>"  totalAmplitude = (0D0,0D0)\n";
		exportFortranDeriv=exportFortranDeriv<>" do j=1,"<>ToString[Length[ToExpression[selfies[[counter]]]]]<>"\n";
		exportFortranDeriv=exportFortranDeriv<>"  totalAmplitude = totalAmplitude + amplitudes(j)\n";
		exportFortranDeriv=exportFortranDeriv<>" end do\n";
		exportFortranDeriv=exportFortranDeriv<>" D" <> selfies[[counter]] <> fermionParts[[n]] <> " = totalAmplitude\n";
		exportFortranDeriv=exportFortranDeriv<>"end function D" <> selfies[[counter]] <> fermionParts[[n]] <> "\n\n";
		If[argv[[1]]=="3",
			(*If[counterIsZeroMom,
				filename = StringReplace[saver[[counter]],".txt"->"ZeroMom.txt"];
			,
				filename = saver[[3*(counter-1)+n]];
			];*)
			filename = saver[[3*(counter-1)+n]];
		,
			If[counterIsZeroMom,
				filename = StringReplace[saver[[counter]],".txt"->"ZeroMom.txt"];
			,
				filename = saver[[counter]];
			];
		];
		Export[filename, exportFortran]; 
		If[FileExistsQ[StringReplace[filename,".txt"->".F90"]],DeleteFile[StringReplace[filename,".txt"->".F90"]];];
		RenameFile[filename,StringReplace[filename,".txt"->".F90"]];
		If[argv[[2]]=="1",
			If[argv[[1]]=="3",
				filenameDeriv = saverDeriv[[3*(counter-1)+n]];
			,
				filenameDeriv = saverDeriv[[counter]];
			];
			Export[filenameDeriv, exportFortranDeriv];
			If[FileExistsQ[StringReplace[filenameDeriv,".txt"->".F90"]],DeleteFile[StringReplace[filenameDeriv,".txt"->".F90"]];];
			RenameFile[filenameDeriv,StringReplace[filenameDeriv,".txt"->".F90"]];
		,
			False;
		];
		Print["Self-energy "<> selfies[[counter]] <> ampVersion <>" finished."];
		Clear[FAamp, simpler, result, export, exportDeriv, filename, filenameDeriv, exportFortran, exportFortranDeriv]; 
		If[counterIsZeroMom,
			If[argv[[1]]=="3",
				If[3*(counter-1)+n>=9,
					counterIsZeroMom=False;
				,
					False;
				];
			,
				False;
			];
		,
			If[argv[[1]]=="2",
				If[counter==4,counter--;counterIsZeroMom=True;,False;];,
				False;
			];
		];
	];
]
DeleteFile[ImportFile];
Exit[]; 
