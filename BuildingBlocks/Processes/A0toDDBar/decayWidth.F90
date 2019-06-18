program decayWidth
    use constants
    use counterterms
    implicit none
#include "looptools.h"
    character(len=20) :: tempVal
    character(len=32) :: arg
    character(len=50) :: fileName, fileNameFilled, targetName
    character(len=600000) :: outputFileContent
    character(300), parameter :: pathToOutputFiles = 'Temp\\Results\\'
    integer arguments(5)
    integer, parameter :: maxNumberSchemes = 10
    logical :: debugModeOn = .false.
    logical :: isIRDivergent = .false.
    logical :: isUVDivergent = .false.
    logical :: isGaugeDependent = .false.
    logical :: isGaugeDependentAList(maxNumberSchemes), isGaugeDependentWList(maxNumberSchemes), &
        & isGaugeDependentZList(maxNumberSchemes)
    logical :: isIRDivergentList(maxNumberSchemes), isUVDivergentList(maxNumberSchemes)
    character isIRDivergentContinue, isUVDivergentContinue, isGaugeDependentContinue
    double precision, parameter :: GaugeDependenceThreshold = 1D-5
    double precision, parameter :: IRDivergenceThreshold = 1D-11
    double precision, parameter :: UVDivergenceThreshold = 1D-3
    double precision prefactor, treeLevelWidth, NLOWidth(maxNumberSchemes), fullamplitude(maxNumberSchemes)
    double precision NLOVCwidth, NLOVCwoIRwidth, NLOIRonlywidth
    double precision GaugeXiAChecks(maxNumberSchemes), GaugeXiWChecks(maxNumberSchemes), GaugeXiZChecks(maxNumberSchemes), &
        & IRChecks(maxNumberSchemes), UVChecks(maxNumberSchemes)
    double precision A0toDDBarTree, A0toDDBarCT, A0toDDBarReal, treeLevelTemp, realCorrectionsTemp
    double complex A0toDDBarVC, A0toDDBarTad, vertexCorrectionsTemp, vertexTadpolesTemp
    integer m, n, o, p, q, r, fileNameLength, point, statWrite

    ! Get the command line arguments standing for the different running options
    ! Argument 1: perform UV divergence check (1: true, 0: false; default: 0)
    ! Argument 2: perform IR divergence check (1: true, 0: false; default: 0)
    ! Argument 3: not implemented
    ! Argument 4: perform numerical evaluation (1: true, 0: false; default: 1)
    ! Argument 5: relative path to the 2HDM input parameter file, starting from the Parameters directory of 2HDMCalc
    ! Argument 6: relative path to the target file containing the results of the calculation, starting from the Temp/Results directory of 2HDMCalc
    do o = 1, iargc()
        call getarg(o, arg)
        if (arg == '1') then
            arguments(o) = 1
        else if (arg == '2') then
            arguments(o) = 2
        else if (arg == '0') then
            arguments(o) = 0
        else
            if (o == 5) then
                fileName = arg
            else if (o == 6) then
                targetName = arg
            end if
        end if
    end do

    if (arguments(3) == 1) then
        print *, "ERROR: Bad parameter. Ending N2HDECAY now."
        stop
    end if

    ! Perform the check for UV divergence
    if (arguments(1) == 1) then
        print *, "Checking for UV divergences ..."

        ! Use this hack to "fill up" the string to the maximum length with whitespace characters so that it can be passed to the subroutine call
        fileName = fileName // ' '

        ! Get all parameters
        call getParameters(fileName)

        ! Set the 2HDM parameters according to the first point in phase-space (random choice)
        MH2 = MH2List(1)
        MH3 = MH3List(1)
        MA0 = MA0List(1)
        MHp = MHpList(1)
        alpha1 = alpha1List(1)
        alpha2 = alpha2List(1)
        alpha3 = alpha3List(1)
        beta = betaList(1)
        CA1 = CA1List(1)
        CA2 = CA2List(1)
        CA3 = CA3List(1)
        SA1 = SA1List(1)
        SA2 = SA2List(1)
        SA3 = SA3List(1)
        CB = CBList(1)
        SB = SBList(1)
        TB = TBList(1)
        YukS1Lep1 = YukS1Lep1List(1)
        YukS1Lep2 = YukS1Lep2List(1)
        YukS1Lep3 = YukS1Lep3List(1)
        YukS2Lep1 = YukS2Lep1List(1)
        YukS2Lep2 = YukS2Lep2List(1)
        YukS3Lep1 = YukS3Lep1List(1)
        YukS3Lep2 = YukS3Lep2List(1)
        YukS1Quark1 = YukS1Quark1List(1)
        YukS1Quark2 = YukS1Quark2List(1)
        YukS1Quark3 = YukS1Quark3List(1)
        YukS2Quark1 = YukS2Quark1List(1)
        YukS2Quark2 = YukS2Quark2List(1)
        YukS3Quark1 = YukS3Quark1List(1)
        YukS3Quark2 = YukS3Quark2List(1)
        m12squared = m12squaredList(1)
        vS = vSList(1)

        ! Calculate the square of the 2HDM input parameters
        MH22 = MH2**2
        MH32 = MH3**2
        MA02 = MA0**2
        MHp2 = MHp**2
        CA12 = CA1**2
        CA22 = CA2**2
        CA32 = CA3**2
        SA12 = SA1**2
        SA22 = SA2**2
        SA32 = SA3**2
        TB2 = TB**2
        SB2 = SB**2
        CB2 = CB**2

        ! Set the scalar couplings
RR11 = CA1*CA2
RR12 = CA2*SA1
RR13 = SA2
RR21 = -1.D0*CA3*SA1 - 1.D0*CA1*SA2*SA3
RR22 = CA1*CA3 - 1.D0*SA1*SA2*SA3
RR23 = CA2*SA3
RR31 = -1.D0*CA1*CA3*SA2 + SA1*SA3
RR32 = -1.D0*CA3*SA1*SA2 - 1.D0*CA1*SA3
RR33 = CA2*CA3

Lam1 = (0.25D0*EL2*((-1.D0*m12squared*SB)/CB + MH12*DBLE(RR11**INT(2.D0)) + MH22*DBLE(RR21**INT(2.D0)) + MH32*DBLE(RR31**INT(2.D0&
  &))))/(CB2*MW2*SW2)
Lam2 = (0.25D0*EL2*((-1.D0*CB*m12squared)/SB + MH12*DBLE(RR12**INT(2.D0)) + MH22*DBLE(RR22**INT(2.D0)) + MH32*DBLE(RR32**INT(2.D0&
  &))))/(MW2*SB2*SW2)
Lam3 = (0.25D0*EL2*(2.D0*MHp2 - (1.D0*m12squared)/(CB*SB) + (MH12*RR11*RR12 + MH22*RR21*RR22 + MH32*RR31*RR32)/(CB*SB)))/(MW2*SW2&
  &)
Lam4 = (0.25D0*EL2*(MA02 - 2.D0*MHp2 + m12squared/(CB*SB)))/(MW2*SW2)
Lam5 = (0.25D0*EL2*(-1.D0*MA02 + m12squared/(CB*SB)))/(MW2*SW2)
Lam6 = (MH12*DBLE(RR13**INT(2.D0)) + MH22*DBLE(RR23**INT(2.D0)) + MH32*DBLE(RR33**INT(2.D0)))*DBLE(vS**INT(-2.D0))
Lam7 = (0.5D0*EL*(MH12*RR11*RR13 + MH22*RR21*RR23 + MH32*RR31*RR33))/(CB*MW*SW*vS)
Lam8 = (0.5D0*EL*(MH12*RR12*RR13 + MH22*RR22*RR23 + MH32*RR32*RR33))/(MW*SB*SW*vS)

CS1S1S1f111 = -1.D0*RR13*(Lam7*RR11*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR12*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR11*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f112 = -1.D0*RR13*(Lam7*RR11*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR12*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR11*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f113 = -1.D0*RR13*(Lam7*RR11*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR12*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR11*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f121 = -1.D0*RR13*(Lam7*RR21*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR22*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR21*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f122 = -1.D0*RR13*(Lam7*RR21*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR22*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR21*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f123 = -1.D0*RR13*(Lam7*RR21*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR22*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR21*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f131 = -1.D0*RR13*(Lam7*RR31*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR32*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR31*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f132 = -1.D0*RR13*(Lam7*RR31*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR32*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR31*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f133 = -1.D0*RR13*(Lam7*RR31*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR32*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR31*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f211 = -1.D0*RR23*(Lam7*RR11*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR12*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR11*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f212 = -1.D0*RR23*(Lam7*RR11*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR12*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR11*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f213 = -1.D0*RR23*(Lam7*RR11*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR12*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR11*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f221 = -1.D0*RR23*(Lam7*RR21*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR22*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR21*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f222 = -1.D0*RR23*(Lam7*RR21*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR22*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR21*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f223 = -1.D0*RR23*(Lam7*RR21*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR22*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR21*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f231 = -1.D0*RR23*(Lam7*RR31*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR32*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR31*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f232 = -1.D0*RR23*(Lam7*RR31*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR32*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR31*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f233 = -1.D0*RR23*(Lam7*RR31*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR32*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR31*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f311 = -1.D0*RR33*(Lam7*RR11*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR12*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR11*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f312 = -1.D0*RR33*(Lam7*RR11*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR12*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR11*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f313 = -1.D0*RR33*(Lam7*RR11*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR12*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR11*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f321 = -1.D0*RR33*(Lam7*RR21*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR22*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR21*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f322 = -1.D0*RR33*(Lam7*RR21*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR22*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR21*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f323 = -1.D0*RR33*(Lam7*RR21*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR22*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR21*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f331 = -1.D0*RR33*(Lam7*RR31*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR32*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR31*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f332 = -1.D0*RR33*(Lam7*RR31*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR32*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR31*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f333 = -1.D0*RR33*(Lam7*RR31*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR32*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR31*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))

CS2S2S1f111 = -1.D0*CB*(Lam5*SB*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + CB*((2.D0*CB*Lam1*MW*RR11*SW)/EL + (2.D0*(L&
  &am3 + Lam4 - 1.D0*Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*SB*(CB*Lam5*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/&
  &EL) + SB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR11*SW)/EL + (2.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS2S2S1f112 = -1.D0*CB*(Lam5*SB*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + CB*((2.D0*CB*Lam1*MW*RR21*SW)/EL + (2.D0*(L&
  &am3 + Lam4 - 1.D0*Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*SB*(CB*Lam5*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/&
  &EL) + SB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR21*SW)/EL + (2.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS2S2S1f113 = -1.D0*CB*(Lam5*SB*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + CB*((2.D0*CB*Lam1*MW*RR31*SW)/EL + (2.D0*(L&
  &am3 + Lam4 - 1.D0*Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*SB*(CB*Lam5*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/&
  &EL) + SB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR31*SW)/EL + (2.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS2S2S1f121 = -1.D0*CB*(CB*Lam5*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) - 1.D0*SB*((2.D0*CB*Lam1*MW*RR11*SW)/EL + (2.&
  &D0*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*SB*(-1.D0*Lam5*SB*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*R&
  &R11*SB*SW)/EL) + CB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR11*SW)/EL + (2.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS2S2S1f122 = -1.D0*CB*(CB*Lam5*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) - 1.D0*SB*((2.D0*CB*Lam1*MW*RR21*SW)/EL + (2.&
  &D0*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*SB*(-1.D0*Lam5*SB*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*R&
  &R21*SB*SW)/EL) + CB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR21*SW)/EL + (2.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS2S2S1f123 = -1.D0*CB*(CB*Lam5*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) - 1.D0*SB*((2.D0*CB*Lam1*MW*RR31*SW)/EL + (2.&
  &D0*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*SB*(-1.D0*Lam5*SB*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*R&
  &R31*SB*SW)/EL) + CB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR31*SW)/EL + (2.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS2S2S1f211 = SB*(Lam5*SB*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + CB*((2.D0*CB*Lam1*MW*RR11*SW)/EL + (2.D0*(Lam3 + &
  &Lam4 - 1.D0*Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*CB*(CB*Lam5*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + &
  &SB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR11*SW)/EL + (2.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS2S2S1f212 = SB*(Lam5*SB*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + CB*((2.D0*CB*Lam1*MW*RR21*SW)/EL + (2.D0*(Lam3 + &
  &Lam4 - 1.D0*Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*CB*(CB*Lam5*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + &
  &SB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR21*SW)/EL + (2.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS2S2S1f213 = SB*(Lam5*SB*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + CB*((2.D0*CB*Lam1*MW*RR31*SW)/EL + (2.D0*(Lam3 + &
  &Lam4 - 1.D0*Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*CB*(CB*Lam5*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + &
  &SB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR31*SW)/EL + (2.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS2S2S1f221 = SB*(CB*Lam5*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) - 1.D0*SB*((2.D0*CB*Lam1*MW*RR11*SW)/EL + (2.D0*(La&
  &m3 + Lam4 - 1.D0*Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*CB*(-1.D0*Lam5*SB*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB&
  &*SW)/EL) + CB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR11*SW)/EL + (2.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS2S2S1f222 = SB*(CB*Lam5*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) - 1.D0*SB*((2.D0*CB*Lam1*MW*RR21*SW)/EL + (2.D0*(La&
  &m3 + Lam4 - 1.D0*Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*CB*(-1.D0*Lam5*SB*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB&
  &*SW)/EL) + CB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR21*SW)/EL + (2.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS2S2S1f223 = SB*(CB*Lam5*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) - 1.D0*SB*((2.D0*CB*Lam1*MW*RR31*SW)/EL + (2.D0*(La&
  &m3 + Lam4 - 1.D0*Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*CB*(-1.D0*Lam5*SB*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB&
  &*SW)/EL) + CB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR31*SW)/EL + (2.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))

CS1S3S3f111 = 0.5D0*(-1.D0*RR12*(CB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + SB*((2.D0*CB2*(Lam4 + La&
  &m5)*MW*SW)/EL + (4.D0*Lam2*MW*SB2*SW)/EL)) - 1.D0*RR11*(SB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) &
  &+ CB*((4.D0*CB2*Lam1*MW*SW)/EL + (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR13*(CB2*Lam7 + Lam8*SB2)*vs)
CS1S3S3f112 = 0.5D0*(-1.D0*RR12*(SB*((4.D0*CB*Lam2*MW*SB*SW)/EL - (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + CB*((2.D0*CB2*(Lam4 + La&
  &m5)*MW*SW)/EL - (4.D0*Lam3*MW*SB2*SW)/EL)) - 1.D0*RR11*(CB*((-4.D0*CB*Lam1*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL)&
  & + SB*((4.D0*CB2*Lam3*MW*SW)/EL - (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR13*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*vs)
CS1S3S3f121 = 0.5D0*(-1.D0*RR12*(-1.D0*SB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + CB*((2.D0*CB2*(Lam&
  &4 + Lam5)*MW*SW)/EL + (4.D0*Lam2*MW*SB2*SW)/EL)) - 1.D0*RR11*(CB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW&
  &)/EL) - 1.D0*SB*((4.D0*CB2*Lam1*MW*SW)/EL + (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR13*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*vs&
  &)
CS1S3S3f122 = 0.5D0*(-1.D0*RR12*(CB*((4.D0*CB*Lam2*MW*SB*SW)/EL - (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) - 1.D0*SB*((2.D0*CB2*(Lam4&
  & + Lam5)*MW*SW)/EL - (4.D0*Lam3*MW*SB2*SW)/EL)) - 1.D0*RR11*(-1.D0*SB*((-4.D0*CB*Lam1*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW&
  &*SB*SW)/EL) + CB*((4.D0*CB2*Lam3*MW*SW)/EL - (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR13*(CB2*Lam8 + Lam7*SB2)*vs)
CS1S3S3f211 = 0.5D0*(-1.D0*RR22*(CB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + SB*((2.D0*CB2*(Lam4 + La&
  &m5)*MW*SW)/EL + (4.D0*Lam2*MW*SB2*SW)/EL)) - 1.D0*RR21*(SB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) &
  &+ CB*((4.D0*CB2*Lam1*MW*SW)/EL + (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR23*(CB2*Lam7 + Lam8*SB2)*vs)
CS1S3S3f212 = 0.5D0*(-1.D0*RR22*(SB*((4.D0*CB*Lam2*MW*SB*SW)/EL - (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + CB*((2.D0*CB2*(Lam4 + La&
  &m5)*MW*SW)/EL - (4.D0*Lam3*MW*SB2*SW)/EL)) - 1.D0*RR21*(CB*((-4.D0*CB*Lam1*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL)&
  & + SB*((4.D0*CB2*Lam3*MW*SW)/EL - (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR23*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*vs)
CS1S3S3f221 = 0.5D0*(-1.D0*RR22*(-1.D0*SB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + CB*((2.D0*CB2*(Lam&
  &4 + Lam5)*MW*SW)/EL + (4.D0*Lam2*MW*SB2*SW)/EL)) - 1.D0*RR21*(CB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW&
  &)/EL) - 1.D0*SB*((4.D0*CB2*Lam1*MW*SW)/EL + (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR23*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*vs&
  &)
CS1S3S3f222 = 0.5D0*(-1.D0*RR22*(CB*((4.D0*CB*Lam2*MW*SB*SW)/EL - (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) - 1.D0*SB*((2.D0*CB2*(Lam4&
  & + Lam5)*MW*SW)/EL - (4.D0*Lam3*MW*SB2*SW)/EL)) - 1.D0*RR21*(-1.D0*SB*((-4.D0*CB*Lam1*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW&
  &*SB*SW)/EL) + CB*((4.D0*CB2*Lam3*MW*SW)/EL - (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR23*(CB2*Lam8 + Lam7*SB2)*vs)
CS1S3S3f311 = 0.5D0*(-1.D0*RR32*(CB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + SB*((2.D0*CB2*(Lam4 + La&
  &m5)*MW*SW)/EL + (4.D0*Lam2*MW*SB2*SW)/EL)) - 1.D0*RR31*(SB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) &
  &+ CB*((4.D0*CB2*Lam1*MW*SW)/EL + (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR33*(CB2*Lam7 + Lam8*SB2)*vs)
CS1S3S3f312 = 0.5D0*(-1.D0*RR32*(SB*((4.D0*CB*Lam2*MW*SB*SW)/EL - (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + CB*((2.D0*CB2*(Lam4 + La&
  &m5)*MW*SW)/EL - (4.D0*Lam3*MW*SB2*SW)/EL)) - 1.D0*RR31*(CB*((-4.D0*CB*Lam1*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL)&
  & + SB*((4.D0*CB2*Lam3*MW*SW)/EL - (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*vs)
CS1S3S3f321 = 0.5D0*(-1.D0*RR32*(-1.D0*SB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + CB*((2.D0*CB2*(Lam&
  &4 + Lam5)*MW*SW)/EL + (4.D0*Lam2*MW*SB2*SW)/EL)) - 1.D0*RR31*(CB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW&
  &)/EL) - 1.D0*SB*((4.D0*CB2*Lam1*MW*SW)/EL + (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*vs&
  &)
CS1S3S3f322 = 0.5D0*(-1.D0*RR32*(CB*((4.D0*CB*Lam2*MW*SB*SW)/EL - (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) - 1.D0*SB*((2.D0*CB2*(Lam4&
  & + Lam5)*MW*SW)/EL - (4.D0*Lam3*MW*SB2*SW)/EL)) - 1.D0*RR31*(-1.D0*SB*((-4.D0*CB*Lam1*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW&
  &*SB*SW)/EL) + CB*((4.D0*CB2*Lam3*MW*SW)/EL - (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR33*(CB2*Lam8 + Lam7*SB2)*vs)

CS1S1S1S1f1111 = -1.D0*RR13*(2.D0*Lam7*RR13*DBLE(RR11**INT(2.D0)) + 2.D0*Lam8*RR13*DBLE(RR12**INT(2.D0)) + RR13*(Lam7*DBLE(RR11**&
  &INT(2.D0)) + Lam8*DBLE(RR12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*DBLE(R&
  &R12**INT(2.D0)) + 2.D0*Lam7*RR11*DBLE(RR13**INT(2.D0)) + RR11*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR1&
  &2**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*DBLE(RR11**INT(2.D0)) + 2.D0*Lam8*RR1&
  &2*DBLE(RR13**INT(2.D0)) + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13*&
  &*INT(2.D0))))
CS1S1S1S1f1112 = -1.D0*RR13*(Lam7*RR11*(RR13*RR21 + RR11*RR23) + Lam8*RR12*(RR13*RR22 + RR12*RR23) + RR13*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR21 + RR11*RR22) + Lam7*RR13*(RR13*RR21 + RR11&
  &*RR23) + RR11*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR13*(RR13*RR22 + RR12*RR23) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f1113 = -1.D0*RR13*(Lam7*RR11*(RR13*RR31 + RR11*RR33) + Lam8*RR12*(RR13*RR32 + RR12*RR33) + RR13*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR31 + RR11*RR32) + Lam7*RR13*(RR13*RR31 + RR11&
  &*RR33) + RR11*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR13*(RR13*RR32 + RR12*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f1121 = -1.D0*RR13*(Lam7*RR11*(RR13*RR21 + RR11*RR23) + Lam8*RR12*(RR13*RR22 + RR12*RR23) + RR13*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR21 + RR11*RR22) + Lam7*RR13*(RR13*RR21 + RR11&
  &*RR23) + RR11*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR13*(RR13*RR22 + RR12*RR23) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f1122 = -1.D0*RR13*(2.D0*Lam7*RR11*RR21*RR23 + 2.D0*Lam8*RR12*RR22*RR23 + RR13*(Lam7*DBLE(RR21**INT(2.D0)) + Lam8*DBLE(R&
  &R22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*RR21*RR22 + 2.D0*Lam7*RR13*RR2&
  &1*RR23 + RR11*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR22**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1&
  &.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR21*RR22 + 2.D0*Lam8*RR13*RR22*RR23 + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS1S1S1S1f1123 = -1.D0*RR13*(Lam7*RR11*(RR23*RR31 + RR21*RR33) + Lam8*RR12*(RR23*RR32 + RR22*RR33) + RR13*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*(RR22*RR31 + RR21*RR32) + Lam7*RR13*(RR23*RR31 + RR21&
  &*RR33) + RR11*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR13*(RR23*RR32 + RR22*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f1131 = -1.D0*RR13*(Lam7*RR11*(RR13*RR31 + RR11*RR33) + Lam8*RR12*(RR13*RR32 + RR12*RR33) + RR13*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR31 + RR11*RR32) + Lam7*RR13*(RR13*RR31 + RR11&
  &*RR33) + RR11*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR13*(RR13*RR32 + RR12*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f1132 = -1.D0*RR13*(Lam7*RR11*(RR23*RR31 + RR21*RR33) + Lam8*RR12*(RR23*RR32 + RR22*RR33) + RR13*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*(RR22*RR31 + RR21*RR32) + Lam7*RR13*(RR23*RR31 + RR21&
  &*RR33) + RR11*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR13*(RR23*RR32 + RR22*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f1133 = -1.D0*RR13*(2.D0*Lam7*RR11*RR31*RR33 + 2.D0*Lam8*RR12*RR32*RR33 + RR13*(Lam7*DBLE(RR31**INT(2.D0)) + Lam8*DBLE(R&
  &R32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*RR31*RR32 + 2.D0*Lam7*RR13*RR3&
  &1*RR33 + RR11*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR32**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1&
  &.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR31*RR32 + 2.D0*Lam8*RR13*RR32*RR33 + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS1S1S1S1f1211 = -1.D0*RR13*(2.D0*Lam7*RR11*RR13*RR21 + 2.D0*Lam8*RR12*RR13*RR22 + RR23*(Lam7*DBLE(RR11**INT(2.D0)) + Lam8*DBLE(R&
  &R12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR22 + 2.D0*Lam7*RR11*RR1&
  &3*RR23 + RR21*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR12**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1&
  &.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR21 + 2.D0*Lam8*RR12*RR13*RR23 + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS1S1S1S1f1212 = -1.D0*RR13*(Lam7*RR21*(RR13*RR21 + RR11*RR23) + Lam8*RR22*(RR13*RR22 + RR12*RR23) + RR23*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR21 + RR11*RR22) + Lam7*RR23*(RR13*RR21 + RR11&
  &*RR23) + RR21*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR23*(RR13*RR22 + RR12*RR23) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f1213 = -1.D0*RR13*(Lam7*RR21*(RR13*RR31 + RR11*RR33) + Lam8*RR22*(RR13*RR32 + RR12*RR33) + RR23*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR31 + RR11*RR32) + Lam7*RR23*(RR13*RR31 + RR11&
  &*RR33) + RR21*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR23*(RR13*RR32 + RR12*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f1221 = -1.D0*RR13*(Lam7*RR21*(RR13*RR21 + RR11*RR23) + Lam8*RR22*(RR13*RR22 + RR12*RR23) + RR23*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR21 + RR11*RR22) + Lam7*RR23*(RR13*RR21 + RR11&
  &*RR23) + RR21*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR23*(RR13*RR22 + RR12*RR23) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f1222 = -1.D0*RR13*(2.D0*Lam7*RR23*DBLE(RR21**INT(2.D0)) + 2.D0*Lam8*RR23*DBLE(RR22**INT(2.D0)) + RR23*(Lam7*DBLE(RR21**&
  &INT(2.D0)) + Lam8*DBLE(RR22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*DBLE(R&
  &R22**INT(2.D0)) + 2.D0*Lam7*RR21*DBLE(RR23**INT(2.D0)) + RR21*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR2&
  &2**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR22*DBLE(RR21**INT(2.D0)) + 2.D0*Lam8*RR2&
  &2*DBLE(RR23**INT(2.D0)) + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23*&
  &*INT(2.D0))))
CS1S1S1S1f1223 = -1.D0*RR13*(Lam7*RR21*(RR23*RR31 + RR21*RR33) + Lam8*RR22*(RR23*RR32 + RR22*RR33) + RR23*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*(RR22*RR31 + RR21*RR32) + Lam7*RR23*(RR23*RR31 + RR21&
  &*RR33) + RR21*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR23*(RR23*RR32 + RR22*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f1231 = -1.D0*RR13*(Lam7*RR21*(RR13*RR31 + RR11*RR33) + Lam8*RR22*(RR13*RR32 + RR12*RR33) + RR23*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR31 + RR11*RR32) + Lam7*RR23*(RR13*RR31 + RR11&
  &*RR33) + RR21*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR23*(RR13*RR32 + RR12*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f1232 = -1.D0*RR13*(Lam7*RR21*(RR23*RR31 + RR21*RR33) + Lam8*RR22*(RR23*RR32 + RR22*RR33) + RR23*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*(RR22*RR31 + RR21*RR32) + Lam7*RR23*(RR23*RR31 + RR21&
  &*RR33) + RR21*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR23*(RR23*RR32 + RR22*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f1233 = -1.D0*RR13*(2.D0*Lam7*RR21*RR31*RR33 + 2.D0*Lam8*RR22*RR32*RR33 + RR23*(Lam7*DBLE(RR31**INT(2.D0)) + Lam8*DBLE(R&
  &R32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR22*RR31*RR32 + 2.D0*Lam7*RR23*RR3&
  &1*RR33 + RR21*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR32**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1&
  &.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR31*RR32 + 2.D0*Lam8*RR23*RR32*RR33 + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS1S1S1S1f1311 = -1.D0*RR13*(2.D0*Lam7*RR11*RR13*RR31 + 2.D0*Lam8*RR12*RR13*RR32 + RR33*(Lam7*DBLE(RR11**INT(2.D0)) + Lam8*DBLE(R&
  &R12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR32 + 2.D0*Lam7*RR11*RR1&
  &3*RR33 + RR31*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR12**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1&
  &.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR31 + 2.D0*Lam8*RR12*RR13*RR33 + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS1S1S1S1f1312 = -1.D0*RR11*((3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)*RR31 + (Lam3 + Lam4 + Lam5)*&
  &(RR12*RR21 + RR11*RR22)*RR32 + Lam7*(RR13*RR21 + RR11*RR23)*RR33) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*(RR12*RR21 + RR11*RR22)*RR&
  &31 + ((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23)*RR32 + Lam8*(RR13*RR22 + RR12*RR23)*RR33) - 1.D0*&
  &RR13*(Lam7*(RR13*RR21 + RR11*RR23)*RR31 + Lam8*(RR13*RR22 + RR12*RR23)*RR32 + (Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23)*RR33)
CS1S1S1S1f1313 = -1.D0*RR13*(Lam7*RR31*(RR13*RR31 + RR11*RR33) + Lam8*RR32*(RR13*RR32 + RR12*RR33) + RR33*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR32*(RR12*RR31 + RR11*RR32) + Lam7*RR33*(RR13*RR31 + RR11&
  &*RR33) + RR31*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR33*(RR13*RR32 + RR12*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f1321 = -1.D0*RR11*((3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)*RR31 + (Lam3 + Lam4 + Lam5)*&
  &(RR12*RR21 + RR11*RR22)*RR32 + Lam7*(RR13*RR21 + RR11*RR23)*RR33) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*(RR12*RR21 + RR11*RR22)*RR&
  &31 + ((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23)*RR32 + Lam8*(RR13*RR22 + RR12*RR23)*RR33) - 1.D0*&
  &RR13*(Lam7*(RR13*RR21 + RR11*RR23)*RR31 + Lam8*(RR13*RR22 + RR12*RR23)*RR32 + (Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23)*RR33)
CS1S1S1S1f1322 = -1.D0*RR13*(2.D0*Lam7*RR21*RR23*RR31 + 2.D0*Lam8*RR22*RR23*RR32 + RR33*(Lam7*DBLE(RR21**INT(2.D0)) + Lam8*DBLE(R&
  &R22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR22*RR32 + 2.D0*Lam7*RR21*RR2&
  &3*RR33 + RR31*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR22**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1&
  &.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR22*RR31 + 2.D0*Lam8*RR22*RR23*RR33 + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS1S1S1S1f1323 = -1.D0*RR13*(Lam7*RR31*(RR23*RR31 + RR21*RR33) + Lam8*RR32*(RR23*RR32 + RR22*RR33) + RR33*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR32*(RR22*RR31 + RR21*RR32) + Lam7*RR33*(RR23*RR31 + RR21&
  &*RR33) + RR31*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR33*(RR23*RR32 + RR22*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f1331 = -1.D0*RR13*(Lam7*RR31*(RR13*RR31 + RR11*RR33) + Lam8*RR32*(RR13*RR32 + RR12*RR33) + RR33*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR32*(RR12*RR31 + RR11*RR32) + Lam7*RR33*(RR13*RR31 + RR11&
  &*RR33) + RR31*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR33*(RR13*RR32 + RR12*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f1332 = -1.D0*RR13*(Lam7*RR31*(RR23*RR31 + RR21*RR33) + Lam8*RR32*(RR23*RR32 + RR22*RR33) + RR33*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR32*(RR22*RR31 + RR21*RR32) + Lam7*RR33*(RR23*RR31 + RR21&
  &*RR33) + RR31*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR33*(RR23*RR32 + RR22*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f1333 = -1.D0*RR13*(2.D0*Lam7*RR33*DBLE(RR31**INT(2.D0)) + 2.D0*Lam8*RR33*DBLE(RR32**INT(2.D0)) + RR33*(Lam7*DBLE(RR31**&
  &INT(2.D0)) + Lam8*DBLE(RR32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR31*DBLE(R&
  &R32**INT(2.D0)) + 2.D0*Lam7*RR31*DBLE(RR33**INT(2.D0)) + RR31*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR3&
  &2**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR32*DBLE(RR31**INT(2.D0)) + 2.D0*Lam8*RR3&
  &2*DBLE(RR33**INT(2.D0)) + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33*&
  &*INT(2.D0))))
CS1S1S1S1f2111 = -1.D0*RR23*(2.D0*Lam7*RR13*DBLE(RR11**INT(2.D0)) + 2.D0*Lam8*RR13*DBLE(RR12**INT(2.D0)) + RR13*(Lam7*DBLE(RR11**&
  &INT(2.D0)) + Lam8*DBLE(RR12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*DBLE(R&
  &R12**INT(2.D0)) + 2.D0*Lam7*RR11*DBLE(RR13**INT(2.D0)) + RR11*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR1&
  &2**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*DBLE(RR11**INT(2.D0)) + 2.D0*Lam8*RR1&
  &2*DBLE(RR13**INT(2.D0)) + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13*&
  &*INT(2.D0))))
CS1S1S1S1f2112 = -1.D0*RR23*(Lam7*RR11*(RR13*RR21 + RR11*RR23) + Lam8*RR12*(RR13*RR22 + RR12*RR23) + RR13*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR21 + RR11*RR22) + Lam7*RR13*(RR13*RR21 + RR11&
  &*RR23) + RR11*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR13*(RR13*RR22 + RR12*RR23) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f2113 = -1.D0*RR23*(Lam7*RR11*(RR13*RR31 + RR11*RR33) + Lam8*RR12*(RR13*RR32 + RR12*RR33) + RR13*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR31 + RR11*RR32) + Lam7*RR13*(RR13*RR31 + RR11&
  &*RR33) + RR11*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR13*(RR13*RR32 + RR12*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f2121 = -1.D0*RR23*(Lam7*RR11*(RR13*RR21 + RR11*RR23) + Lam8*RR12*(RR13*RR22 + RR12*RR23) + RR13*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR21 + RR11*RR22) + Lam7*RR13*(RR13*RR21 + RR11&
  &*RR23) + RR11*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR13*(RR13*RR22 + RR12*RR23) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f2122 = -1.D0*RR23*(2.D0*Lam7*RR11*RR21*RR23 + 2.D0*Lam8*RR12*RR22*RR23 + RR13*(Lam7*DBLE(RR21**INT(2.D0)) + Lam8*DBLE(R&
  &R22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*RR21*RR22 + 2.D0*Lam7*RR13*RR2&
  &1*RR23 + RR11*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR22**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1&
  &.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR21*RR22 + 2.D0*Lam8*RR13*RR22*RR23 + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS1S1S1S1f2123 = -1.D0*RR23*(Lam7*RR11*(RR23*RR31 + RR21*RR33) + Lam8*RR12*(RR23*RR32 + RR22*RR33) + RR13*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*(RR22*RR31 + RR21*RR32) + Lam7*RR13*(RR23*RR31 + RR21&
  &*RR33) + RR11*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR13*(RR23*RR32 + RR22*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f2131 = -1.D0*RR23*(Lam7*RR11*(RR13*RR31 + RR11*RR33) + Lam8*RR12*(RR13*RR32 + RR12*RR33) + RR13*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR31 + RR11*RR32) + Lam7*RR13*(RR13*RR31 + RR11&
  &*RR33) + RR11*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR13*(RR13*RR32 + RR12*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f2132 = -1.D0*RR23*(Lam7*RR11*(RR23*RR31 + RR21*RR33) + Lam8*RR12*(RR23*RR32 + RR22*RR33) + RR13*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*(RR22*RR31 + RR21*RR32) + Lam7*RR13*(RR23*RR31 + RR21&
  &*RR33) + RR11*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR13*(RR23*RR32 + RR22*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f2133 = -1.D0*RR23*(2.D0*Lam7*RR11*RR31*RR33 + 2.D0*Lam8*RR12*RR32*RR33 + RR13*(Lam7*DBLE(RR31**INT(2.D0)) + Lam8*DBLE(R&
  &R32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*RR31*RR32 + 2.D0*Lam7*RR13*RR3&
  &1*RR33 + RR11*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR32**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1&
  &.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR31*RR32 + 2.D0*Lam8*RR13*RR32*RR33 + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS1S1S1S1f2211 = -1.D0*RR23*(2.D0*Lam7*RR11*RR13*RR21 + 2.D0*Lam8*RR12*RR13*RR22 + RR23*(Lam7*DBLE(RR11**INT(2.D0)) + Lam8*DBLE(R&
  &R12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR22 + 2.D0*Lam7*RR11*RR1&
  &3*RR23 + RR21*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR12**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1&
  &.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR21 + 2.D0*Lam8*RR12*RR13*RR23 + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS1S1S1S1f2212 = -1.D0*RR23*(Lam7*RR21*(RR13*RR21 + RR11*RR23) + Lam8*RR22*(RR13*RR22 + RR12*RR23) + RR23*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR21 + RR11*RR22) + Lam7*RR23*(RR13*RR21 + RR11&
  &*RR23) + RR21*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR23*(RR13*RR22 + RR12*RR23) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f2213 = -1.D0*RR23*(Lam7*RR21*(RR13*RR31 + RR11*RR33) + Lam8*RR22*(RR13*RR32 + RR12*RR33) + RR23*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR31 + RR11*RR32) + Lam7*RR23*(RR13*RR31 + RR11&
  &*RR33) + RR21*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR23*(RR13*RR32 + RR12*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f2221 = -1.D0*RR23*(Lam7*RR21*(RR13*RR21 + RR11*RR23) + Lam8*RR22*(RR13*RR22 + RR12*RR23) + RR23*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR21 + RR11*RR22) + Lam7*RR23*(RR13*RR21 + RR11&
  &*RR23) + RR21*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR23*(RR13*RR22 + RR12*RR23) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f2222 = -1.D0*RR23*(2.D0*Lam7*RR23*DBLE(RR21**INT(2.D0)) + 2.D0*Lam8*RR23*DBLE(RR22**INT(2.D0)) + RR23*(Lam7*DBLE(RR21**&
  &INT(2.D0)) + Lam8*DBLE(RR22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*DBLE(R&
  &R22**INT(2.D0)) + 2.D0*Lam7*RR21*DBLE(RR23**INT(2.D0)) + RR21*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR2&
  &2**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR22*DBLE(RR21**INT(2.D0)) + 2.D0*Lam8*RR2&
  &2*DBLE(RR23**INT(2.D0)) + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23*&
  &*INT(2.D0))))
CS1S1S1S1f2223 = -1.D0*RR23*(Lam7*RR21*(RR23*RR31 + RR21*RR33) + Lam8*RR22*(RR23*RR32 + RR22*RR33) + RR23*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*(RR22*RR31 + RR21*RR32) + Lam7*RR23*(RR23*RR31 + RR21&
  &*RR33) + RR21*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR23*(RR23*RR32 + RR22*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f2231 = -1.D0*RR23*(Lam7*RR21*(RR13*RR31 + RR11*RR33) + Lam8*RR22*(RR13*RR32 + RR12*RR33) + RR23*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR31 + RR11*RR32) + Lam7*RR23*(RR13*RR31 + RR11&
  &*RR33) + RR21*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR23*(RR13*RR32 + RR12*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f2232 = -1.D0*RR23*(Lam7*RR21*(RR23*RR31 + RR21*RR33) + Lam8*RR22*(RR23*RR32 + RR22*RR33) + RR23*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*(RR22*RR31 + RR21*RR32) + Lam7*RR23*(RR23*RR31 + RR21&
  &*RR33) + RR21*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR23*(RR23*RR32 + RR22*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f2233 = -1.D0*RR23*(2.D0*Lam7*RR21*RR31*RR33 + 2.D0*Lam8*RR22*RR32*RR33 + RR23*(Lam7*DBLE(RR31**INT(2.D0)) + Lam8*DBLE(R&
  &R32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR22*RR31*RR32 + 2.D0*Lam7*RR23*RR3&
  &1*RR33 + RR21*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR32**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1&
  &.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR31*RR32 + 2.D0*Lam8*RR23*RR32*RR33 + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS1S1S1S1f2311 = -1.D0*RR23*(2.D0*Lam7*RR11*RR13*RR31 + 2.D0*Lam8*RR12*RR13*RR32 + RR33*(Lam7*DBLE(RR11**INT(2.D0)) + Lam8*DBLE(R&
  &R12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR32 + 2.D0*Lam7*RR11*RR1&
  &3*RR33 + RR31*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR12**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1&
  &.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR31 + 2.D0*Lam8*RR12*RR13*RR33 + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS1S1S1S1f2312 = -1.D0*RR21*((3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)*RR31 + (Lam3 + Lam4 + Lam5)*&
  &(RR12*RR21 + RR11*RR22)*RR32 + Lam7*(RR13*RR21 + RR11*RR23)*RR33) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*(RR12*RR21 + RR11*RR22)*RR&
  &31 + ((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23)*RR32 + Lam8*(RR13*RR22 + RR12*RR23)*RR33) - 1.D0*&
  &RR23*(Lam7*(RR13*RR21 + RR11*RR23)*RR31 + Lam8*(RR13*RR22 + RR12*RR23)*RR32 + (Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23)*RR33)
CS1S1S1S1f2313 = -1.D0*RR23*(Lam7*RR31*(RR13*RR31 + RR11*RR33) + Lam8*RR32*(RR13*RR32 + RR12*RR33) + RR33*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR32*(RR12*RR31 + RR11*RR32) + Lam7*RR33*(RR13*RR31 + RR11&
  &*RR33) + RR31*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR33*(RR13*RR32 + RR12*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f2321 = -1.D0*RR21*((3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)*RR31 + (Lam3 + Lam4 + Lam5)*&
  &(RR12*RR21 + RR11*RR22)*RR32 + Lam7*(RR13*RR21 + RR11*RR23)*RR33) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*(RR12*RR21 + RR11*RR22)*RR&
  &31 + ((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23)*RR32 + Lam8*(RR13*RR22 + RR12*RR23)*RR33) - 1.D0*&
  &RR23*(Lam7*(RR13*RR21 + RR11*RR23)*RR31 + Lam8*(RR13*RR22 + RR12*RR23)*RR32 + (Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23)*RR33)
CS1S1S1S1f2322 = -1.D0*RR23*(2.D0*Lam7*RR21*RR23*RR31 + 2.D0*Lam8*RR22*RR23*RR32 + RR33*(Lam7*DBLE(RR21**INT(2.D0)) + Lam8*DBLE(R&
  &R22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR22*RR32 + 2.D0*Lam7*RR21*RR2&
  &3*RR33 + RR31*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR22**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1&
  &.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR22*RR31 + 2.D0*Lam8*RR22*RR23*RR33 + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS1S1S1S1f2323 = -1.D0*RR23*(Lam7*RR31*(RR23*RR31 + RR21*RR33) + Lam8*RR32*(RR23*RR32 + RR22*RR33) + RR33*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR32*(RR22*RR31 + RR21*RR32) + Lam7*RR33*(RR23*RR31 + RR21&
  &*RR33) + RR31*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR33*(RR23*RR32 + RR22*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f2331 = -1.D0*RR23*(Lam7*RR31*(RR13*RR31 + RR11*RR33) + Lam8*RR32*(RR13*RR32 + RR12*RR33) + RR33*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR32*(RR12*RR31 + RR11*RR32) + Lam7*RR33*(RR13*RR31 + RR11&
  &*RR33) + RR31*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR33*(RR13*RR32 + RR12*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f2332 = -1.D0*RR23*(Lam7*RR31*(RR23*RR31 + RR21*RR33) + Lam8*RR32*(RR23*RR32 + RR22*RR33) + RR33*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR32*(RR22*RR31 + RR21*RR32) + Lam7*RR33*(RR23*RR31 + RR21&
  &*RR33) + RR31*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR33*(RR23*RR32 + RR22*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f2333 = -1.D0*RR23*(2.D0*Lam7*RR33*DBLE(RR31**INT(2.D0)) + 2.D0*Lam8*RR33*DBLE(RR32**INT(2.D0)) + RR33*(Lam7*DBLE(RR31**&
  &INT(2.D0)) + Lam8*DBLE(RR32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR31*DBLE(R&
  &R32**INT(2.D0)) + 2.D0*Lam7*RR31*DBLE(RR33**INT(2.D0)) + RR31*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR3&
  &2**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR32*DBLE(RR31**INT(2.D0)) + 2.D0*Lam8*RR3&
  &2*DBLE(RR33**INT(2.D0)) + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33*&
  &*INT(2.D0))))
CS1S1S1S1f3111 = -1.D0*RR33*(2.D0*Lam7*RR13*DBLE(RR11**INT(2.D0)) + 2.D0*Lam8*RR13*DBLE(RR12**INT(2.D0)) + RR13*(Lam7*DBLE(RR11**&
  &INT(2.D0)) + Lam8*DBLE(RR12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*DBLE(R&
  &R12**INT(2.D0)) + 2.D0*Lam7*RR11*DBLE(RR13**INT(2.D0)) + RR11*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR1&
  &2**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*DBLE(RR11**INT(2.D0)) + 2.D0*Lam8*RR1&
  &2*DBLE(RR13**INT(2.D0)) + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13*&
  &*INT(2.D0))))
CS1S1S1S1f3112 = -1.D0*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR21 + RR11*RR22) + Lam7*RR13*(RR13*RR21 + RR11*RR23) + RR11*(3.D0*Lam1*R&
  &R11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23))*RR31 - 1.D0*((Lam3 + Lam4 + Lam5)*RR11*(RR12*RR21 + RR11*RR22) + &
  &Lam8*RR13*(RR13*RR22 + RR12*RR23) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23))*RR32 - 1.D0*&
  &(Lam7*RR11*(RR13*RR21 + RR11*RR23) + Lam8*RR12*(RR13*RR22 + RR12*RR23) + RR13*(Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23))*RR33
CS1S1S1S1f3113 = -1.D0*RR33*(Lam7*RR11*(RR13*RR31 + RR11*RR33) + Lam8*RR12*(RR13*RR32 + RR12*RR33) + RR13*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR31 + RR11*RR32) + Lam7*RR13*(RR13*RR31 + RR11&
  &*RR33) + RR11*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR13*(RR13*RR32 + RR12*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f3121 = -1.D0*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR21 + RR11*RR22) + Lam7*RR13*(RR13*RR21 + RR11*RR23) + RR11*(3.D0*Lam1*R&
  &R11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23))*RR31 - 1.D0*((Lam3 + Lam4 + Lam5)*RR11*(RR12*RR21 + RR11*RR22) + &
  &Lam8*RR13*(RR13*RR22 + RR12*RR23) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23))*RR32 - 1.D0*&
  &(Lam7*RR11*(RR13*RR21 + RR11*RR23) + Lam8*RR12*(RR13*RR22 + RR12*RR23) + RR13*(Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23))*RR33
CS1S1S1S1f3122 = -1.D0*RR33*(2.D0*Lam7*RR11*RR21*RR23 + 2.D0*Lam8*RR12*RR22*RR23 + RR13*(Lam7*DBLE(RR21**INT(2.D0)) + Lam8*DBLE(R&
  &R22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*RR21*RR22 + 2.D0*Lam7*RR13*RR2&
  &1*RR23 + RR11*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR22**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1&
  &.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR21*RR22 + 2.D0*Lam8*RR13*RR22*RR23 + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS1S1S1S1f3123 = -1.D0*RR33*(Lam7*RR11*(RR23*RR31 + RR21*RR33) + Lam8*RR12*(RR23*RR32 + RR22*RR33) + RR13*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR12*(RR22*RR31 + RR21*RR32) + Lam7*RR13*(RR23*RR31 + RR21&
  &*RR33) + RR11*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR13*(RR23*RR32 + RR22*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f3131 = -1.D0*RR33*(Lam7*RR11*(RR13*RR31 + RR11*RR33) + Lam8*RR12*(RR13*RR32 + RR12*RR33) + RR13*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR31 + RR11*RR32) + Lam7*RR13*(RR13*RR31 + RR11&
  &*RR33) + RR11*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR13*(RR13*RR32 + RR12*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f3132 = -1.D0*RR33*(Lam7*RR11*(RR23*RR31 + RR21*RR33) + Lam8*RR12*(RR23*RR32 + RR22*RR33) + RR13*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR12*(RR22*RR31 + RR21*RR32) + Lam7*RR13*(RR23*RR31 + RR21&
  &*RR33) + RR11*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR13*(RR23*RR32 + RR22*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f3133 = -1.D0*RR33*(2.D0*Lam7*RR11*RR31*RR33 + 2.D0*Lam8*RR12*RR32*RR33 + RR13*(Lam7*DBLE(RR31**INT(2.D0)) + Lam8*DBLE(R&
  &R32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*RR31*RR32 + 2.D0*Lam7*RR13*RR3&
  &1*RR33 + RR11*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR32**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1&
  &.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR31*RR32 + 2.D0*Lam8*RR13*RR32*RR33 + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS1S1S1S1f3211 = -1.D0*RR33*(2.D0*Lam7*RR11*RR13*RR21 + 2.D0*Lam8*RR12*RR13*RR22 + RR23*(Lam7*DBLE(RR11**INT(2.D0)) + Lam8*DBLE(R&
  &R12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR22 + 2.D0*Lam7*RR11*RR1&
  &3*RR23 + RR21*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR12**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1&
  &.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR21 + 2.D0*Lam8*RR12*RR13*RR23 + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS1S1S1S1f3212 = -1.D0*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR21 + RR11*RR22) + Lam7*RR23*(RR13*RR21 + RR11*RR23) + RR21*(3.D0*Lam1*R&
  &R11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23))*RR31 - 1.D0*((Lam3 + Lam4 + Lam5)*RR21*(RR12*RR21 + RR11*RR22) + &
  &Lam8*RR23*(RR13*RR22 + RR12*RR23) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23))*RR32 - 1.D0*&
  &(Lam7*RR21*(RR13*RR21 + RR11*RR23) + Lam8*RR22*(RR13*RR22 + RR12*RR23) + RR23*(Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23))*RR33
CS1S1S1S1f3213 = -1.D0*RR33*(Lam7*RR21*(RR13*RR31 + RR11*RR33) + Lam8*RR22*(RR13*RR32 + RR12*RR33) + RR23*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR31 + RR11*RR32) + Lam7*RR23*(RR13*RR31 + RR11&
  &*RR33) + RR21*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR23*(RR13*RR32 + RR12*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f3221 = -1.D0*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR21 + RR11*RR22) + Lam7*RR23*(RR13*RR21 + RR11*RR23) + RR21*(3.D0*Lam1*R&
  &R11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23))*RR31 - 1.D0*((Lam3 + Lam4 + Lam5)*RR21*(RR12*RR21 + RR11*RR22) + &
  &Lam8*RR23*(RR13*RR22 + RR12*RR23) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23))*RR32 - 1.D0*&
  &(Lam7*RR21*(RR13*RR21 + RR11*RR23) + Lam8*RR22*(RR13*RR22 + RR12*RR23) + RR23*(Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23))*RR33
CS1S1S1S1f3222 = -1.D0*RR33*(2.D0*Lam7*RR23*DBLE(RR21**INT(2.D0)) + 2.D0*Lam8*RR23*DBLE(RR22**INT(2.D0)) + RR23*(Lam7*DBLE(RR21**&
  &INT(2.D0)) + Lam8*DBLE(RR22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*DBLE(R&
  &R22**INT(2.D0)) + 2.D0*Lam7*RR21*DBLE(RR23**INT(2.D0)) + RR21*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR2&
  &2**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR22*DBLE(RR21**INT(2.D0)) + 2.D0*Lam8*RR2&
  &2*DBLE(RR23**INT(2.D0)) + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23*&
  &*INT(2.D0))))
CS1S1S1S1f3223 = -1.D0*RR33*(Lam7*RR21*(RR23*RR31 + RR21*RR33) + Lam8*RR22*(RR23*RR32 + RR22*RR33) + RR23*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR22*(RR22*RR31 + RR21*RR32) + Lam7*RR23*(RR23*RR31 + RR21&
  &*RR33) + RR21*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR23*(RR23*RR32 + RR22*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f3231 = -1.D0*RR33*(Lam7*RR21*(RR13*RR31 + RR11*RR33) + Lam8*RR22*(RR13*RR32 + RR12*RR33) + RR23*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR31 + RR11*RR32) + Lam7*RR23*(RR13*RR31 + RR11&
  &*RR33) + RR21*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR23*(RR13*RR32 + RR12*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f3232 = -1.D0*RR33*(Lam7*RR21*(RR23*RR31 + RR21*RR33) + Lam8*RR22*(RR23*RR32 + RR22*RR33) + RR23*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR22*(RR22*RR31 + RR21*RR32) + Lam7*RR23*(RR23*RR31 + RR21&
  &*RR33) + RR21*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR23*(RR23*RR32 + RR22*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f3233 = -1.D0*RR33*(2.D0*Lam7*RR21*RR31*RR33 + 2.D0*Lam8*RR22*RR32*RR33 + RR23*(Lam7*DBLE(RR31**INT(2.D0)) + Lam8*DBLE(R&
  &R32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR22*RR31*RR32 + 2.D0*Lam7*RR23*RR3&
  &1*RR33 + RR21*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR32**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1&
  &.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR31*RR32 + 2.D0*Lam8*RR23*RR32*RR33 + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS1S1S1S1f3311 = -1.D0*RR33*(2.D0*Lam7*RR11*RR13*RR31 + 2.D0*Lam8*RR12*RR13*RR32 + RR33*(Lam7*DBLE(RR11**INT(2.D0)) + Lam8*DBLE(R&
  &R12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR32 + 2.D0*Lam7*RR11*RR1&
  &3*RR33 + RR31*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR12**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1&
  &.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR31 + 2.D0*Lam8*RR12*RR13*RR33 + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS1S1S1S1f3312 = -1.D0*RR31*((3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)*RR31 + (Lam3 + Lam4 + Lam5)*&
  &(RR12*RR21 + RR11*RR22)*RR32 + Lam7*(RR13*RR21 + RR11*RR23)*RR33) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*(RR12*RR21 + RR11*RR22)*RR&
  &31 + ((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23)*RR32 + Lam8*(RR13*RR22 + RR12*RR23)*RR33) - 1.D0*&
  &RR33*(Lam7*(RR13*RR21 + RR11*RR23)*RR31 + Lam8*(RR13*RR22 + RR12*RR23)*RR32 + (Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23)*RR33)
CS1S1S1S1f3313 = -1.D0*RR33*(Lam7*RR31*(RR13*RR31 + RR11*RR33) + Lam8*RR32*(RR13*RR32 + RR12*RR33) + RR33*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR32*(RR12*RR31 + RR11*RR32) + Lam7*RR33*(RR13*RR31 + RR11&
  &*RR33) + RR31*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR33*(RR13*RR32 + RR12*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f3321 = -1.D0*RR31*((3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)*RR31 + (Lam3 + Lam4 + Lam5)*&
  &(RR12*RR21 + RR11*RR22)*RR32 + Lam7*(RR13*RR21 + RR11*RR23)*RR33) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*(RR12*RR21 + RR11*RR22)*RR&
  &31 + ((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23)*RR32 + Lam8*(RR13*RR22 + RR12*RR23)*RR33) - 1.D0*&
  &RR33*(Lam7*(RR13*RR21 + RR11*RR23)*RR31 + Lam8*(RR13*RR22 + RR12*RR23)*RR32 + (Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23)*RR33)
CS1S1S1S1f3322 = -1.D0*RR33*(2.D0*Lam7*RR21*RR23*RR31 + 2.D0*Lam8*RR22*RR23*RR32 + RR33*(Lam7*DBLE(RR21**INT(2.D0)) + Lam8*DBLE(R&
  &R22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR22*RR32 + 2.D0*Lam7*RR21*RR2&
  &3*RR33 + RR31*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR22**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1&
  &.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR22*RR31 + 2.D0*Lam8*RR22*RR23*RR33 + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS1S1S1S1f3323 = -1.D0*RR33*(Lam7*RR31*(RR23*RR31 + RR21*RR33) + Lam8*RR32*(RR23*RR32 + RR22*RR33) + RR33*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR32*(RR22*RR31 + RR21*RR32) + Lam7*RR33*(RR23*RR31 + RR21&
  &*RR33) + RR31*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR33*(RR23*RR32 + RR22*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f3331 = -1.D0*RR33*(Lam7*RR31*(RR13*RR31 + RR11*RR33) + Lam8*RR32*(RR13*RR32 + RR12*RR33) + RR33*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR32*(RR12*RR31 + RR11*RR32) + Lam7*RR33*(RR13*RR31 + RR11&
  &*RR33) + RR31*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR33*(RR13*RR32 + RR12*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f3332 = -1.D0*RR33*(Lam7*RR31*(RR23*RR31 + RR21*RR33) + Lam8*RR32*(RR23*RR32 + RR22*RR33) + RR33*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR32*(RR22*RR31 + RR21*RR32) + Lam7*RR33*(RR23*RR31 + RR21&
  &*RR33) + RR31*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR33*(RR23*RR32 + RR22*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f3333 = -1.D0*RR33*(2.D0*Lam7*RR33*DBLE(RR31**INT(2.D0)) + 2.D0*Lam8*RR33*DBLE(RR32**INT(2.D0)) + RR33*(Lam7*DBLE(RR31**&
  &INT(2.D0)) + Lam8*DBLE(RR32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR31*DBLE(R&
  &R32**INT(2.D0)) + 2.D0*Lam7*RR31*DBLE(RR33**INT(2.D0)) + RR31*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR3&
  &2**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR32*DBLE(RR31**INT(2.D0)) + 2.D0*Lam8*RR3&
  &2*DBLE(RR33**INT(2.D0)) + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33*&
  &*INT(2.D0))))

CS2S2S2S2f1111 = -1.D0*SB*(2.D0*CB2*(Lam3 + Lam4 + Lam5)*SB + SB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam2*SB2)) - 1.D0*CB*(2.D0*CB*(&
  &Lam3 + Lam4 + Lam5)*SB2 + CB*(3.D0*CB2*Lam1 + (Lam3 + Lam4 + Lam5)*SB2))
CS2S2S2S2f1112 = -1.D0*SB*(SB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1&
  &.D0*CB*(CB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + (Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f1121 = -1.D0*SB*(SB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1&
  &.D0*CB*(CB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + (Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f1122 = -1.D0*CB*(-2.D0*CB*(Lam3 + Lam4 + Lam5)*SB2 + CB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam1*SB2)) - 1.D0*SB*(-2.D0*CB&
  &2*(Lam3 + Lam4 + Lam5)*SB + SB*(3.D0*CB2*Lam2 + (Lam3 + Lam4 + Lam5)*SB2))
CS2S2S2S2f1211 = -1.D0*SB*(-2.D0*CB*(Lam3 + Lam4 + Lam5)*SB2 + CB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam2*SB2)) - 1.D0*CB*(2.D0*CB2&
  &*(Lam3 + Lam4 + Lam5)*SB - 1.D0*SB*(3.D0*CB2*Lam1 + (Lam3 + Lam4 + Lam5)*SB2))
CS2S2S2S2f1212 = -1.D0*CB*(-1.D0*SB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) -&
  & 1.D0*SB*(CB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) - 1.D0*(Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f1221 = -1.D0*CB*(-1.D0*SB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) -&
  & 1.D0*SB*(CB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) - 1.D0*(Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f1222 = -1.D0*CB*(-2.D0*CB2*(Lam3 + Lam4 + Lam5)*SB - 1.D0*SB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam1*SB2)) - 1.D0*SB*(2.D&
  &0*CB*(Lam3 + Lam4 + Lam5)*SB2 + CB*(3.D0*CB2*Lam2 + (Lam3 + Lam4 + Lam5)*SB2))
CS2S2S2S2f2111 = -1.D0*CB*(2.D0*CB2*(Lam3 + Lam4 + Lam5)*SB + SB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam2*SB2)) + SB*(2.D0*CB*(Lam3 &
  &+ Lam4 + Lam5)*SB2 + CB*(3.D0*CB2*Lam1 + (Lam3 + Lam4 + Lam5)*SB2))
CS2S2S2S2f2112 = -1.D0*CB*(SB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) + S&
  &B*(CB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + (Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f2121 = -1.D0*CB*(SB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) + S&
  &B*(CB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + (Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f2122 = SB*(-2.D0*CB*(Lam3 + Lam4 + Lam5)*SB2 + CB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam1*SB2)) - 1.D0*CB*(-2.D0*CB2*(Lam&
  &3 + Lam4 + Lam5)*SB + SB*(3.D0*CB2*Lam2 + (Lam3 + Lam4 + Lam5)*SB2))
CS2S2S2S2f2211 = -1.D0*CB*(-2.D0*CB*(Lam3 + Lam4 + Lam5)*SB2 + CB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam2*SB2)) + SB*(2.D0*CB2*(Lam&
  &3 + Lam4 + Lam5)*SB - 1.D0*SB*(3.D0*CB2*Lam1 + (Lam3 + Lam4 + Lam5)*SB2))
CS2S2S2S2f2212 = SB*(-1.D0*SB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*&
  &CB*(CB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) - 1.D0*(Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f2221 = SB*(-1.D0*SB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*&
  &CB*(CB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) - 1.D0*(Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f2222 = SB*(-2.D0*CB2*(Lam3 + Lam4 + Lam5)*SB - 1.D0*SB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam1*SB2)) - 1.D0*CB*(2.D0*CB*(&
  &Lam3 + Lam4 + Lam5)*SB2 + CB*(3.D0*CB2*Lam2 + (Lam3 + Lam4 + Lam5)*SB2))

CS3S3S3S3f1111 = -1.D0*SB*(2.D0*CB2*(Lam3 + Lam4)*SB + 2.D0*SB*(CB2*Lam5 + Lam2*SB2)) - 1.D0*CB*(2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D0&
  &*CB*(CB2*Lam1 + Lam5*SB2))
CS3S3S3S3f1112 = -1.D0*SB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*(-1.D0*&
  &CB*Lam1*SB + CB*Lam5*SB) + (Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f1121 = -1.D0*SB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*(-1.D0*&
  &CB*Lam1*SB + CB*Lam5*SB) + (Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f1122 = -1.D0*CB*(-2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D0*CB*(CB2*Lam5 + Lam1*SB2)) - 1.D0*SB*(-2.D0*CB2*(Lam3 + Lam4)*SB + 2.&
  &D0*SB*(CB2*Lam2 + Lam5*SB2))
CS3S3S3S3f1211 = -1.D0*SB*(-2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D0*CB*(CB2*Lam5 + Lam2*SB2)) - 1.D0*CB*(2.D0*CB2*(Lam3 + Lam4)*SB - 2.D&
  &0*SB*(CB2*Lam1 + Lam5*SB2))
CS3S3S3S3f1212 = -1.D0*CB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) - 1.D0*SB*(2.D0*CB*(CB*L&
  &am2*SB - 1.D0*CB*Lam5*SB) - 1.D0*(Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f1221 = -1.D0*CB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) - 1.D0*SB*(2.D0*CB*(CB*L&
  &am2*SB - 1.D0*CB*Lam5*SB) - 1.D0*(Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f1222 = -1.D0*CB*(-2.D0*CB2*(Lam3 + Lam4)*SB - 2.D0*SB*(CB2*Lam5 + Lam1*SB2)) - 1.D0*SB*(2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D&
  &0*CB*(CB2*Lam2 + Lam5*SB2))
CS3S3S3S3f2111 = -1.D0*CB*(2.D0*CB2*(Lam3 + Lam4)*SB + 2.D0*SB*(CB2*Lam5 + Lam2*SB2)) + SB*(2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D0*CB*(&
  &CB2*Lam1 + Lam5*SB2))
CS3S3S3S3f2112 = -1.D0*CB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) + SB*(2.D0*CB*(-1.D0*CB*La&
  &m1*SB + CB*Lam5*SB) + (Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f2121 = -1.D0*CB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) + SB*(2.D0*CB*(-1.D0*CB*La&
  &m1*SB + CB*Lam5*SB) + (Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f2122 = SB*(-2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D0*CB*(CB2*Lam5 + Lam1*SB2)) - 1.D0*CB*(-2.D0*CB2*(Lam3 + Lam4)*SB + 2.D0*SB*&
  &(CB2*Lam2 + Lam5*SB2))
CS3S3S3S3f2211 = -1.D0*CB*(-2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D0*CB*(CB2*Lam5 + Lam2*SB2)) + SB*(2.D0*CB2*(Lam3 + Lam4)*SB - 2.D0*SB*&
  &(CB2*Lam1 + Lam5*SB2))
CS3S3S3S3f2212 = SB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*(CB*Lam2*SB&
  & - 1.D0*CB*Lam5*SB) - 1.D0*(Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f2221 = SB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*(CB*Lam2*SB&
  & - 1.D0*CB*Lam5*SB) - 1.D0*(Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f2222 = SB*(-2.D0*CB2*(Lam3 + Lam4)*SB - 2.D0*SB*(CB2*Lam5 + Lam1*SB2)) - 1.D0*CB*(2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D0*CB*(&
  &CB2*Lam2 + Lam5*SB2))

CS2S2S1S1f1111 = -1.D0*CB*(2.D0*Lam5*RR11*RR12*SB + CB*(Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR12**INT(2.D&
  &0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1.D0*SB*(2.D0*CB*Lam5*RR11*RR12 + SB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR11**INT(2.D0)) + L&
  &am2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS2S2S1S1f1112 = -1.D0*CB*(CB*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR13*RR23) + Lam5*(RR12*RR21 + RR11*RR&
  &22)*SB) - 1.D0*SB*(CB*Lam5*(RR12*RR21 + RR11*RR22) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23)*S&
  &B)
CS2S2S1S1f1113 = -1.D0*CB*(CB*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR13*RR33) + Lam5*(RR12*RR31 + RR11*RR&
  &32)*SB) - 1.D0*SB*(CB*Lam5*(RR12*RR31 + RR11*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33)*S&
  &B)
CS2S2S1S1f1121 = -1.D0*CB*(CB*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR13*RR23) + Lam5*(RR12*RR21 + RR11*RR&
  &22)*SB) - 1.D0*SB*(CB*Lam5*(RR12*RR21 + RR11*RR22) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23)*S&
  &B)
CS2S2S1S1f1122 = -1.D0*CB*(2.D0*Lam5*RR21*RR22*SB + CB*(Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR22**INT(2.D&
  &0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1.D0*SB*(2.D0*CB*Lam5*RR21*RR22 + SB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR21**INT(2.D0)) + L&
  &am2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS2S2S1S1f1123 = -1.D0*CB*(CB*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR23*RR33) + Lam5*(RR22*RR31 + RR21*RR&
  &32)*SB) - 1.D0*SB*(CB*Lam5*(RR22*RR31 + RR21*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33)*S&
  &B)
CS2S2S1S1f1131 = -1.D0*CB*(CB*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR13*RR33) + Lam5*(RR12*RR31 + RR11*RR&
  &32)*SB) - 1.D0*SB*(CB*Lam5*(RR12*RR31 + RR11*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33)*S&
  &B)
CS2S2S1S1f1132 = -1.D0*CB*(CB*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR23*RR33) + Lam5*(RR22*RR31 + RR21*RR&
  &32)*SB) - 1.D0*SB*(CB*Lam5*(RR22*RR31 + RR21*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33)*S&
  &B)
CS2S2S1S1f1133 = -1.D0*CB*(2.D0*Lam5*RR31*RR32*SB + CB*(Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR32**INT(2.D&
  &0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1.D0*SB*(2.D0*CB*Lam5*RR31*RR32 + SB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR31**INT(2.D0)) + L&
  &am2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS2S2S1S1f1211 = -1.D0*CB*(2.D0*CB*Lam5*RR11*RR12 - 1.D0*SB*(Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR12**IN&
  &T(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1.D0*SB*(-2.D0*Lam5*RR11*RR12*SB + CB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR11**INT(2.D0&
  &)) + Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS2S2S1S1f1212 = -1.D0*SB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23) - 1.D0*Lam5*(RR12*RR21 + RR&
  &11*RR22)*SB) - 1.D0*CB*(CB*Lam5*(RR12*RR21 + RR11*RR22) - 1.D0*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR&
  &13*RR23)*SB)
CS2S2S1S1f1213 = -1.D0*SB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33) - 1.D0*Lam5*(RR12*RR31 + RR&
  &11*RR32)*SB) - 1.D0*CB*(CB*Lam5*(RR12*RR31 + RR11*RR32) - 1.D0*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR&
  &13*RR33)*SB)
CS2S2S1S1f1221 = -1.D0*SB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23) - 1.D0*Lam5*(RR12*RR21 + RR&
  &11*RR22)*SB) - 1.D0*CB*(CB*Lam5*(RR12*RR21 + RR11*RR22) - 1.D0*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR&
  &13*RR23)*SB)
CS2S2S1S1f1222 = -1.D0*CB*(2.D0*CB*Lam5*RR21*RR22 - 1.D0*SB*(Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR22**IN&
  &T(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1.D0*SB*(-2.D0*Lam5*RR21*RR22*SB + CB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR21**INT(2.D0&
  &)) + Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS2S2S1S1f1223 = -1.D0*SB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33) - 1.D0*Lam5*(RR22*RR31 + RR&
  &21*RR32)*SB) - 1.D0*CB*(CB*Lam5*(RR22*RR31 + RR21*RR32) - 1.D0*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR&
  &23*RR33)*SB)
CS2S2S1S1f1231 = -1.D0*SB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33) - 1.D0*Lam5*(RR12*RR31 + RR&
  &11*RR32)*SB) - 1.D0*CB*(CB*Lam5*(RR12*RR31 + RR11*RR32) - 1.D0*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR&
  &13*RR33)*SB)
CS2S2S1S1f1232 = -1.D0*SB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33) - 1.D0*Lam5*(RR22*RR31 + RR&
  &21*RR32)*SB) - 1.D0*CB*(CB*Lam5*(RR22*RR31 + RR21*RR32) - 1.D0*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR&
  &23*RR33)*SB)
CS2S2S1S1f1233 = -1.D0*CB*(2.D0*CB*Lam5*RR31*RR32 - 1.D0*SB*(Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR32**IN&
  &T(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1.D0*SB*(-2.D0*Lam5*RR31*RR32*SB + CB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR31**INT(2.D0&
  &)) + Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS2S2S1S1f2111 = SB*(2.D0*Lam5*RR11*RR12*SB + CB*(Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR12**INT(2.D0)) + &
  &Lam7*DBLE(RR13**INT(2.D0)))) - 1.D0*CB*(2.D0*CB*Lam5*RR11*RR12 + SB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR11**INT(2.D0)) + Lam2*DB&
  &LE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS2S2S1S1f2112 = SB*(CB*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR13*RR23) + Lam5*(RR12*RR21 + RR11*RR22)*SB&
  &) - 1.D0*CB*(CB*Lam5*(RR12*RR21 + RR11*RR22) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23)*SB)
CS2S2S1S1f2113 = SB*(CB*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR13*RR33) + Lam5*(RR12*RR31 + RR11*RR32)*SB&
  &) - 1.D0*CB*(CB*Lam5*(RR12*RR31 + RR11*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33)*SB)
CS2S2S1S1f2121 = SB*(CB*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR13*RR23) + Lam5*(RR12*RR21 + RR11*RR22)*SB&
  &) - 1.D0*CB*(CB*Lam5*(RR12*RR21 + RR11*RR22) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23)*SB)
CS2S2S1S1f2122 = SB*(2.D0*Lam5*RR21*RR22*SB + CB*(Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR22**INT(2.D0)) + &
  &Lam7*DBLE(RR23**INT(2.D0)))) - 1.D0*CB*(2.D0*CB*Lam5*RR21*RR22 + SB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR21**INT(2.D0)) + Lam2*DB&
  &LE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS2S2S1S1f2123 = SB*(CB*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR23*RR33) + Lam5*(RR22*RR31 + RR21*RR32)*SB&
  &) - 1.D0*CB*(CB*Lam5*(RR22*RR31 + RR21*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33)*SB)
CS2S2S1S1f2131 = SB*(CB*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR13*RR33) + Lam5*(RR12*RR31 + RR11*RR32)*SB&
  &) - 1.D0*CB*(CB*Lam5*(RR12*RR31 + RR11*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33)*SB)
CS2S2S1S1f2132 = SB*(CB*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR23*RR33) + Lam5*(RR22*RR31 + RR21*RR32)*SB&
  &) - 1.D0*CB*(CB*Lam5*(RR22*RR31 + RR21*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33)*SB)
CS2S2S1S1f2133 = SB*(2.D0*Lam5*RR31*RR32*SB + CB*(Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR32**INT(2.D0)) + &
  &Lam7*DBLE(RR33**INT(2.D0)))) - 1.D0*CB*(2.D0*CB*Lam5*RR31*RR32 + SB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR31**INT(2.D0)) + Lam2*DB&
  &LE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS2S2S1S1f2211 = SB*(2.D0*CB*Lam5*RR11*RR12 - 1.D0*SB*(Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR12**INT(2.D0&
  &)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1.D0*CB*(-2.D0*Lam5*RR11*RR12*SB + CB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR11**INT(2.D0)) + L&
  &am2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS2S2S1S1f2212 = -1.D0*CB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23) - 1.D0*Lam5*(RR12*RR21 + RR&
  &11*RR22)*SB) + SB*(CB*Lam5*(RR12*RR21 + RR11*RR22) - 1.D0*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR13*RR&
  &23)*SB)
CS2S2S1S1f2213 = -1.D0*CB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33) - 1.D0*Lam5*(RR12*RR31 + RR&
  &11*RR32)*SB) + SB*(CB*Lam5*(RR12*RR31 + RR11*RR32) - 1.D0*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR13*RR&
  &33)*SB)
CS2S2S1S1f2221 = -1.D0*CB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23) - 1.D0*Lam5*(RR12*RR21 + RR&
  &11*RR22)*SB) + SB*(CB*Lam5*(RR12*RR21 + RR11*RR22) - 1.D0*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR13*RR&
  &23)*SB)
CS2S2S1S1f2222 = SB*(2.D0*CB*Lam5*RR21*RR22 - 1.D0*SB*(Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR22**INT(2.D0&
  &)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1.D0*CB*(-2.D0*Lam5*RR21*RR22*SB + CB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR21**INT(2.D0)) + L&
  &am2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS2S2S1S1f2223 = -1.D0*CB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33) - 1.D0*Lam5*(RR22*RR31 + RR&
  &21*RR32)*SB) + SB*(CB*Lam5*(RR22*RR31 + RR21*RR32) - 1.D0*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR23*RR&
  &33)*SB)
CS2S2S1S1f2231 = -1.D0*CB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33) - 1.D0*Lam5*(RR12*RR31 + RR&
  &11*RR32)*SB) + SB*(CB*Lam5*(RR12*RR31 + RR11*RR32) - 1.D0*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR13*RR&
  &33)*SB)
CS2S2S1S1f2232 = -1.D0*CB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33) - 1.D0*Lam5*(RR22*RR31 + RR&
  &21*RR32)*SB) + SB*(CB*Lam5*(RR22*RR31 + RR21*RR32) - 1.D0*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR23*RR&
  &33)*SB)
CS2S2S1S1f2233 = SB*(2.D0*CB*Lam5*RR31*RR32 - 1.D0*SB*(Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR32**INT(2.D0&
  &)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1.D0*CB*(-2.D0*Lam5*RR31*RR32*SB + CB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR31**INT(2.D0)) + L&
  &am2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))

CS2S2S3S3f1111 = 0.5D0*(-1.D0*SB*(2.D0*CB2*(Lam4 + Lam5)*SB + 2.D0*SB*(CB2*Lam3 + Lam2*SB2)) - 1.D0*CB*(2.D0*CB*(Lam4 + Lam5)*SB2&
  & + 2.D0*CB*(CB2*Lam1 + Lam3*SB2)))
CS2S2S3S3f1112 = 0.5D0*(-1.D0*SB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*&
  &(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f1121 = 0.5D0*(-1.D0*SB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*&
  &(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f1122 = 0.5D0*(-1.D0*CB*(-2.D0*CB*(Lam4 + Lam5)*SB2 + 2.D0*CB*(CB2*Lam3 + Lam1*SB2)) - 1.D0*SB*(-2.D0*CB2*(Lam4 + Lam5)*&
  &SB + 2.D0*SB*(CB2*Lam2 + Lam3*SB2)))
CS2S2S3S3f1211 = 0.5D0*(-1.D0*SB*(-2.D0*CB*(Lam4 + Lam5)*SB2 + 2.D0*CB*(CB2*Lam3 + Lam2*SB2)) - 1.D0*CB*(2.D0*CB2*(Lam4 + Lam5)*S&
  &B - 2.D0*SB*(CB2*Lam1 + Lam3*SB2)))
CS2S2S3S3f1212 = 0.5D0*(-1.D0*CB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*SB*(2.D0*C&
  &B*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) - 1.D0*(Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f1221 = 0.5D0*(-1.D0*CB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*SB*(2.D0*C&
  &B*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) - 1.D0*(Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f1222 = 0.5D0*(-1.D0*CB*(-2.D0*CB2*(Lam4 + Lam5)*SB - 2.D0*SB*(CB2*Lam3 + Lam1*SB2)) - 1.D0*SB*(2.D0*CB*(Lam4 + Lam5)*SB&
  &2 + 2.D0*CB*(CB2*Lam2 + Lam3*SB2)))
CS2S2S3S3f2111 = 0.5D0*(-1.D0*CB*(2.D0*CB2*(Lam4 + Lam5)*SB + 2.D0*SB*(CB2*Lam3 + Lam2*SB2)) + SB*(2.D0*CB*(Lam4 + Lam5)*SB2 + 2.&
  &D0*CB*(CB2*Lam1 + Lam3*SB2)))
CS2S2S3S3f2112 = 0.5D0*(-1.D0*CB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) + SB*(2.D0*CB*(-1.D&
  &0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f2121 = 0.5D0*(-1.D0*CB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) + SB*(2.D0*CB*(-1.D&
  &0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f2122 = 0.5D0*(SB*(-2.D0*CB*(Lam4 + Lam5)*SB2 + 2.D0*CB*(CB2*Lam3 + Lam1*SB2)) - 1.D0*CB*(-2.D0*CB2*(Lam4 + Lam5)*SB + 2&
  &.D0*SB*(CB2*Lam2 + Lam3*SB2)))
CS2S2S3S3f2211 = 0.5D0*(-1.D0*CB*(-2.D0*CB*(Lam4 + Lam5)*SB2 + 2.D0*CB*(CB2*Lam3 + Lam2*SB2)) + SB*(2.D0*CB2*(Lam4 + Lam5)*SB - 2&
  &.D0*SB*(CB2*Lam1 + Lam3*SB2)))
CS2S2S3S3f2212 = 0.5D0*(SB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*(CB*&
  &Lam2*SB - 1.D0*CB*Lam3*SB) - 1.D0*(Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f2221 = 0.5D0*(SB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*(CB*&
  &Lam2*SB - 1.D0*CB*Lam3*SB) - 1.D0*(Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f2222 = 0.5D0*(SB*(-2.D0*CB2*(Lam4 + Lam5)*SB - 2.D0*SB*(CB2*Lam3 + Lam1*SB2)) - 1.D0*CB*(2.D0*CB*(Lam4 + Lam5)*SB2 + 2.&
  &D0*CB*(CB2*Lam2 + Lam3*SB2)))

CS1S1S3S3f1111 = 0.5D0*(-1.D0*RR12*(2.D0*CB*(Lam4 + Lam5)*RR11*SB + 2.D0*RR12*(CB2*Lam3 + Lam2*SB2)) - 1.D0*RR11*(2.D0*CB*(Lam4 +&
  & Lam5)*RR12*SB + 2.D0*RR11*(CB2*Lam1 + Lam3*SB2)) - 2.D0*(CB2*Lam7 + Lam8*SB2)*DBLE(RR13**INT(2.D0)))
CS1S1S3S3f1112 = 0.5D0*(-1.D0*RR12*(2.D0*RR12*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (Lam4 + Lam5)*RR11*(CB2 - 1.D0*SB2)) - 1.D0*RR11*(&
  &2.D0*RR11*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR12*(CB2 - 1.D0*SB2)) - 2.D0*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*DBLE(R&
  &R13**INT(2.D0)))
CS1S1S3S3f1121 = 0.5D0*(-1.D0*RR12*(2.D0*RR12*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (Lam4 + Lam5)*RR11*(CB2 - 1.D0*SB2)) - 1.D0*RR11*(&
  &2.D0*RR11*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR12*(CB2 - 1.D0*SB2)) - 2.D0*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*DBLE(R&
  &R13**INT(2.D0)))
CS1S1S3S3f1122 = 0.5D0*(-1.D0*RR11*(-2.D0*CB*(Lam4 + Lam5)*RR12*SB + 2.D0*RR11*(CB2*Lam3 + Lam1*SB2)) - 1.D0*RR12*(-2.D0*CB*(Lam4&
  & + Lam5)*RR11*SB + 2.D0*RR12*(CB2*Lam2 + Lam3*SB2)) - 2.D0*(CB2*Lam8 + Lam7*SB2)*DBLE(RR13**INT(2.D0)))
CS1S1S3S3f1211 = 0.5D0*(-2.D0*RR13*RR23*(CB2*Lam7 + Lam8*SB2) - 1.D0*RR12*(2.D0*CB*(Lam4 + Lam5)*RR21*SB + 2.D0*RR22*(CB2*Lam3 + &
  &Lam2*SB2)) - 1.D0*RR11*(2.D0*CB*(Lam4 + Lam5)*RR22*SB + 2.D0*RR21*(CB2*Lam1 + Lam3*SB2)))
CS1S1S3S3f1212 = 0.5D0*(-2.D0*RR13*RR23*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR12*(2.D0*RR22*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR21*(CB2 - 1.D0*SB2)) - 1.D0*RR11*(2.D0*RR21*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR22*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f1221 = 0.5D0*(-2.D0*RR13*RR23*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR12*(2.D0*RR22*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR21*(CB2 - 1.D0*SB2)) - 1.D0*RR11*(2.D0*RR21*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR22*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f1222 = 0.5D0*(-2.D0*RR13*RR23*(CB2*Lam8 + Lam7*SB2) - 1.D0*RR11*(-2.D0*CB*(Lam4 + Lam5)*RR22*SB + 2.D0*RR21*(CB2*Lam3 +&
  & Lam1*SB2)) - 1.D0*RR12*(-2.D0*CB*(Lam4 + Lam5)*RR21*SB + 2.D0*RR22*(CB2*Lam2 + Lam3*SB2)))
CS1S1S3S3f1311 = 0.5D0*(-2.D0*RR13*RR33*(CB2*Lam7 + Lam8*SB2) - 1.D0*RR12*(2.D0*CB*(Lam4 + Lam5)*RR31*SB + 2.D0*RR32*(CB2*Lam3 + &
  &Lam2*SB2)) - 1.D0*RR11*(2.D0*CB*(Lam4 + Lam5)*RR32*SB + 2.D0*RR31*(CB2*Lam1 + Lam3*SB2)))
CS1S1S3S3f1312 = 0.5D0*(-2.D0*RR13*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR12*(2.D0*RR32*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR31*(CB2 - 1.D0*SB2)) - 1.D0*RR11*(2.D0*RR31*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR32*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f1321 = 0.5D0*(-2.D0*RR13*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR12*(2.D0*RR32*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR31*(CB2 - 1.D0*SB2)) - 1.D0*RR11*(2.D0*RR31*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR32*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f1322 = 0.5D0*(-2.D0*RR13*RR33*(CB2*Lam8 + Lam7*SB2) - 1.D0*RR11*(-2.D0*CB*(Lam4 + Lam5)*RR32*SB + 2.D0*RR31*(CB2*Lam3 +&
  & Lam1*SB2)) - 1.D0*RR12*(-2.D0*CB*(Lam4 + Lam5)*RR31*SB + 2.D0*RR32*(CB2*Lam2 + Lam3*SB2)))
CS1S1S3S3f2111 = 0.5D0*(-2.D0*RR13*RR23*(CB2*Lam7 + Lam8*SB2) - 1.D0*RR22*(2.D0*CB*(Lam4 + Lam5)*RR11*SB + 2.D0*RR12*(CB2*Lam3 + &
  &Lam2*SB2)) - 1.D0*RR21*(2.D0*CB*(Lam4 + Lam5)*RR12*SB + 2.D0*RR11*(CB2*Lam1 + Lam3*SB2)))
CS1S1S3S3f2112 = 0.5D0*(-2.D0*RR13*RR23*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR22*(2.D0*RR12*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR11*(CB2 - 1.D0*SB2)) - 1.D0*RR21*(2.D0*RR11*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR12*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f2121 = 0.5D0*(-2.D0*RR13*RR23*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR22*(2.D0*RR12*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR11*(CB2 - 1.D0*SB2)) - 1.D0*RR21*(2.D0*RR11*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR12*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f2122 = 0.5D0*(-2.D0*RR13*RR23*(CB2*Lam8 + Lam7*SB2) - 1.D0*RR21*(-2.D0*CB*(Lam4 + Lam5)*RR12*SB + 2.D0*RR11*(CB2*Lam3 +&
  & Lam1*SB2)) - 1.D0*RR22*(-2.D0*CB*(Lam4 + Lam5)*RR11*SB + 2.D0*RR12*(CB2*Lam2 + Lam3*SB2)))
CS1S1S3S3f2211 = 0.5D0*(-1.D0*RR22*(2.D0*CB*(Lam4 + Lam5)*RR21*SB + 2.D0*RR22*(CB2*Lam3 + Lam2*SB2)) - 1.D0*RR21*(2.D0*CB*(Lam4 +&
  & Lam5)*RR22*SB + 2.D0*RR21*(CB2*Lam1 + Lam3*SB2)) - 2.D0*(CB2*Lam7 + Lam8*SB2)*DBLE(RR23**INT(2.D0)))
CS1S1S3S3f2212 = 0.5D0*(-1.D0*RR22*(2.D0*RR22*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (Lam4 + Lam5)*RR21*(CB2 - 1.D0*SB2)) - 1.D0*RR21*(&
  &2.D0*RR21*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR22*(CB2 - 1.D0*SB2)) - 2.D0*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*DBLE(R&
  &R23**INT(2.D0)))
CS1S1S3S3f2221 = 0.5D0*(-1.D0*RR22*(2.D0*RR22*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (Lam4 + Lam5)*RR21*(CB2 - 1.D0*SB2)) - 1.D0*RR21*(&
  &2.D0*RR21*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR22*(CB2 - 1.D0*SB2)) - 2.D0*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*DBLE(R&
  &R23**INT(2.D0)))
CS1S1S3S3f2222 = 0.5D0*(-1.D0*RR21*(-2.D0*CB*(Lam4 + Lam5)*RR22*SB + 2.D0*RR21*(CB2*Lam3 + Lam1*SB2)) - 1.D0*RR22*(-2.D0*CB*(Lam4&
  & + Lam5)*RR21*SB + 2.D0*RR22*(CB2*Lam2 + Lam3*SB2)) - 2.D0*(CB2*Lam8 + Lam7*SB2)*DBLE(RR23**INT(2.D0)))
CS1S1S3S3f2311 = 0.5D0*(-2.D0*RR23*RR33*(CB2*Lam7 + Lam8*SB2) - 1.D0*RR22*(2.D0*CB*(Lam4 + Lam5)*RR31*SB + 2.D0*RR32*(CB2*Lam3 + &
  &Lam2*SB2)) - 1.D0*RR21*(2.D0*CB*(Lam4 + Lam5)*RR32*SB + 2.D0*RR31*(CB2*Lam1 + Lam3*SB2)))
CS1S1S3S3f2312 = 0.5D0*(-2.D0*RR23*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR22*(2.D0*RR32*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR31*(CB2 - 1.D0*SB2)) - 1.D0*RR21*(2.D0*RR31*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR32*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f2321 = 0.5D0*(-2.D0*RR23*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR22*(2.D0*RR32*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR31*(CB2 - 1.D0*SB2)) - 1.D0*RR21*(2.D0*RR31*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR32*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f2322 = 0.5D0*(-2.D0*RR23*RR33*(CB2*Lam8 + Lam7*SB2) - 1.D0*RR21*(-2.D0*CB*(Lam4 + Lam5)*RR32*SB + 2.D0*RR31*(CB2*Lam3 +&
  & Lam1*SB2)) - 1.D0*RR22*(-2.D0*CB*(Lam4 + Lam5)*RR31*SB + 2.D0*RR32*(CB2*Lam2 + Lam3*SB2)))
CS1S1S3S3f3111 = 0.5D0*(-2.D0*RR13*RR33*(CB2*Lam7 + Lam8*SB2) - 1.D0*RR32*(2.D0*CB*(Lam4 + Lam5)*RR11*SB + 2.D0*RR12*(CB2*Lam3 + &
  &Lam2*SB2)) - 1.D0*RR31*(2.D0*CB*(Lam4 + Lam5)*RR12*SB + 2.D0*RR11*(CB2*Lam1 + Lam3*SB2)))
CS1S1S3S3f3112 = 0.5D0*(-2.D0*RR13*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR32*(2.D0*RR12*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR11*(CB2 - 1.D0*SB2)) - 1.D0*RR31*(2.D0*RR11*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR12*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f3121 = 0.5D0*(-2.D0*RR13*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR32*(2.D0*RR12*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR11*(CB2 - 1.D0*SB2)) - 1.D0*RR31*(2.D0*RR11*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR12*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f3122 = 0.5D0*(-2.D0*RR13*RR33*(CB2*Lam8 + Lam7*SB2) - 1.D0*RR31*(-2.D0*CB*(Lam4 + Lam5)*RR12*SB + 2.D0*RR11*(CB2*Lam3 +&
  & Lam1*SB2)) - 1.D0*RR32*(-2.D0*CB*(Lam4 + Lam5)*RR11*SB + 2.D0*RR12*(CB2*Lam2 + Lam3*SB2)))
CS1S1S3S3f3211 = 0.5D0*(-2.D0*RR23*RR33*(CB2*Lam7 + Lam8*SB2) - 1.D0*RR32*(2.D0*CB*(Lam4 + Lam5)*RR21*SB + 2.D0*RR22*(CB2*Lam3 + &
  &Lam2*SB2)) - 1.D0*RR31*(2.D0*CB*(Lam4 + Lam5)*RR22*SB + 2.D0*RR21*(CB2*Lam1 + Lam3*SB2)))
CS1S1S3S3f3212 = 0.5D0*(-2.D0*RR23*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR32*(2.D0*RR22*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR21*(CB2 - 1.D0*SB2)) - 1.D0*RR31*(2.D0*RR21*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR22*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f3221 = 0.5D0*(-2.D0*RR23*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR32*(2.D0*RR22*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR21*(CB2 - 1.D0*SB2)) - 1.D0*RR31*(2.D0*RR21*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR22*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f3222 = 0.5D0*(-2.D0*RR23*RR33*(CB2*Lam8 + Lam7*SB2) - 1.D0*RR31*(-2.D0*CB*(Lam4 + Lam5)*RR22*SB + 2.D0*RR21*(CB2*Lam3 +&
  & Lam1*SB2)) - 1.D0*RR32*(-2.D0*CB*(Lam4 + Lam5)*RR21*SB + 2.D0*RR22*(CB2*Lam2 + Lam3*SB2)))
CS1S1S3S3f3311 = 0.5D0*(-1.D0*RR32*(2.D0*CB*(Lam4 + Lam5)*RR31*SB + 2.D0*RR32*(CB2*Lam3 + Lam2*SB2)) - 1.D0*RR31*(2.D0*CB*(Lam4 +&
  & Lam5)*RR32*SB + 2.D0*RR31*(CB2*Lam1 + Lam3*SB2)) - 2.D0*(CB2*Lam7 + Lam8*SB2)*DBLE(RR33**INT(2.D0)))
CS1S1S3S3f3312 = 0.5D0*(-1.D0*RR32*(2.D0*RR32*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (Lam4 + Lam5)*RR31*(CB2 - 1.D0*SB2)) - 1.D0*RR31*(&
  &2.D0*RR31*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR32*(CB2 - 1.D0*SB2)) - 2.D0*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*DBLE(R&
  &R33**INT(2.D0)))
CS1S1S3S3f3321 = 0.5D0*(-1.D0*RR32*(2.D0*RR32*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (Lam4 + Lam5)*RR31*(CB2 - 1.D0*SB2)) - 1.D0*RR31*(&
  &2.D0*RR31*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR32*(CB2 - 1.D0*SB2)) - 2.D0*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*DBLE(R&
  &R33**INT(2.D0)))
CS1S1S3S3f3322 = 0.5D0*(-1.D0*RR31*(-2.D0*CB*(Lam4 + Lam5)*RR32*SB + 2.D0*RR31*(CB2*Lam3 + Lam1*SB2)) - 1.D0*RR32*(-2.D0*CB*(Lam4&
  & + Lam5)*RR31*SB + 2.D0*RR32*(CB2*Lam2 + Lam3*SB2)) - 2.D0*(CB2*Lam8 + Lam7*SB2)*DBLE(RR33**INT(2.D0)))

        call ltini
            call setmudim(1D0)
            call setlambda(1D0)
            call setdelta(0D0)
            IRLambda = getlambda()

            ! Set the UV scale 
            UVDelta = getdelta()

            ! Perform the check for UV divergence for all schemes
            do n = 1, maxNumberSchemes, 1
                ! Reset the checking arrays for all schemes
                isUVDivergentList(n) = .false.
                isUVDivergent = .false.

                ! Info header
                print *, " "
                print *, "==========="
                print *, "UV divergence of scheme ", n
                print *, "==========="

                ! Calculate the amplitude for several different Delta values
                ! Schemes 1, 2, and 9 are without tadpoles
                if ((n == 1) .OR. (n == 2) .OR. (n == 9)) then
                    UVChecks(1) = DBLE(A0toDDBarVC()) + A0toDDBarCT(n)
                else
                    UVChecks(1) = DBLE(A0toDDBarVC() + A0toDDBarTad()) + &
                                & A0toDDBarCT(n)
                end if
                do m = 2, 10, 1
                    call clearcache
                    call setdelta(DBLE(m-1))
                    UVDelta = getdelta()
                    ! Schemes 1, 2, and 9 are without tadpoles
                    if ((n == 1) .OR. (n == 2) .OR. (n == 9)) then
                        UVChecks(m) = ( DBLE(A0toDDBarVC()) + A0toDDBarCT(n) - &
                                    & UVChecks(1) )/UVChecks(1)
                    else
                        UVChecks(m) = ( DBLE(A0toDDBarVC() + A0toDDBarTad()) + &
                                    & A0toDDBarCT(n) - UVChecks(1) &
                                    & )/UVChecks(1)
                    end if
                    print *, "Delta: ", UVDelta, ", Difference from Delta=0: ", UVChecks(m)
                    if (abs(UVChecks(m)) > UVDivergenceThreshold) then
                        isUVDivergent = .true.
                        isUVDivergentList(n) = .true.
                    end if
                end do
                call clearcache
                call setdelta(0D0)
                UVDelta = getdelta()
                print *, "==========="

                ! Exception handling: if UV divergences are found, inform the user and ask if the program shall be terminated
                if (isUVDivergent) then
                    write (*, '(/,A)') "WARNING: potential UV divergence found! Please check the vertex corrections and &
                        &the counterterm!"
                    do
                        print *, ">>> Do you want to continue with the evaluation of the program? [y/n]"
                        read (*,*) isUVDivergentContinue
                        if (isUVDivergentContinue == 'n') then
                            print *, "Termination requested by user. 2HDMCalc will be terminated now."
                            stop
                        else if (isUVDivergentContinue == 'y') then
                            exit
                        else
                            print *, "Invalid character. Enter y or n."
                        end if
                    end do
                end if
            end do
        call ltexi

        ! Print the summary of the checks for UV divergence
        print *, "==========="
        print *, "Results of the check on UV divergence"
        print *, "==========="
        do n = 1, maxNumberSchemes, 1
            if (isUVDivergentList(n)) then
                print *, "Scheme ", n, ": potentially UV-divergent."
            else
                print *, "Scheme ", n, ": most likely UV-finite."
            end if
        end do
        print *, "==========="

    end if

    ! Perform the check for IR divergence
    if (arguments(2) == 1) then
        print *, "Checking for IR divergences ..."

        ! Use this hack to "fill up" the string to the maximum length with whitespace characters so that it can be passed to the subroutine call
        fileName = fileName // ' '

        ! Get all parameters
        call getParameters(fileName)

        ! Set the 2HDM parameters according to the first point in phase-space (random choice)
        MH2 = MH2List(1)
        MH3 = MH3List(1)
        MA0 = MA0List(1)
        MHp = MHpList(1)
        alpha1 = alpha1List(1)
        alpha2 = alpha2List(1)
        alpha3 = alpha3List(1)
        beta = betaList(1)
        CA1 = CA1List(1)
        CA2 = CA2List(1)
        CA3 = CA3List(1)
        SA1 = SA1List(1)
        SA2 = SA2List(1)
        SA3 = SA3List(1)
        CB = CBList(1)
        SB = SBList(1)
        TB = TBList(1)
        YukS1Lep1 = YukS1Lep1List(1)
        YukS1Lep2 = YukS1Lep2List(1)
        YukS1Lep3 = YukS1Lep3List(1)
        YukS2Lep1 = YukS2Lep1List(1)
        YukS2Lep2 = YukS2Lep2List(1)
        YukS3Lep1 = YukS3Lep1List(1)
        YukS3Lep2 = YukS3Lep2List(1)
        YukS1Quark1 = YukS1Quark1List(1)
        YukS1Quark2 = YukS1Quark2List(1)
        YukS1Quark3 = YukS1Quark3List(1)
        YukS2Quark1 = YukS2Quark1List(1)
        YukS2Quark2 = YukS2Quark2List(1)
        YukS3Quark1 = YukS3Quark1List(1)
        YukS3Quark2 = YukS3Quark2List(1)
        m12squared = m12squaredList(1)
        vS = vSList(1)

        ! Calculate the square of the 2HDM input parameters
        MH22 = MH2**2
        MH32 = MH3**2
        MA02 = MA0**2
        MHp2 = MHp**2
        CA12 = CA1**2
        CA22 = CA2**2
        CA32 = CA3**2
        SA12 = SA1**2
        SA22 = SA2**2
        SA32 = SA3**2
        TB2 = TB**2
        SB2 = SB**2
        CB2 = CB**2

        ! Set the scalar couplings
RR11 = CA1*CA2
RR12 = CA2*SA1
RR13 = SA2
RR21 = -1.D0*CA3*SA1 - 1.D0*CA1*SA2*SA3
RR22 = CA1*CA3 - 1.D0*SA1*SA2*SA3
RR23 = CA2*SA3
RR31 = -1.D0*CA1*CA3*SA2 + SA1*SA3
RR32 = -1.D0*CA3*SA1*SA2 - 1.D0*CA1*SA3
RR33 = CA2*CA3

Lam1 = (0.25D0*EL2*((-1.D0*m12squared*SB)/CB + MH12*DBLE(RR11**INT(2.D0)) + MH22*DBLE(RR21**INT(2.D0)) + MH32*DBLE(RR31**INT(2.D0&
  &))))/(CB2*MW2*SW2)
Lam2 = (0.25D0*EL2*((-1.D0*CB*m12squared)/SB + MH12*DBLE(RR12**INT(2.D0)) + MH22*DBLE(RR22**INT(2.D0)) + MH32*DBLE(RR32**INT(2.D0&
  &))))/(MW2*SB2*SW2)
Lam3 = (0.25D0*EL2*(2.D0*MHp2 - (1.D0*m12squared)/(CB*SB) + (MH12*RR11*RR12 + MH22*RR21*RR22 + MH32*RR31*RR32)/(CB*SB)))/(MW2*SW2&
  &)
Lam4 = (0.25D0*EL2*(MA02 - 2.D0*MHp2 + m12squared/(CB*SB)))/(MW2*SW2)
Lam5 = (0.25D0*EL2*(-1.D0*MA02 + m12squared/(CB*SB)))/(MW2*SW2)
Lam6 = (MH12*DBLE(RR13**INT(2.D0)) + MH22*DBLE(RR23**INT(2.D0)) + MH32*DBLE(RR33**INT(2.D0)))*DBLE(vS**INT(-2.D0))
Lam7 = (0.5D0*EL*(MH12*RR11*RR13 + MH22*RR21*RR23 + MH32*RR31*RR33))/(CB*MW*SW*vS)
Lam8 = (0.5D0*EL*(MH12*RR12*RR13 + MH22*RR22*RR23 + MH32*RR32*RR33))/(MW*SB*SW*vS)

CS1S1S1f111 = -1.D0*RR13*(Lam7*RR11*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR12*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR11*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f112 = -1.D0*RR13*(Lam7*RR11*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR12*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR11*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f113 = -1.D0*RR13*(Lam7*RR11*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR12*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR11*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f121 = -1.D0*RR13*(Lam7*RR21*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR22*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR21*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f122 = -1.D0*RR13*(Lam7*RR21*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR22*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR21*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f123 = -1.D0*RR13*(Lam7*RR21*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR22*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR21*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f131 = -1.D0*RR13*(Lam7*RR31*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR32*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR31*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f132 = -1.D0*RR13*(Lam7*RR31*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR32*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR31*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f133 = -1.D0*RR13*(Lam7*RR31*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR32*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR31*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f211 = -1.D0*RR23*(Lam7*RR11*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR12*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR11*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f212 = -1.D0*RR23*(Lam7*RR11*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR12*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR11*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f213 = -1.D0*RR23*(Lam7*RR11*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR12*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR11*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f221 = -1.D0*RR23*(Lam7*RR21*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR22*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR21*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f222 = -1.D0*RR23*(Lam7*RR21*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR22*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR21*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f223 = -1.D0*RR23*(Lam7*RR21*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR22*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR21*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f231 = -1.D0*RR23*(Lam7*RR31*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR32*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR31*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f232 = -1.D0*RR23*(Lam7*RR31*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR32*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR31*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f233 = -1.D0*RR23*(Lam7*RR31*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR32*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR31*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f311 = -1.D0*RR33*(Lam7*RR11*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR12*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR11*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f312 = -1.D0*RR33*(Lam7*RR11*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR12*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR11*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f313 = -1.D0*RR33*(Lam7*RR11*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR12*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR11*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f321 = -1.D0*RR33*(Lam7*RR21*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR22*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR21*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f322 = -1.D0*RR33*(Lam7*RR21*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR22*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR21*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f323 = -1.D0*RR33*(Lam7*RR21*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR22*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR21*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f331 = -1.D0*RR33*(Lam7*RR31*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR32*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR31*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f332 = -1.D0*RR33*(Lam7*RR31*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR32*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR31*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f333 = -1.D0*RR33*(Lam7*RR31*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR32*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR31*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))

CS2S2S1f111 = -1.D0*CB*(Lam5*SB*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + CB*((2.D0*CB*Lam1*MW*RR11*SW)/EL + (2.D0*(L&
  &am3 + Lam4 - 1.D0*Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*SB*(CB*Lam5*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/&
  &EL) + SB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR11*SW)/EL + (2.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS2S2S1f112 = -1.D0*CB*(Lam5*SB*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + CB*((2.D0*CB*Lam1*MW*RR21*SW)/EL + (2.D0*(L&
  &am3 + Lam4 - 1.D0*Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*SB*(CB*Lam5*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/&
  &EL) + SB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR21*SW)/EL + (2.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS2S2S1f113 = -1.D0*CB*(Lam5*SB*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + CB*((2.D0*CB*Lam1*MW*RR31*SW)/EL + (2.D0*(L&
  &am3 + Lam4 - 1.D0*Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*SB*(CB*Lam5*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/&
  &EL) + SB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR31*SW)/EL + (2.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS2S2S1f121 = -1.D0*CB*(CB*Lam5*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) - 1.D0*SB*((2.D0*CB*Lam1*MW*RR11*SW)/EL + (2.&
  &D0*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*SB*(-1.D0*Lam5*SB*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*R&
  &R11*SB*SW)/EL) + CB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR11*SW)/EL + (2.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS2S2S1f122 = -1.D0*CB*(CB*Lam5*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) - 1.D0*SB*((2.D0*CB*Lam1*MW*RR21*SW)/EL + (2.&
  &D0*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*SB*(-1.D0*Lam5*SB*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*R&
  &R21*SB*SW)/EL) + CB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR21*SW)/EL + (2.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS2S2S1f123 = -1.D0*CB*(CB*Lam5*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) - 1.D0*SB*((2.D0*CB*Lam1*MW*RR31*SW)/EL + (2.&
  &D0*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*SB*(-1.D0*Lam5*SB*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*R&
  &R31*SB*SW)/EL) + CB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR31*SW)/EL + (2.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS2S2S1f211 = SB*(Lam5*SB*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + CB*((2.D0*CB*Lam1*MW*RR11*SW)/EL + (2.D0*(Lam3 + &
  &Lam4 - 1.D0*Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*CB*(CB*Lam5*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + &
  &SB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR11*SW)/EL + (2.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS2S2S1f212 = SB*(Lam5*SB*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + CB*((2.D0*CB*Lam1*MW*RR21*SW)/EL + (2.D0*(Lam3 + &
  &Lam4 - 1.D0*Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*CB*(CB*Lam5*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + &
  &SB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR21*SW)/EL + (2.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS2S2S1f213 = SB*(Lam5*SB*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + CB*((2.D0*CB*Lam1*MW*RR31*SW)/EL + (2.D0*(Lam3 + &
  &Lam4 - 1.D0*Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*CB*(CB*Lam5*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + &
  &SB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR31*SW)/EL + (2.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS2S2S1f221 = SB*(CB*Lam5*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) - 1.D0*SB*((2.D0*CB*Lam1*MW*RR11*SW)/EL + (2.D0*(La&
  &m3 + Lam4 - 1.D0*Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*CB*(-1.D0*Lam5*SB*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB&
  &*SW)/EL) + CB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR11*SW)/EL + (2.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS2S2S1f222 = SB*(CB*Lam5*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) - 1.D0*SB*((2.D0*CB*Lam1*MW*RR21*SW)/EL + (2.D0*(La&
  &m3 + Lam4 - 1.D0*Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*CB*(-1.D0*Lam5*SB*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB&
  &*SW)/EL) + CB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR21*SW)/EL + (2.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS2S2S1f223 = SB*(CB*Lam5*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) - 1.D0*SB*((2.D0*CB*Lam1*MW*RR31*SW)/EL + (2.D0*(La&
  &m3 + Lam4 - 1.D0*Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*CB*(-1.D0*Lam5*SB*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB&
  &*SW)/EL) + CB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR31*SW)/EL + (2.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))

CS1S3S3f111 = 0.5D0*(-1.D0*RR12*(CB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + SB*((2.D0*CB2*(Lam4 + La&
  &m5)*MW*SW)/EL + (4.D0*Lam2*MW*SB2*SW)/EL)) - 1.D0*RR11*(SB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) &
  &+ CB*((4.D0*CB2*Lam1*MW*SW)/EL + (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR13*(CB2*Lam7 + Lam8*SB2)*vs)
CS1S3S3f112 = 0.5D0*(-1.D0*RR12*(SB*((4.D0*CB*Lam2*MW*SB*SW)/EL - (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + CB*((2.D0*CB2*(Lam4 + La&
  &m5)*MW*SW)/EL - (4.D0*Lam3*MW*SB2*SW)/EL)) - 1.D0*RR11*(CB*((-4.D0*CB*Lam1*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL)&
  & + SB*((4.D0*CB2*Lam3*MW*SW)/EL - (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR13*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*vs)
CS1S3S3f121 = 0.5D0*(-1.D0*RR12*(-1.D0*SB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + CB*((2.D0*CB2*(Lam&
  &4 + Lam5)*MW*SW)/EL + (4.D0*Lam2*MW*SB2*SW)/EL)) - 1.D0*RR11*(CB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW&
  &)/EL) - 1.D0*SB*((4.D0*CB2*Lam1*MW*SW)/EL + (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR13*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*vs&
  &)
CS1S3S3f122 = 0.5D0*(-1.D0*RR12*(CB*((4.D0*CB*Lam2*MW*SB*SW)/EL - (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) - 1.D0*SB*((2.D0*CB2*(Lam4&
  & + Lam5)*MW*SW)/EL - (4.D0*Lam3*MW*SB2*SW)/EL)) - 1.D0*RR11*(-1.D0*SB*((-4.D0*CB*Lam1*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW&
  &*SB*SW)/EL) + CB*((4.D0*CB2*Lam3*MW*SW)/EL - (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR13*(CB2*Lam8 + Lam7*SB2)*vs)
CS1S3S3f211 = 0.5D0*(-1.D0*RR22*(CB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + SB*((2.D0*CB2*(Lam4 + La&
  &m5)*MW*SW)/EL + (4.D0*Lam2*MW*SB2*SW)/EL)) - 1.D0*RR21*(SB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) &
  &+ CB*((4.D0*CB2*Lam1*MW*SW)/EL + (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR23*(CB2*Lam7 + Lam8*SB2)*vs)
CS1S3S3f212 = 0.5D0*(-1.D0*RR22*(SB*((4.D0*CB*Lam2*MW*SB*SW)/EL - (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + CB*((2.D0*CB2*(Lam4 + La&
  &m5)*MW*SW)/EL - (4.D0*Lam3*MW*SB2*SW)/EL)) - 1.D0*RR21*(CB*((-4.D0*CB*Lam1*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL)&
  & + SB*((4.D0*CB2*Lam3*MW*SW)/EL - (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR23*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*vs)
CS1S3S3f221 = 0.5D0*(-1.D0*RR22*(-1.D0*SB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + CB*((2.D0*CB2*(Lam&
  &4 + Lam5)*MW*SW)/EL + (4.D0*Lam2*MW*SB2*SW)/EL)) - 1.D0*RR21*(CB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW&
  &)/EL) - 1.D0*SB*((4.D0*CB2*Lam1*MW*SW)/EL + (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR23*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*vs&
  &)
CS1S3S3f222 = 0.5D0*(-1.D0*RR22*(CB*((4.D0*CB*Lam2*MW*SB*SW)/EL - (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) - 1.D0*SB*((2.D0*CB2*(Lam4&
  & + Lam5)*MW*SW)/EL - (4.D0*Lam3*MW*SB2*SW)/EL)) - 1.D0*RR21*(-1.D0*SB*((-4.D0*CB*Lam1*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW&
  &*SB*SW)/EL) + CB*((4.D0*CB2*Lam3*MW*SW)/EL - (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR23*(CB2*Lam8 + Lam7*SB2)*vs)
CS1S3S3f311 = 0.5D0*(-1.D0*RR32*(CB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + SB*((2.D0*CB2*(Lam4 + La&
  &m5)*MW*SW)/EL + (4.D0*Lam2*MW*SB2*SW)/EL)) - 1.D0*RR31*(SB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) &
  &+ CB*((4.D0*CB2*Lam1*MW*SW)/EL + (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR33*(CB2*Lam7 + Lam8*SB2)*vs)
CS1S3S3f312 = 0.5D0*(-1.D0*RR32*(SB*((4.D0*CB*Lam2*MW*SB*SW)/EL - (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + CB*((2.D0*CB2*(Lam4 + La&
  &m5)*MW*SW)/EL - (4.D0*Lam3*MW*SB2*SW)/EL)) - 1.D0*RR31*(CB*((-4.D0*CB*Lam1*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL)&
  & + SB*((4.D0*CB2*Lam3*MW*SW)/EL - (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*vs)
CS1S3S3f321 = 0.5D0*(-1.D0*RR32*(-1.D0*SB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + CB*((2.D0*CB2*(Lam&
  &4 + Lam5)*MW*SW)/EL + (4.D0*Lam2*MW*SB2*SW)/EL)) - 1.D0*RR31*(CB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW&
  &)/EL) - 1.D0*SB*((4.D0*CB2*Lam1*MW*SW)/EL + (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*vs&
  &)
CS1S3S3f322 = 0.5D0*(-1.D0*RR32*(CB*((4.D0*CB*Lam2*MW*SB*SW)/EL - (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) - 1.D0*SB*((2.D0*CB2*(Lam4&
  & + Lam5)*MW*SW)/EL - (4.D0*Lam3*MW*SB2*SW)/EL)) - 1.D0*RR31*(-1.D0*SB*((-4.D0*CB*Lam1*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW&
  &*SB*SW)/EL) + CB*((4.D0*CB2*Lam3*MW*SW)/EL - (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR33*(CB2*Lam8 + Lam7*SB2)*vs)

CS1S1S1S1f1111 = -1.D0*RR13*(2.D0*Lam7*RR13*DBLE(RR11**INT(2.D0)) + 2.D0*Lam8*RR13*DBLE(RR12**INT(2.D0)) + RR13*(Lam7*DBLE(RR11**&
  &INT(2.D0)) + Lam8*DBLE(RR12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*DBLE(R&
  &R12**INT(2.D0)) + 2.D0*Lam7*RR11*DBLE(RR13**INT(2.D0)) + RR11*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR1&
  &2**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*DBLE(RR11**INT(2.D0)) + 2.D0*Lam8*RR1&
  &2*DBLE(RR13**INT(2.D0)) + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13*&
  &*INT(2.D0))))
CS1S1S1S1f1112 = -1.D0*RR13*(Lam7*RR11*(RR13*RR21 + RR11*RR23) + Lam8*RR12*(RR13*RR22 + RR12*RR23) + RR13*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR21 + RR11*RR22) + Lam7*RR13*(RR13*RR21 + RR11&
  &*RR23) + RR11*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR13*(RR13*RR22 + RR12*RR23) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f1113 = -1.D0*RR13*(Lam7*RR11*(RR13*RR31 + RR11*RR33) + Lam8*RR12*(RR13*RR32 + RR12*RR33) + RR13*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR31 + RR11*RR32) + Lam7*RR13*(RR13*RR31 + RR11&
  &*RR33) + RR11*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR13*(RR13*RR32 + RR12*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f1121 = -1.D0*RR13*(Lam7*RR11*(RR13*RR21 + RR11*RR23) + Lam8*RR12*(RR13*RR22 + RR12*RR23) + RR13*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR21 + RR11*RR22) + Lam7*RR13*(RR13*RR21 + RR11&
  &*RR23) + RR11*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR13*(RR13*RR22 + RR12*RR23) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f1122 = -1.D0*RR13*(2.D0*Lam7*RR11*RR21*RR23 + 2.D0*Lam8*RR12*RR22*RR23 + RR13*(Lam7*DBLE(RR21**INT(2.D0)) + Lam8*DBLE(R&
  &R22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*RR21*RR22 + 2.D0*Lam7*RR13*RR2&
  &1*RR23 + RR11*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR22**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1&
  &.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR21*RR22 + 2.D0*Lam8*RR13*RR22*RR23 + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS1S1S1S1f1123 = -1.D0*RR13*(Lam7*RR11*(RR23*RR31 + RR21*RR33) + Lam8*RR12*(RR23*RR32 + RR22*RR33) + RR13*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*(RR22*RR31 + RR21*RR32) + Lam7*RR13*(RR23*RR31 + RR21&
  &*RR33) + RR11*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR13*(RR23*RR32 + RR22*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f1131 = -1.D0*RR13*(Lam7*RR11*(RR13*RR31 + RR11*RR33) + Lam8*RR12*(RR13*RR32 + RR12*RR33) + RR13*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR31 + RR11*RR32) + Lam7*RR13*(RR13*RR31 + RR11&
  &*RR33) + RR11*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR13*(RR13*RR32 + RR12*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f1132 = -1.D0*RR13*(Lam7*RR11*(RR23*RR31 + RR21*RR33) + Lam8*RR12*(RR23*RR32 + RR22*RR33) + RR13*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*(RR22*RR31 + RR21*RR32) + Lam7*RR13*(RR23*RR31 + RR21&
  &*RR33) + RR11*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR13*(RR23*RR32 + RR22*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f1133 = -1.D0*RR13*(2.D0*Lam7*RR11*RR31*RR33 + 2.D0*Lam8*RR12*RR32*RR33 + RR13*(Lam7*DBLE(RR31**INT(2.D0)) + Lam8*DBLE(R&
  &R32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*RR31*RR32 + 2.D0*Lam7*RR13*RR3&
  &1*RR33 + RR11*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR32**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1&
  &.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR31*RR32 + 2.D0*Lam8*RR13*RR32*RR33 + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS1S1S1S1f1211 = -1.D0*RR13*(2.D0*Lam7*RR11*RR13*RR21 + 2.D0*Lam8*RR12*RR13*RR22 + RR23*(Lam7*DBLE(RR11**INT(2.D0)) + Lam8*DBLE(R&
  &R12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR22 + 2.D0*Lam7*RR11*RR1&
  &3*RR23 + RR21*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR12**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1&
  &.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR21 + 2.D0*Lam8*RR12*RR13*RR23 + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS1S1S1S1f1212 = -1.D0*RR13*(Lam7*RR21*(RR13*RR21 + RR11*RR23) + Lam8*RR22*(RR13*RR22 + RR12*RR23) + RR23*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR21 + RR11*RR22) + Lam7*RR23*(RR13*RR21 + RR11&
  &*RR23) + RR21*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR23*(RR13*RR22 + RR12*RR23) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f1213 = -1.D0*RR13*(Lam7*RR21*(RR13*RR31 + RR11*RR33) + Lam8*RR22*(RR13*RR32 + RR12*RR33) + RR23*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR31 + RR11*RR32) + Lam7*RR23*(RR13*RR31 + RR11&
  &*RR33) + RR21*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR23*(RR13*RR32 + RR12*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f1221 = -1.D0*RR13*(Lam7*RR21*(RR13*RR21 + RR11*RR23) + Lam8*RR22*(RR13*RR22 + RR12*RR23) + RR23*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR21 + RR11*RR22) + Lam7*RR23*(RR13*RR21 + RR11&
  &*RR23) + RR21*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR23*(RR13*RR22 + RR12*RR23) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f1222 = -1.D0*RR13*(2.D0*Lam7*RR23*DBLE(RR21**INT(2.D0)) + 2.D0*Lam8*RR23*DBLE(RR22**INT(2.D0)) + RR23*(Lam7*DBLE(RR21**&
  &INT(2.D0)) + Lam8*DBLE(RR22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*DBLE(R&
  &R22**INT(2.D0)) + 2.D0*Lam7*RR21*DBLE(RR23**INT(2.D0)) + RR21*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR2&
  &2**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR22*DBLE(RR21**INT(2.D0)) + 2.D0*Lam8*RR2&
  &2*DBLE(RR23**INT(2.D0)) + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23*&
  &*INT(2.D0))))
CS1S1S1S1f1223 = -1.D0*RR13*(Lam7*RR21*(RR23*RR31 + RR21*RR33) + Lam8*RR22*(RR23*RR32 + RR22*RR33) + RR23*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*(RR22*RR31 + RR21*RR32) + Lam7*RR23*(RR23*RR31 + RR21&
  &*RR33) + RR21*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR23*(RR23*RR32 + RR22*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f1231 = -1.D0*RR13*(Lam7*RR21*(RR13*RR31 + RR11*RR33) + Lam8*RR22*(RR13*RR32 + RR12*RR33) + RR23*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR31 + RR11*RR32) + Lam7*RR23*(RR13*RR31 + RR11&
  &*RR33) + RR21*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR23*(RR13*RR32 + RR12*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f1232 = -1.D0*RR13*(Lam7*RR21*(RR23*RR31 + RR21*RR33) + Lam8*RR22*(RR23*RR32 + RR22*RR33) + RR23*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*(RR22*RR31 + RR21*RR32) + Lam7*RR23*(RR23*RR31 + RR21&
  &*RR33) + RR21*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR23*(RR23*RR32 + RR22*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f1233 = -1.D0*RR13*(2.D0*Lam7*RR21*RR31*RR33 + 2.D0*Lam8*RR22*RR32*RR33 + RR23*(Lam7*DBLE(RR31**INT(2.D0)) + Lam8*DBLE(R&
  &R32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR22*RR31*RR32 + 2.D0*Lam7*RR23*RR3&
  &1*RR33 + RR21*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR32**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1&
  &.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR31*RR32 + 2.D0*Lam8*RR23*RR32*RR33 + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS1S1S1S1f1311 = -1.D0*RR13*(2.D0*Lam7*RR11*RR13*RR31 + 2.D0*Lam8*RR12*RR13*RR32 + RR33*(Lam7*DBLE(RR11**INT(2.D0)) + Lam8*DBLE(R&
  &R12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR32 + 2.D0*Lam7*RR11*RR1&
  &3*RR33 + RR31*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR12**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1&
  &.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR31 + 2.D0*Lam8*RR12*RR13*RR33 + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS1S1S1S1f1312 = -1.D0*RR11*((3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)*RR31 + (Lam3 + Lam4 + Lam5)*&
  &(RR12*RR21 + RR11*RR22)*RR32 + Lam7*(RR13*RR21 + RR11*RR23)*RR33) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*(RR12*RR21 + RR11*RR22)*RR&
  &31 + ((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23)*RR32 + Lam8*(RR13*RR22 + RR12*RR23)*RR33) - 1.D0*&
  &RR13*(Lam7*(RR13*RR21 + RR11*RR23)*RR31 + Lam8*(RR13*RR22 + RR12*RR23)*RR32 + (Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23)*RR33)
CS1S1S1S1f1313 = -1.D0*RR13*(Lam7*RR31*(RR13*RR31 + RR11*RR33) + Lam8*RR32*(RR13*RR32 + RR12*RR33) + RR33*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR32*(RR12*RR31 + RR11*RR32) + Lam7*RR33*(RR13*RR31 + RR11&
  &*RR33) + RR31*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR33*(RR13*RR32 + RR12*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f1321 = -1.D0*RR11*((3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)*RR31 + (Lam3 + Lam4 + Lam5)*&
  &(RR12*RR21 + RR11*RR22)*RR32 + Lam7*(RR13*RR21 + RR11*RR23)*RR33) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*(RR12*RR21 + RR11*RR22)*RR&
  &31 + ((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23)*RR32 + Lam8*(RR13*RR22 + RR12*RR23)*RR33) - 1.D0*&
  &RR13*(Lam7*(RR13*RR21 + RR11*RR23)*RR31 + Lam8*(RR13*RR22 + RR12*RR23)*RR32 + (Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23)*RR33)
CS1S1S1S1f1322 = -1.D0*RR13*(2.D0*Lam7*RR21*RR23*RR31 + 2.D0*Lam8*RR22*RR23*RR32 + RR33*(Lam7*DBLE(RR21**INT(2.D0)) + Lam8*DBLE(R&
  &R22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR22*RR32 + 2.D0*Lam7*RR21*RR2&
  &3*RR33 + RR31*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR22**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1&
  &.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR22*RR31 + 2.D0*Lam8*RR22*RR23*RR33 + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS1S1S1S1f1323 = -1.D0*RR13*(Lam7*RR31*(RR23*RR31 + RR21*RR33) + Lam8*RR32*(RR23*RR32 + RR22*RR33) + RR33*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR32*(RR22*RR31 + RR21*RR32) + Lam7*RR33*(RR23*RR31 + RR21&
  &*RR33) + RR31*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR33*(RR23*RR32 + RR22*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f1331 = -1.D0*RR13*(Lam7*RR31*(RR13*RR31 + RR11*RR33) + Lam8*RR32*(RR13*RR32 + RR12*RR33) + RR33*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR32*(RR12*RR31 + RR11*RR32) + Lam7*RR33*(RR13*RR31 + RR11&
  &*RR33) + RR31*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR33*(RR13*RR32 + RR12*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f1332 = -1.D0*RR13*(Lam7*RR31*(RR23*RR31 + RR21*RR33) + Lam8*RR32*(RR23*RR32 + RR22*RR33) + RR33*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR32*(RR22*RR31 + RR21*RR32) + Lam7*RR33*(RR23*RR31 + RR21&
  &*RR33) + RR31*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR33*(RR23*RR32 + RR22*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f1333 = -1.D0*RR13*(2.D0*Lam7*RR33*DBLE(RR31**INT(2.D0)) + 2.D0*Lam8*RR33*DBLE(RR32**INT(2.D0)) + RR33*(Lam7*DBLE(RR31**&
  &INT(2.D0)) + Lam8*DBLE(RR32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR31*DBLE(R&
  &R32**INT(2.D0)) + 2.D0*Lam7*RR31*DBLE(RR33**INT(2.D0)) + RR31*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR3&
  &2**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR32*DBLE(RR31**INT(2.D0)) + 2.D0*Lam8*RR3&
  &2*DBLE(RR33**INT(2.D0)) + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33*&
  &*INT(2.D0))))
CS1S1S1S1f2111 = -1.D0*RR23*(2.D0*Lam7*RR13*DBLE(RR11**INT(2.D0)) + 2.D0*Lam8*RR13*DBLE(RR12**INT(2.D0)) + RR13*(Lam7*DBLE(RR11**&
  &INT(2.D0)) + Lam8*DBLE(RR12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*DBLE(R&
  &R12**INT(2.D0)) + 2.D0*Lam7*RR11*DBLE(RR13**INT(2.D0)) + RR11*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR1&
  &2**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*DBLE(RR11**INT(2.D0)) + 2.D0*Lam8*RR1&
  &2*DBLE(RR13**INT(2.D0)) + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13*&
  &*INT(2.D0))))
CS1S1S1S1f2112 = -1.D0*RR23*(Lam7*RR11*(RR13*RR21 + RR11*RR23) + Lam8*RR12*(RR13*RR22 + RR12*RR23) + RR13*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR21 + RR11*RR22) + Lam7*RR13*(RR13*RR21 + RR11&
  &*RR23) + RR11*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR13*(RR13*RR22 + RR12*RR23) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f2113 = -1.D0*RR23*(Lam7*RR11*(RR13*RR31 + RR11*RR33) + Lam8*RR12*(RR13*RR32 + RR12*RR33) + RR13*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR31 + RR11*RR32) + Lam7*RR13*(RR13*RR31 + RR11&
  &*RR33) + RR11*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR13*(RR13*RR32 + RR12*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f2121 = -1.D0*RR23*(Lam7*RR11*(RR13*RR21 + RR11*RR23) + Lam8*RR12*(RR13*RR22 + RR12*RR23) + RR13*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR21 + RR11*RR22) + Lam7*RR13*(RR13*RR21 + RR11&
  &*RR23) + RR11*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR13*(RR13*RR22 + RR12*RR23) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f2122 = -1.D0*RR23*(2.D0*Lam7*RR11*RR21*RR23 + 2.D0*Lam8*RR12*RR22*RR23 + RR13*(Lam7*DBLE(RR21**INT(2.D0)) + Lam8*DBLE(R&
  &R22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*RR21*RR22 + 2.D0*Lam7*RR13*RR2&
  &1*RR23 + RR11*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR22**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1&
  &.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR21*RR22 + 2.D0*Lam8*RR13*RR22*RR23 + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS1S1S1S1f2123 = -1.D0*RR23*(Lam7*RR11*(RR23*RR31 + RR21*RR33) + Lam8*RR12*(RR23*RR32 + RR22*RR33) + RR13*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*(RR22*RR31 + RR21*RR32) + Lam7*RR13*(RR23*RR31 + RR21&
  &*RR33) + RR11*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR13*(RR23*RR32 + RR22*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f2131 = -1.D0*RR23*(Lam7*RR11*(RR13*RR31 + RR11*RR33) + Lam8*RR12*(RR13*RR32 + RR12*RR33) + RR13*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR31 + RR11*RR32) + Lam7*RR13*(RR13*RR31 + RR11&
  &*RR33) + RR11*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR13*(RR13*RR32 + RR12*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f2132 = -1.D0*RR23*(Lam7*RR11*(RR23*RR31 + RR21*RR33) + Lam8*RR12*(RR23*RR32 + RR22*RR33) + RR13*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*(RR22*RR31 + RR21*RR32) + Lam7*RR13*(RR23*RR31 + RR21&
  &*RR33) + RR11*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR13*(RR23*RR32 + RR22*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f2133 = -1.D0*RR23*(2.D0*Lam7*RR11*RR31*RR33 + 2.D0*Lam8*RR12*RR32*RR33 + RR13*(Lam7*DBLE(RR31**INT(2.D0)) + Lam8*DBLE(R&
  &R32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*RR31*RR32 + 2.D0*Lam7*RR13*RR3&
  &1*RR33 + RR11*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR32**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1&
  &.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR31*RR32 + 2.D0*Lam8*RR13*RR32*RR33 + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS1S1S1S1f2211 = -1.D0*RR23*(2.D0*Lam7*RR11*RR13*RR21 + 2.D0*Lam8*RR12*RR13*RR22 + RR23*(Lam7*DBLE(RR11**INT(2.D0)) + Lam8*DBLE(R&
  &R12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR22 + 2.D0*Lam7*RR11*RR1&
  &3*RR23 + RR21*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR12**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1&
  &.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR21 + 2.D0*Lam8*RR12*RR13*RR23 + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS1S1S1S1f2212 = -1.D0*RR23*(Lam7*RR21*(RR13*RR21 + RR11*RR23) + Lam8*RR22*(RR13*RR22 + RR12*RR23) + RR23*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR21 + RR11*RR22) + Lam7*RR23*(RR13*RR21 + RR11&
  &*RR23) + RR21*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR23*(RR13*RR22 + RR12*RR23) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f2213 = -1.D0*RR23*(Lam7*RR21*(RR13*RR31 + RR11*RR33) + Lam8*RR22*(RR13*RR32 + RR12*RR33) + RR23*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR31 + RR11*RR32) + Lam7*RR23*(RR13*RR31 + RR11&
  &*RR33) + RR21*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR23*(RR13*RR32 + RR12*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f2221 = -1.D0*RR23*(Lam7*RR21*(RR13*RR21 + RR11*RR23) + Lam8*RR22*(RR13*RR22 + RR12*RR23) + RR23*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR21 + RR11*RR22) + Lam7*RR23*(RR13*RR21 + RR11&
  &*RR23) + RR21*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR23*(RR13*RR22 + RR12*RR23) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f2222 = -1.D0*RR23*(2.D0*Lam7*RR23*DBLE(RR21**INT(2.D0)) + 2.D0*Lam8*RR23*DBLE(RR22**INT(2.D0)) + RR23*(Lam7*DBLE(RR21**&
  &INT(2.D0)) + Lam8*DBLE(RR22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*DBLE(R&
  &R22**INT(2.D0)) + 2.D0*Lam7*RR21*DBLE(RR23**INT(2.D0)) + RR21*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR2&
  &2**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR22*DBLE(RR21**INT(2.D0)) + 2.D0*Lam8*RR2&
  &2*DBLE(RR23**INT(2.D0)) + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23*&
  &*INT(2.D0))))
CS1S1S1S1f2223 = -1.D0*RR23*(Lam7*RR21*(RR23*RR31 + RR21*RR33) + Lam8*RR22*(RR23*RR32 + RR22*RR33) + RR23*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*(RR22*RR31 + RR21*RR32) + Lam7*RR23*(RR23*RR31 + RR21&
  &*RR33) + RR21*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR23*(RR23*RR32 + RR22*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f2231 = -1.D0*RR23*(Lam7*RR21*(RR13*RR31 + RR11*RR33) + Lam8*RR22*(RR13*RR32 + RR12*RR33) + RR23*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR31 + RR11*RR32) + Lam7*RR23*(RR13*RR31 + RR11&
  &*RR33) + RR21*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR23*(RR13*RR32 + RR12*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f2232 = -1.D0*RR23*(Lam7*RR21*(RR23*RR31 + RR21*RR33) + Lam8*RR22*(RR23*RR32 + RR22*RR33) + RR23*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*(RR22*RR31 + RR21*RR32) + Lam7*RR23*(RR23*RR31 + RR21&
  &*RR33) + RR21*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR23*(RR23*RR32 + RR22*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f2233 = -1.D0*RR23*(2.D0*Lam7*RR21*RR31*RR33 + 2.D0*Lam8*RR22*RR32*RR33 + RR23*(Lam7*DBLE(RR31**INT(2.D0)) + Lam8*DBLE(R&
  &R32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR22*RR31*RR32 + 2.D0*Lam7*RR23*RR3&
  &1*RR33 + RR21*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR32**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1&
  &.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR31*RR32 + 2.D0*Lam8*RR23*RR32*RR33 + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS1S1S1S1f2311 = -1.D0*RR23*(2.D0*Lam7*RR11*RR13*RR31 + 2.D0*Lam8*RR12*RR13*RR32 + RR33*(Lam7*DBLE(RR11**INT(2.D0)) + Lam8*DBLE(R&
  &R12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR32 + 2.D0*Lam7*RR11*RR1&
  &3*RR33 + RR31*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR12**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1&
  &.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR31 + 2.D0*Lam8*RR12*RR13*RR33 + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS1S1S1S1f2312 = -1.D0*RR21*((3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)*RR31 + (Lam3 + Lam4 + Lam5)*&
  &(RR12*RR21 + RR11*RR22)*RR32 + Lam7*(RR13*RR21 + RR11*RR23)*RR33) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*(RR12*RR21 + RR11*RR22)*RR&
  &31 + ((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23)*RR32 + Lam8*(RR13*RR22 + RR12*RR23)*RR33) - 1.D0*&
  &RR23*(Lam7*(RR13*RR21 + RR11*RR23)*RR31 + Lam8*(RR13*RR22 + RR12*RR23)*RR32 + (Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23)*RR33)
CS1S1S1S1f2313 = -1.D0*RR23*(Lam7*RR31*(RR13*RR31 + RR11*RR33) + Lam8*RR32*(RR13*RR32 + RR12*RR33) + RR33*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR32*(RR12*RR31 + RR11*RR32) + Lam7*RR33*(RR13*RR31 + RR11&
  &*RR33) + RR31*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR33*(RR13*RR32 + RR12*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f2321 = -1.D0*RR21*((3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)*RR31 + (Lam3 + Lam4 + Lam5)*&
  &(RR12*RR21 + RR11*RR22)*RR32 + Lam7*(RR13*RR21 + RR11*RR23)*RR33) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*(RR12*RR21 + RR11*RR22)*RR&
  &31 + ((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23)*RR32 + Lam8*(RR13*RR22 + RR12*RR23)*RR33) - 1.D0*&
  &RR23*(Lam7*(RR13*RR21 + RR11*RR23)*RR31 + Lam8*(RR13*RR22 + RR12*RR23)*RR32 + (Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23)*RR33)
CS1S1S1S1f2322 = -1.D0*RR23*(2.D0*Lam7*RR21*RR23*RR31 + 2.D0*Lam8*RR22*RR23*RR32 + RR33*(Lam7*DBLE(RR21**INT(2.D0)) + Lam8*DBLE(R&
  &R22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR22*RR32 + 2.D0*Lam7*RR21*RR2&
  &3*RR33 + RR31*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR22**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1&
  &.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR22*RR31 + 2.D0*Lam8*RR22*RR23*RR33 + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS1S1S1S1f2323 = -1.D0*RR23*(Lam7*RR31*(RR23*RR31 + RR21*RR33) + Lam8*RR32*(RR23*RR32 + RR22*RR33) + RR33*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR32*(RR22*RR31 + RR21*RR32) + Lam7*RR33*(RR23*RR31 + RR21&
  &*RR33) + RR31*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR33*(RR23*RR32 + RR22*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f2331 = -1.D0*RR23*(Lam7*RR31*(RR13*RR31 + RR11*RR33) + Lam8*RR32*(RR13*RR32 + RR12*RR33) + RR33*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR32*(RR12*RR31 + RR11*RR32) + Lam7*RR33*(RR13*RR31 + RR11&
  &*RR33) + RR31*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR33*(RR13*RR32 + RR12*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f2332 = -1.D0*RR23*(Lam7*RR31*(RR23*RR31 + RR21*RR33) + Lam8*RR32*(RR23*RR32 + RR22*RR33) + RR33*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR32*(RR22*RR31 + RR21*RR32) + Lam7*RR33*(RR23*RR31 + RR21&
  &*RR33) + RR31*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR33*(RR23*RR32 + RR22*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f2333 = -1.D0*RR23*(2.D0*Lam7*RR33*DBLE(RR31**INT(2.D0)) + 2.D0*Lam8*RR33*DBLE(RR32**INT(2.D0)) + RR33*(Lam7*DBLE(RR31**&
  &INT(2.D0)) + Lam8*DBLE(RR32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR31*DBLE(R&
  &R32**INT(2.D0)) + 2.D0*Lam7*RR31*DBLE(RR33**INT(2.D0)) + RR31*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR3&
  &2**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR32*DBLE(RR31**INT(2.D0)) + 2.D0*Lam8*RR3&
  &2*DBLE(RR33**INT(2.D0)) + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33*&
  &*INT(2.D0))))
CS1S1S1S1f3111 = -1.D0*RR33*(2.D0*Lam7*RR13*DBLE(RR11**INT(2.D0)) + 2.D0*Lam8*RR13*DBLE(RR12**INT(2.D0)) + RR13*(Lam7*DBLE(RR11**&
  &INT(2.D0)) + Lam8*DBLE(RR12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*DBLE(R&
  &R12**INT(2.D0)) + 2.D0*Lam7*RR11*DBLE(RR13**INT(2.D0)) + RR11*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR1&
  &2**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*DBLE(RR11**INT(2.D0)) + 2.D0*Lam8*RR1&
  &2*DBLE(RR13**INT(2.D0)) + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13*&
  &*INT(2.D0))))
CS1S1S1S1f3112 = -1.D0*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR21 + RR11*RR22) + Lam7*RR13*(RR13*RR21 + RR11*RR23) + RR11*(3.D0*Lam1*R&
  &R11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23))*RR31 - 1.D0*((Lam3 + Lam4 + Lam5)*RR11*(RR12*RR21 + RR11*RR22) + &
  &Lam8*RR13*(RR13*RR22 + RR12*RR23) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23))*RR32 - 1.D0*&
  &(Lam7*RR11*(RR13*RR21 + RR11*RR23) + Lam8*RR12*(RR13*RR22 + RR12*RR23) + RR13*(Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23))*RR33
CS1S1S1S1f3113 = -1.D0*RR33*(Lam7*RR11*(RR13*RR31 + RR11*RR33) + Lam8*RR12*(RR13*RR32 + RR12*RR33) + RR13*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR31 + RR11*RR32) + Lam7*RR13*(RR13*RR31 + RR11&
  &*RR33) + RR11*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR13*(RR13*RR32 + RR12*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f3121 = -1.D0*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR21 + RR11*RR22) + Lam7*RR13*(RR13*RR21 + RR11*RR23) + RR11*(3.D0*Lam1*R&
  &R11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23))*RR31 - 1.D0*((Lam3 + Lam4 + Lam5)*RR11*(RR12*RR21 + RR11*RR22) + &
  &Lam8*RR13*(RR13*RR22 + RR12*RR23) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23))*RR32 - 1.D0*&
  &(Lam7*RR11*(RR13*RR21 + RR11*RR23) + Lam8*RR12*(RR13*RR22 + RR12*RR23) + RR13*(Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23))*RR33
CS1S1S1S1f3122 = -1.D0*RR33*(2.D0*Lam7*RR11*RR21*RR23 + 2.D0*Lam8*RR12*RR22*RR23 + RR13*(Lam7*DBLE(RR21**INT(2.D0)) + Lam8*DBLE(R&
  &R22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*RR21*RR22 + 2.D0*Lam7*RR13*RR2&
  &1*RR23 + RR11*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR22**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1&
  &.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR21*RR22 + 2.D0*Lam8*RR13*RR22*RR23 + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS1S1S1S1f3123 = -1.D0*RR33*(Lam7*RR11*(RR23*RR31 + RR21*RR33) + Lam8*RR12*(RR23*RR32 + RR22*RR33) + RR13*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR12*(RR22*RR31 + RR21*RR32) + Lam7*RR13*(RR23*RR31 + RR21&
  &*RR33) + RR11*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR13*(RR23*RR32 + RR22*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f3131 = -1.D0*RR33*(Lam7*RR11*(RR13*RR31 + RR11*RR33) + Lam8*RR12*(RR13*RR32 + RR12*RR33) + RR13*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR31 + RR11*RR32) + Lam7*RR13*(RR13*RR31 + RR11&
  &*RR33) + RR11*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR13*(RR13*RR32 + RR12*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f3132 = -1.D0*RR33*(Lam7*RR11*(RR23*RR31 + RR21*RR33) + Lam8*RR12*(RR23*RR32 + RR22*RR33) + RR13*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR12*(RR22*RR31 + RR21*RR32) + Lam7*RR13*(RR23*RR31 + RR21&
  &*RR33) + RR11*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR13*(RR23*RR32 + RR22*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f3133 = -1.D0*RR33*(2.D0*Lam7*RR11*RR31*RR33 + 2.D0*Lam8*RR12*RR32*RR33 + RR13*(Lam7*DBLE(RR31**INT(2.D0)) + Lam8*DBLE(R&
  &R32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*RR31*RR32 + 2.D0*Lam7*RR13*RR3&
  &1*RR33 + RR11*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR32**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1&
  &.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR31*RR32 + 2.D0*Lam8*RR13*RR32*RR33 + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS1S1S1S1f3211 = -1.D0*RR33*(2.D0*Lam7*RR11*RR13*RR21 + 2.D0*Lam8*RR12*RR13*RR22 + RR23*(Lam7*DBLE(RR11**INT(2.D0)) + Lam8*DBLE(R&
  &R12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR22 + 2.D0*Lam7*RR11*RR1&
  &3*RR23 + RR21*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR12**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1&
  &.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR21 + 2.D0*Lam8*RR12*RR13*RR23 + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS1S1S1S1f3212 = -1.D0*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR21 + RR11*RR22) + Lam7*RR23*(RR13*RR21 + RR11*RR23) + RR21*(3.D0*Lam1*R&
  &R11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23))*RR31 - 1.D0*((Lam3 + Lam4 + Lam5)*RR21*(RR12*RR21 + RR11*RR22) + &
  &Lam8*RR23*(RR13*RR22 + RR12*RR23) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23))*RR32 - 1.D0*&
  &(Lam7*RR21*(RR13*RR21 + RR11*RR23) + Lam8*RR22*(RR13*RR22 + RR12*RR23) + RR23*(Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23))*RR33
CS1S1S1S1f3213 = -1.D0*RR33*(Lam7*RR21*(RR13*RR31 + RR11*RR33) + Lam8*RR22*(RR13*RR32 + RR12*RR33) + RR23*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR31 + RR11*RR32) + Lam7*RR23*(RR13*RR31 + RR11&
  &*RR33) + RR21*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR23*(RR13*RR32 + RR12*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f3221 = -1.D0*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR21 + RR11*RR22) + Lam7*RR23*(RR13*RR21 + RR11*RR23) + RR21*(3.D0*Lam1*R&
  &R11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23))*RR31 - 1.D0*((Lam3 + Lam4 + Lam5)*RR21*(RR12*RR21 + RR11*RR22) + &
  &Lam8*RR23*(RR13*RR22 + RR12*RR23) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23))*RR32 - 1.D0*&
  &(Lam7*RR21*(RR13*RR21 + RR11*RR23) + Lam8*RR22*(RR13*RR22 + RR12*RR23) + RR23*(Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23))*RR33
CS1S1S1S1f3222 = -1.D0*RR33*(2.D0*Lam7*RR23*DBLE(RR21**INT(2.D0)) + 2.D0*Lam8*RR23*DBLE(RR22**INT(2.D0)) + RR23*(Lam7*DBLE(RR21**&
  &INT(2.D0)) + Lam8*DBLE(RR22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*DBLE(R&
  &R22**INT(2.D0)) + 2.D0*Lam7*RR21*DBLE(RR23**INT(2.D0)) + RR21*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR2&
  &2**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR22*DBLE(RR21**INT(2.D0)) + 2.D0*Lam8*RR2&
  &2*DBLE(RR23**INT(2.D0)) + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23*&
  &*INT(2.D0))))
CS1S1S1S1f3223 = -1.D0*RR33*(Lam7*RR21*(RR23*RR31 + RR21*RR33) + Lam8*RR22*(RR23*RR32 + RR22*RR33) + RR23*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR22*(RR22*RR31 + RR21*RR32) + Lam7*RR23*(RR23*RR31 + RR21&
  &*RR33) + RR21*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR23*(RR23*RR32 + RR22*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f3231 = -1.D0*RR33*(Lam7*RR21*(RR13*RR31 + RR11*RR33) + Lam8*RR22*(RR13*RR32 + RR12*RR33) + RR23*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR31 + RR11*RR32) + Lam7*RR23*(RR13*RR31 + RR11&
  &*RR33) + RR21*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR23*(RR13*RR32 + RR12*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f3232 = -1.D0*RR33*(Lam7*RR21*(RR23*RR31 + RR21*RR33) + Lam8*RR22*(RR23*RR32 + RR22*RR33) + RR23*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR22*(RR22*RR31 + RR21*RR32) + Lam7*RR23*(RR23*RR31 + RR21&
  &*RR33) + RR21*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR23*(RR23*RR32 + RR22*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f3233 = -1.D0*RR33*(2.D0*Lam7*RR21*RR31*RR33 + 2.D0*Lam8*RR22*RR32*RR33 + RR23*(Lam7*DBLE(RR31**INT(2.D0)) + Lam8*DBLE(R&
  &R32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR22*RR31*RR32 + 2.D0*Lam7*RR23*RR3&
  &1*RR33 + RR21*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR32**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1&
  &.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR31*RR32 + 2.D0*Lam8*RR23*RR32*RR33 + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS1S1S1S1f3311 = -1.D0*RR33*(2.D0*Lam7*RR11*RR13*RR31 + 2.D0*Lam8*RR12*RR13*RR32 + RR33*(Lam7*DBLE(RR11**INT(2.D0)) + Lam8*DBLE(R&
  &R12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR32 + 2.D0*Lam7*RR11*RR1&
  &3*RR33 + RR31*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR12**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1&
  &.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR31 + 2.D0*Lam8*RR12*RR13*RR33 + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS1S1S1S1f3312 = -1.D0*RR31*((3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)*RR31 + (Lam3 + Lam4 + Lam5)*&
  &(RR12*RR21 + RR11*RR22)*RR32 + Lam7*(RR13*RR21 + RR11*RR23)*RR33) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*(RR12*RR21 + RR11*RR22)*RR&
  &31 + ((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23)*RR32 + Lam8*(RR13*RR22 + RR12*RR23)*RR33) - 1.D0*&
  &RR33*(Lam7*(RR13*RR21 + RR11*RR23)*RR31 + Lam8*(RR13*RR22 + RR12*RR23)*RR32 + (Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23)*RR33)
CS1S1S1S1f3313 = -1.D0*RR33*(Lam7*RR31*(RR13*RR31 + RR11*RR33) + Lam8*RR32*(RR13*RR32 + RR12*RR33) + RR33*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR32*(RR12*RR31 + RR11*RR32) + Lam7*RR33*(RR13*RR31 + RR11&
  &*RR33) + RR31*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR33*(RR13*RR32 + RR12*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f3321 = -1.D0*RR31*((3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)*RR31 + (Lam3 + Lam4 + Lam5)*&
  &(RR12*RR21 + RR11*RR22)*RR32 + Lam7*(RR13*RR21 + RR11*RR23)*RR33) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*(RR12*RR21 + RR11*RR22)*RR&
  &31 + ((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23)*RR32 + Lam8*(RR13*RR22 + RR12*RR23)*RR33) - 1.D0*&
  &RR33*(Lam7*(RR13*RR21 + RR11*RR23)*RR31 + Lam8*(RR13*RR22 + RR12*RR23)*RR32 + (Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23)*RR33)
CS1S1S1S1f3322 = -1.D0*RR33*(2.D0*Lam7*RR21*RR23*RR31 + 2.D0*Lam8*RR22*RR23*RR32 + RR33*(Lam7*DBLE(RR21**INT(2.D0)) + Lam8*DBLE(R&
  &R22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR22*RR32 + 2.D0*Lam7*RR21*RR2&
  &3*RR33 + RR31*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR22**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1&
  &.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR22*RR31 + 2.D0*Lam8*RR22*RR23*RR33 + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS1S1S1S1f3323 = -1.D0*RR33*(Lam7*RR31*(RR23*RR31 + RR21*RR33) + Lam8*RR32*(RR23*RR32 + RR22*RR33) + RR33*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR32*(RR22*RR31 + RR21*RR32) + Lam7*RR33*(RR23*RR31 + RR21&
  &*RR33) + RR31*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR33*(RR23*RR32 + RR22*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f3331 = -1.D0*RR33*(Lam7*RR31*(RR13*RR31 + RR11*RR33) + Lam8*RR32*(RR13*RR32 + RR12*RR33) + RR33*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR32*(RR12*RR31 + RR11*RR32) + Lam7*RR33*(RR13*RR31 + RR11&
  &*RR33) + RR31*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR33*(RR13*RR32 + RR12*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f3332 = -1.D0*RR33*(Lam7*RR31*(RR23*RR31 + RR21*RR33) + Lam8*RR32*(RR23*RR32 + RR22*RR33) + RR33*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR32*(RR22*RR31 + RR21*RR32) + Lam7*RR33*(RR23*RR31 + RR21&
  &*RR33) + RR31*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR33*(RR23*RR32 + RR22*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f3333 = -1.D0*RR33*(2.D0*Lam7*RR33*DBLE(RR31**INT(2.D0)) + 2.D0*Lam8*RR33*DBLE(RR32**INT(2.D0)) + RR33*(Lam7*DBLE(RR31**&
  &INT(2.D0)) + Lam8*DBLE(RR32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR31*DBLE(R&
  &R32**INT(2.D0)) + 2.D0*Lam7*RR31*DBLE(RR33**INT(2.D0)) + RR31*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR3&
  &2**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR32*DBLE(RR31**INT(2.D0)) + 2.D0*Lam8*RR3&
  &2*DBLE(RR33**INT(2.D0)) + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33*&
  &*INT(2.D0))))

CS2S2S2S2f1111 = -1.D0*SB*(2.D0*CB2*(Lam3 + Lam4 + Lam5)*SB + SB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam2*SB2)) - 1.D0*CB*(2.D0*CB*(&
  &Lam3 + Lam4 + Lam5)*SB2 + CB*(3.D0*CB2*Lam1 + (Lam3 + Lam4 + Lam5)*SB2))
CS2S2S2S2f1112 = -1.D0*SB*(SB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1&
  &.D0*CB*(CB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + (Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f1121 = -1.D0*SB*(SB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1&
  &.D0*CB*(CB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + (Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f1122 = -1.D0*CB*(-2.D0*CB*(Lam3 + Lam4 + Lam5)*SB2 + CB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam1*SB2)) - 1.D0*SB*(-2.D0*CB&
  &2*(Lam3 + Lam4 + Lam5)*SB + SB*(3.D0*CB2*Lam2 + (Lam3 + Lam4 + Lam5)*SB2))
CS2S2S2S2f1211 = -1.D0*SB*(-2.D0*CB*(Lam3 + Lam4 + Lam5)*SB2 + CB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam2*SB2)) - 1.D0*CB*(2.D0*CB2&
  &*(Lam3 + Lam4 + Lam5)*SB - 1.D0*SB*(3.D0*CB2*Lam1 + (Lam3 + Lam4 + Lam5)*SB2))
CS2S2S2S2f1212 = -1.D0*CB*(-1.D0*SB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) -&
  & 1.D0*SB*(CB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) - 1.D0*(Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f1221 = -1.D0*CB*(-1.D0*SB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) -&
  & 1.D0*SB*(CB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) - 1.D0*(Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f1222 = -1.D0*CB*(-2.D0*CB2*(Lam3 + Lam4 + Lam5)*SB - 1.D0*SB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam1*SB2)) - 1.D0*SB*(2.D&
  &0*CB*(Lam3 + Lam4 + Lam5)*SB2 + CB*(3.D0*CB2*Lam2 + (Lam3 + Lam4 + Lam5)*SB2))
CS2S2S2S2f2111 = -1.D0*CB*(2.D0*CB2*(Lam3 + Lam4 + Lam5)*SB + SB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam2*SB2)) + SB*(2.D0*CB*(Lam3 &
  &+ Lam4 + Lam5)*SB2 + CB*(3.D0*CB2*Lam1 + (Lam3 + Lam4 + Lam5)*SB2))
CS2S2S2S2f2112 = -1.D0*CB*(SB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) + S&
  &B*(CB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + (Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f2121 = -1.D0*CB*(SB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) + S&
  &B*(CB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + (Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f2122 = SB*(-2.D0*CB*(Lam3 + Lam4 + Lam5)*SB2 + CB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam1*SB2)) - 1.D0*CB*(-2.D0*CB2*(Lam&
  &3 + Lam4 + Lam5)*SB + SB*(3.D0*CB2*Lam2 + (Lam3 + Lam4 + Lam5)*SB2))
CS2S2S2S2f2211 = -1.D0*CB*(-2.D0*CB*(Lam3 + Lam4 + Lam5)*SB2 + CB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam2*SB2)) + SB*(2.D0*CB2*(Lam&
  &3 + Lam4 + Lam5)*SB - 1.D0*SB*(3.D0*CB2*Lam1 + (Lam3 + Lam4 + Lam5)*SB2))
CS2S2S2S2f2212 = SB*(-1.D0*SB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*&
  &CB*(CB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) - 1.D0*(Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f2221 = SB*(-1.D0*SB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*&
  &CB*(CB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) - 1.D0*(Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f2222 = SB*(-2.D0*CB2*(Lam3 + Lam4 + Lam5)*SB - 1.D0*SB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam1*SB2)) - 1.D0*CB*(2.D0*CB*(&
  &Lam3 + Lam4 + Lam5)*SB2 + CB*(3.D0*CB2*Lam2 + (Lam3 + Lam4 + Lam5)*SB2))

CS3S3S3S3f1111 = -1.D0*SB*(2.D0*CB2*(Lam3 + Lam4)*SB + 2.D0*SB*(CB2*Lam5 + Lam2*SB2)) - 1.D0*CB*(2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D0&
  &*CB*(CB2*Lam1 + Lam5*SB2))
CS3S3S3S3f1112 = -1.D0*SB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*(-1.D0*&
  &CB*Lam1*SB + CB*Lam5*SB) + (Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f1121 = -1.D0*SB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*(-1.D0*&
  &CB*Lam1*SB + CB*Lam5*SB) + (Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f1122 = -1.D0*CB*(-2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D0*CB*(CB2*Lam5 + Lam1*SB2)) - 1.D0*SB*(-2.D0*CB2*(Lam3 + Lam4)*SB + 2.&
  &D0*SB*(CB2*Lam2 + Lam5*SB2))
CS3S3S3S3f1211 = -1.D0*SB*(-2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D0*CB*(CB2*Lam5 + Lam2*SB2)) - 1.D0*CB*(2.D0*CB2*(Lam3 + Lam4)*SB - 2.D&
  &0*SB*(CB2*Lam1 + Lam5*SB2))
CS3S3S3S3f1212 = -1.D0*CB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) - 1.D0*SB*(2.D0*CB*(CB*L&
  &am2*SB - 1.D0*CB*Lam5*SB) - 1.D0*(Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f1221 = -1.D0*CB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) - 1.D0*SB*(2.D0*CB*(CB*L&
  &am2*SB - 1.D0*CB*Lam5*SB) - 1.D0*(Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f1222 = -1.D0*CB*(-2.D0*CB2*(Lam3 + Lam4)*SB - 2.D0*SB*(CB2*Lam5 + Lam1*SB2)) - 1.D0*SB*(2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D&
  &0*CB*(CB2*Lam2 + Lam5*SB2))
CS3S3S3S3f2111 = -1.D0*CB*(2.D0*CB2*(Lam3 + Lam4)*SB + 2.D0*SB*(CB2*Lam5 + Lam2*SB2)) + SB*(2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D0*CB*(&
  &CB2*Lam1 + Lam5*SB2))
CS3S3S3S3f2112 = -1.D0*CB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) + SB*(2.D0*CB*(-1.D0*CB*La&
  &m1*SB + CB*Lam5*SB) + (Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f2121 = -1.D0*CB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) + SB*(2.D0*CB*(-1.D0*CB*La&
  &m1*SB + CB*Lam5*SB) + (Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f2122 = SB*(-2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D0*CB*(CB2*Lam5 + Lam1*SB2)) - 1.D0*CB*(-2.D0*CB2*(Lam3 + Lam4)*SB + 2.D0*SB*&
  &(CB2*Lam2 + Lam5*SB2))
CS3S3S3S3f2211 = -1.D0*CB*(-2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D0*CB*(CB2*Lam5 + Lam2*SB2)) + SB*(2.D0*CB2*(Lam3 + Lam4)*SB - 2.D0*SB*&
  &(CB2*Lam1 + Lam5*SB2))
CS3S3S3S3f2212 = SB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*(CB*Lam2*SB&
  & - 1.D0*CB*Lam5*SB) - 1.D0*(Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f2221 = SB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*(CB*Lam2*SB&
  & - 1.D0*CB*Lam5*SB) - 1.D0*(Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f2222 = SB*(-2.D0*CB2*(Lam3 + Lam4)*SB - 2.D0*SB*(CB2*Lam5 + Lam1*SB2)) - 1.D0*CB*(2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D0*CB*(&
  &CB2*Lam2 + Lam5*SB2))

CS2S2S1S1f1111 = -1.D0*CB*(2.D0*Lam5*RR11*RR12*SB + CB*(Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR12**INT(2.D&
  &0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1.D0*SB*(2.D0*CB*Lam5*RR11*RR12 + SB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR11**INT(2.D0)) + L&
  &am2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS2S2S1S1f1112 = -1.D0*CB*(CB*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR13*RR23) + Lam5*(RR12*RR21 + RR11*RR&
  &22)*SB) - 1.D0*SB*(CB*Lam5*(RR12*RR21 + RR11*RR22) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23)*S&
  &B)
CS2S2S1S1f1113 = -1.D0*CB*(CB*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR13*RR33) + Lam5*(RR12*RR31 + RR11*RR&
  &32)*SB) - 1.D0*SB*(CB*Lam5*(RR12*RR31 + RR11*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33)*S&
  &B)
CS2S2S1S1f1121 = -1.D0*CB*(CB*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR13*RR23) + Lam5*(RR12*RR21 + RR11*RR&
  &22)*SB) - 1.D0*SB*(CB*Lam5*(RR12*RR21 + RR11*RR22) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23)*S&
  &B)
CS2S2S1S1f1122 = -1.D0*CB*(2.D0*Lam5*RR21*RR22*SB + CB*(Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR22**INT(2.D&
  &0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1.D0*SB*(2.D0*CB*Lam5*RR21*RR22 + SB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR21**INT(2.D0)) + L&
  &am2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS2S2S1S1f1123 = -1.D0*CB*(CB*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR23*RR33) + Lam5*(RR22*RR31 + RR21*RR&
  &32)*SB) - 1.D0*SB*(CB*Lam5*(RR22*RR31 + RR21*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33)*S&
  &B)
CS2S2S1S1f1131 = -1.D0*CB*(CB*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR13*RR33) + Lam5*(RR12*RR31 + RR11*RR&
  &32)*SB) - 1.D0*SB*(CB*Lam5*(RR12*RR31 + RR11*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33)*S&
  &B)
CS2S2S1S1f1132 = -1.D0*CB*(CB*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR23*RR33) + Lam5*(RR22*RR31 + RR21*RR&
  &32)*SB) - 1.D0*SB*(CB*Lam5*(RR22*RR31 + RR21*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33)*S&
  &B)
CS2S2S1S1f1133 = -1.D0*CB*(2.D0*Lam5*RR31*RR32*SB + CB*(Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR32**INT(2.D&
  &0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1.D0*SB*(2.D0*CB*Lam5*RR31*RR32 + SB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR31**INT(2.D0)) + L&
  &am2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS2S2S1S1f1211 = -1.D0*CB*(2.D0*CB*Lam5*RR11*RR12 - 1.D0*SB*(Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR12**IN&
  &T(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1.D0*SB*(-2.D0*Lam5*RR11*RR12*SB + CB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR11**INT(2.D0&
  &)) + Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS2S2S1S1f1212 = -1.D0*SB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23) - 1.D0*Lam5*(RR12*RR21 + RR&
  &11*RR22)*SB) - 1.D0*CB*(CB*Lam5*(RR12*RR21 + RR11*RR22) - 1.D0*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR&
  &13*RR23)*SB)
CS2S2S1S1f1213 = -1.D0*SB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33) - 1.D0*Lam5*(RR12*RR31 + RR&
  &11*RR32)*SB) - 1.D0*CB*(CB*Lam5*(RR12*RR31 + RR11*RR32) - 1.D0*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR&
  &13*RR33)*SB)
CS2S2S1S1f1221 = -1.D0*SB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23) - 1.D0*Lam5*(RR12*RR21 + RR&
  &11*RR22)*SB) - 1.D0*CB*(CB*Lam5*(RR12*RR21 + RR11*RR22) - 1.D0*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR&
  &13*RR23)*SB)
CS2S2S1S1f1222 = -1.D0*CB*(2.D0*CB*Lam5*RR21*RR22 - 1.D0*SB*(Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR22**IN&
  &T(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1.D0*SB*(-2.D0*Lam5*RR21*RR22*SB + CB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR21**INT(2.D0&
  &)) + Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS2S2S1S1f1223 = -1.D0*SB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33) - 1.D0*Lam5*(RR22*RR31 + RR&
  &21*RR32)*SB) - 1.D0*CB*(CB*Lam5*(RR22*RR31 + RR21*RR32) - 1.D0*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR&
  &23*RR33)*SB)
CS2S2S1S1f1231 = -1.D0*SB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33) - 1.D0*Lam5*(RR12*RR31 + RR&
  &11*RR32)*SB) - 1.D0*CB*(CB*Lam5*(RR12*RR31 + RR11*RR32) - 1.D0*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR&
  &13*RR33)*SB)
CS2S2S1S1f1232 = -1.D0*SB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33) - 1.D0*Lam5*(RR22*RR31 + RR&
  &21*RR32)*SB) - 1.D0*CB*(CB*Lam5*(RR22*RR31 + RR21*RR32) - 1.D0*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR&
  &23*RR33)*SB)
CS2S2S1S1f1233 = -1.D0*CB*(2.D0*CB*Lam5*RR31*RR32 - 1.D0*SB*(Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR32**IN&
  &T(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1.D0*SB*(-2.D0*Lam5*RR31*RR32*SB + CB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR31**INT(2.D0&
  &)) + Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS2S2S1S1f2111 = SB*(2.D0*Lam5*RR11*RR12*SB + CB*(Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR12**INT(2.D0)) + &
  &Lam7*DBLE(RR13**INT(2.D0)))) - 1.D0*CB*(2.D0*CB*Lam5*RR11*RR12 + SB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR11**INT(2.D0)) + Lam2*DB&
  &LE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS2S2S1S1f2112 = SB*(CB*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR13*RR23) + Lam5*(RR12*RR21 + RR11*RR22)*SB&
  &) - 1.D0*CB*(CB*Lam5*(RR12*RR21 + RR11*RR22) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23)*SB)
CS2S2S1S1f2113 = SB*(CB*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR13*RR33) + Lam5*(RR12*RR31 + RR11*RR32)*SB&
  &) - 1.D0*CB*(CB*Lam5*(RR12*RR31 + RR11*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33)*SB)
CS2S2S1S1f2121 = SB*(CB*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR13*RR23) + Lam5*(RR12*RR21 + RR11*RR22)*SB&
  &) - 1.D0*CB*(CB*Lam5*(RR12*RR21 + RR11*RR22) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23)*SB)
CS2S2S1S1f2122 = SB*(2.D0*Lam5*RR21*RR22*SB + CB*(Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR22**INT(2.D0)) + &
  &Lam7*DBLE(RR23**INT(2.D0)))) - 1.D0*CB*(2.D0*CB*Lam5*RR21*RR22 + SB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR21**INT(2.D0)) + Lam2*DB&
  &LE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS2S2S1S1f2123 = SB*(CB*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR23*RR33) + Lam5*(RR22*RR31 + RR21*RR32)*SB&
  &) - 1.D0*CB*(CB*Lam5*(RR22*RR31 + RR21*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33)*SB)
CS2S2S1S1f2131 = SB*(CB*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR13*RR33) + Lam5*(RR12*RR31 + RR11*RR32)*SB&
  &) - 1.D0*CB*(CB*Lam5*(RR12*RR31 + RR11*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33)*SB)
CS2S2S1S1f2132 = SB*(CB*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR23*RR33) + Lam5*(RR22*RR31 + RR21*RR32)*SB&
  &) - 1.D0*CB*(CB*Lam5*(RR22*RR31 + RR21*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33)*SB)
CS2S2S1S1f2133 = SB*(2.D0*Lam5*RR31*RR32*SB + CB*(Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR32**INT(2.D0)) + &
  &Lam7*DBLE(RR33**INT(2.D0)))) - 1.D0*CB*(2.D0*CB*Lam5*RR31*RR32 + SB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR31**INT(2.D0)) + Lam2*DB&
  &LE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS2S2S1S1f2211 = SB*(2.D0*CB*Lam5*RR11*RR12 - 1.D0*SB*(Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR12**INT(2.D0&
  &)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1.D0*CB*(-2.D0*Lam5*RR11*RR12*SB + CB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR11**INT(2.D0)) + L&
  &am2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS2S2S1S1f2212 = -1.D0*CB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23) - 1.D0*Lam5*(RR12*RR21 + RR&
  &11*RR22)*SB) + SB*(CB*Lam5*(RR12*RR21 + RR11*RR22) - 1.D0*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR13*RR&
  &23)*SB)
CS2S2S1S1f2213 = -1.D0*CB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33) - 1.D0*Lam5*(RR12*RR31 + RR&
  &11*RR32)*SB) + SB*(CB*Lam5*(RR12*RR31 + RR11*RR32) - 1.D0*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR13*RR&
  &33)*SB)
CS2S2S1S1f2221 = -1.D0*CB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23) - 1.D0*Lam5*(RR12*RR21 + RR&
  &11*RR22)*SB) + SB*(CB*Lam5*(RR12*RR21 + RR11*RR22) - 1.D0*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR13*RR&
  &23)*SB)
CS2S2S1S1f2222 = SB*(2.D0*CB*Lam5*RR21*RR22 - 1.D0*SB*(Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR22**INT(2.D0&
  &)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1.D0*CB*(-2.D0*Lam5*RR21*RR22*SB + CB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR21**INT(2.D0)) + L&
  &am2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS2S2S1S1f2223 = -1.D0*CB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33) - 1.D0*Lam5*(RR22*RR31 + RR&
  &21*RR32)*SB) + SB*(CB*Lam5*(RR22*RR31 + RR21*RR32) - 1.D0*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR23*RR&
  &33)*SB)
CS2S2S1S1f2231 = -1.D0*CB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33) - 1.D0*Lam5*(RR12*RR31 + RR&
  &11*RR32)*SB) + SB*(CB*Lam5*(RR12*RR31 + RR11*RR32) - 1.D0*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR13*RR&
  &33)*SB)
CS2S2S1S1f2232 = -1.D0*CB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33) - 1.D0*Lam5*(RR22*RR31 + RR&
  &21*RR32)*SB) + SB*(CB*Lam5*(RR22*RR31 + RR21*RR32) - 1.D0*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR23*RR&
  &33)*SB)
CS2S2S1S1f2233 = SB*(2.D0*CB*Lam5*RR31*RR32 - 1.D0*SB*(Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR32**INT(2.D0&
  &)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1.D0*CB*(-2.D0*Lam5*RR31*RR32*SB + CB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR31**INT(2.D0)) + L&
  &am2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))

CS2S2S3S3f1111 = 0.5D0*(-1.D0*SB*(2.D0*CB2*(Lam4 + Lam5)*SB + 2.D0*SB*(CB2*Lam3 + Lam2*SB2)) - 1.D0*CB*(2.D0*CB*(Lam4 + Lam5)*SB2&
  & + 2.D0*CB*(CB2*Lam1 + Lam3*SB2)))
CS2S2S3S3f1112 = 0.5D0*(-1.D0*SB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*&
  &(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f1121 = 0.5D0*(-1.D0*SB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*&
  &(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f1122 = 0.5D0*(-1.D0*CB*(-2.D0*CB*(Lam4 + Lam5)*SB2 + 2.D0*CB*(CB2*Lam3 + Lam1*SB2)) - 1.D0*SB*(-2.D0*CB2*(Lam4 + Lam5)*&
  &SB + 2.D0*SB*(CB2*Lam2 + Lam3*SB2)))
CS2S2S3S3f1211 = 0.5D0*(-1.D0*SB*(-2.D0*CB*(Lam4 + Lam5)*SB2 + 2.D0*CB*(CB2*Lam3 + Lam2*SB2)) - 1.D0*CB*(2.D0*CB2*(Lam4 + Lam5)*S&
  &B - 2.D0*SB*(CB2*Lam1 + Lam3*SB2)))
CS2S2S3S3f1212 = 0.5D0*(-1.D0*CB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*SB*(2.D0*C&
  &B*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) - 1.D0*(Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f1221 = 0.5D0*(-1.D0*CB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*SB*(2.D0*C&
  &B*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) - 1.D0*(Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f1222 = 0.5D0*(-1.D0*CB*(-2.D0*CB2*(Lam4 + Lam5)*SB - 2.D0*SB*(CB2*Lam3 + Lam1*SB2)) - 1.D0*SB*(2.D0*CB*(Lam4 + Lam5)*SB&
  &2 + 2.D0*CB*(CB2*Lam2 + Lam3*SB2)))
CS2S2S3S3f2111 = 0.5D0*(-1.D0*CB*(2.D0*CB2*(Lam4 + Lam5)*SB + 2.D0*SB*(CB2*Lam3 + Lam2*SB2)) + SB*(2.D0*CB*(Lam4 + Lam5)*SB2 + 2.&
  &D0*CB*(CB2*Lam1 + Lam3*SB2)))
CS2S2S3S3f2112 = 0.5D0*(-1.D0*CB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) + SB*(2.D0*CB*(-1.D&
  &0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f2121 = 0.5D0*(-1.D0*CB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) + SB*(2.D0*CB*(-1.D&
  &0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f2122 = 0.5D0*(SB*(-2.D0*CB*(Lam4 + Lam5)*SB2 + 2.D0*CB*(CB2*Lam3 + Lam1*SB2)) - 1.D0*CB*(-2.D0*CB2*(Lam4 + Lam5)*SB + 2&
  &.D0*SB*(CB2*Lam2 + Lam3*SB2)))
CS2S2S3S3f2211 = 0.5D0*(-1.D0*CB*(-2.D0*CB*(Lam4 + Lam5)*SB2 + 2.D0*CB*(CB2*Lam3 + Lam2*SB2)) + SB*(2.D0*CB2*(Lam4 + Lam5)*SB - 2&
  &.D0*SB*(CB2*Lam1 + Lam3*SB2)))
CS2S2S3S3f2212 = 0.5D0*(SB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*(CB*&
  &Lam2*SB - 1.D0*CB*Lam3*SB) - 1.D0*(Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f2221 = 0.5D0*(SB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*(CB*&
  &Lam2*SB - 1.D0*CB*Lam3*SB) - 1.D0*(Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f2222 = 0.5D0*(SB*(-2.D0*CB2*(Lam4 + Lam5)*SB - 2.D0*SB*(CB2*Lam3 + Lam1*SB2)) - 1.D0*CB*(2.D0*CB*(Lam4 + Lam5)*SB2 + 2.&
  &D0*CB*(CB2*Lam2 + Lam3*SB2)))

CS1S1S3S3f1111 = 0.5D0*(-1.D0*RR12*(2.D0*CB*(Lam4 + Lam5)*RR11*SB + 2.D0*RR12*(CB2*Lam3 + Lam2*SB2)) - 1.D0*RR11*(2.D0*CB*(Lam4 +&
  & Lam5)*RR12*SB + 2.D0*RR11*(CB2*Lam1 + Lam3*SB2)) - 2.D0*(CB2*Lam7 + Lam8*SB2)*DBLE(RR13**INT(2.D0)))
CS1S1S3S3f1112 = 0.5D0*(-1.D0*RR12*(2.D0*RR12*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (Lam4 + Lam5)*RR11*(CB2 - 1.D0*SB2)) - 1.D0*RR11*(&
  &2.D0*RR11*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR12*(CB2 - 1.D0*SB2)) - 2.D0*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*DBLE(R&
  &R13**INT(2.D0)))
CS1S1S3S3f1121 = 0.5D0*(-1.D0*RR12*(2.D0*RR12*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (Lam4 + Lam5)*RR11*(CB2 - 1.D0*SB2)) - 1.D0*RR11*(&
  &2.D0*RR11*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR12*(CB2 - 1.D0*SB2)) - 2.D0*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*DBLE(R&
  &R13**INT(2.D0)))
CS1S1S3S3f1122 = 0.5D0*(-1.D0*RR11*(-2.D0*CB*(Lam4 + Lam5)*RR12*SB + 2.D0*RR11*(CB2*Lam3 + Lam1*SB2)) - 1.D0*RR12*(-2.D0*CB*(Lam4&
  & + Lam5)*RR11*SB + 2.D0*RR12*(CB2*Lam2 + Lam3*SB2)) - 2.D0*(CB2*Lam8 + Lam7*SB2)*DBLE(RR13**INT(2.D0)))
CS1S1S3S3f1211 = 0.5D0*(-2.D0*RR13*RR23*(CB2*Lam7 + Lam8*SB2) - 1.D0*RR12*(2.D0*CB*(Lam4 + Lam5)*RR21*SB + 2.D0*RR22*(CB2*Lam3 + &
  &Lam2*SB2)) - 1.D0*RR11*(2.D0*CB*(Lam4 + Lam5)*RR22*SB + 2.D0*RR21*(CB2*Lam1 + Lam3*SB2)))
CS1S1S3S3f1212 = 0.5D0*(-2.D0*RR13*RR23*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR12*(2.D0*RR22*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR21*(CB2 - 1.D0*SB2)) - 1.D0*RR11*(2.D0*RR21*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR22*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f1221 = 0.5D0*(-2.D0*RR13*RR23*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR12*(2.D0*RR22*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR21*(CB2 - 1.D0*SB2)) - 1.D0*RR11*(2.D0*RR21*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR22*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f1222 = 0.5D0*(-2.D0*RR13*RR23*(CB2*Lam8 + Lam7*SB2) - 1.D0*RR11*(-2.D0*CB*(Lam4 + Lam5)*RR22*SB + 2.D0*RR21*(CB2*Lam3 +&
  & Lam1*SB2)) - 1.D0*RR12*(-2.D0*CB*(Lam4 + Lam5)*RR21*SB + 2.D0*RR22*(CB2*Lam2 + Lam3*SB2)))
CS1S1S3S3f1311 = 0.5D0*(-2.D0*RR13*RR33*(CB2*Lam7 + Lam8*SB2) - 1.D0*RR12*(2.D0*CB*(Lam4 + Lam5)*RR31*SB + 2.D0*RR32*(CB2*Lam3 + &
  &Lam2*SB2)) - 1.D0*RR11*(2.D0*CB*(Lam4 + Lam5)*RR32*SB + 2.D0*RR31*(CB2*Lam1 + Lam3*SB2)))
CS1S1S3S3f1312 = 0.5D0*(-2.D0*RR13*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR12*(2.D0*RR32*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR31*(CB2 - 1.D0*SB2)) - 1.D0*RR11*(2.D0*RR31*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR32*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f1321 = 0.5D0*(-2.D0*RR13*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR12*(2.D0*RR32*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR31*(CB2 - 1.D0*SB2)) - 1.D0*RR11*(2.D0*RR31*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR32*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f1322 = 0.5D0*(-2.D0*RR13*RR33*(CB2*Lam8 + Lam7*SB2) - 1.D0*RR11*(-2.D0*CB*(Lam4 + Lam5)*RR32*SB + 2.D0*RR31*(CB2*Lam3 +&
  & Lam1*SB2)) - 1.D0*RR12*(-2.D0*CB*(Lam4 + Lam5)*RR31*SB + 2.D0*RR32*(CB2*Lam2 + Lam3*SB2)))
CS1S1S3S3f2111 = 0.5D0*(-2.D0*RR13*RR23*(CB2*Lam7 + Lam8*SB2) - 1.D0*RR22*(2.D0*CB*(Lam4 + Lam5)*RR11*SB + 2.D0*RR12*(CB2*Lam3 + &
  &Lam2*SB2)) - 1.D0*RR21*(2.D0*CB*(Lam4 + Lam5)*RR12*SB + 2.D0*RR11*(CB2*Lam1 + Lam3*SB2)))
CS1S1S3S3f2112 = 0.5D0*(-2.D0*RR13*RR23*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR22*(2.D0*RR12*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR11*(CB2 - 1.D0*SB2)) - 1.D0*RR21*(2.D0*RR11*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR12*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f2121 = 0.5D0*(-2.D0*RR13*RR23*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR22*(2.D0*RR12*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR11*(CB2 - 1.D0*SB2)) - 1.D0*RR21*(2.D0*RR11*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR12*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f2122 = 0.5D0*(-2.D0*RR13*RR23*(CB2*Lam8 + Lam7*SB2) - 1.D0*RR21*(-2.D0*CB*(Lam4 + Lam5)*RR12*SB + 2.D0*RR11*(CB2*Lam3 +&
  & Lam1*SB2)) - 1.D0*RR22*(-2.D0*CB*(Lam4 + Lam5)*RR11*SB + 2.D0*RR12*(CB2*Lam2 + Lam3*SB2)))
CS1S1S3S3f2211 = 0.5D0*(-1.D0*RR22*(2.D0*CB*(Lam4 + Lam5)*RR21*SB + 2.D0*RR22*(CB2*Lam3 + Lam2*SB2)) - 1.D0*RR21*(2.D0*CB*(Lam4 +&
  & Lam5)*RR22*SB + 2.D0*RR21*(CB2*Lam1 + Lam3*SB2)) - 2.D0*(CB2*Lam7 + Lam8*SB2)*DBLE(RR23**INT(2.D0)))
CS1S1S3S3f2212 = 0.5D0*(-1.D0*RR22*(2.D0*RR22*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (Lam4 + Lam5)*RR21*(CB2 - 1.D0*SB2)) - 1.D0*RR21*(&
  &2.D0*RR21*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR22*(CB2 - 1.D0*SB2)) - 2.D0*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*DBLE(R&
  &R23**INT(2.D0)))
CS1S1S3S3f2221 = 0.5D0*(-1.D0*RR22*(2.D0*RR22*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (Lam4 + Lam5)*RR21*(CB2 - 1.D0*SB2)) - 1.D0*RR21*(&
  &2.D0*RR21*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR22*(CB2 - 1.D0*SB2)) - 2.D0*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*DBLE(R&
  &R23**INT(2.D0)))
CS1S1S3S3f2222 = 0.5D0*(-1.D0*RR21*(-2.D0*CB*(Lam4 + Lam5)*RR22*SB + 2.D0*RR21*(CB2*Lam3 + Lam1*SB2)) - 1.D0*RR22*(-2.D0*CB*(Lam4&
  & + Lam5)*RR21*SB + 2.D0*RR22*(CB2*Lam2 + Lam3*SB2)) - 2.D0*(CB2*Lam8 + Lam7*SB2)*DBLE(RR23**INT(2.D0)))
CS1S1S3S3f2311 = 0.5D0*(-2.D0*RR23*RR33*(CB2*Lam7 + Lam8*SB2) - 1.D0*RR22*(2.D0*CB*(Lam4 + Lam5)*RR31*SB + 2.D0*RR32*(CB2*Lam3 + &
  &Lam2*SB2)) - 1.D0*RR21*(2.D0*CB*(Lam4 + Lam5)*RR32*SB + 2.D0*RR31*(CB2*Lam1 + Lam3*SB2)))
CS1S1S3S3f2312 = 0.5D0*(-2.D0*RR23*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR22*(2.D0*RR32*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR31*(CB2 - 1.D0*SB2)) - 1.D0*RR21*(2.D0*RR31*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR32*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f2321 = 0.5D0*(-2.D0*RR23*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR22*(2.D0*RR32*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR31*(CB2 - 1.D0*SB2)) - 1.D0*RR21*(2.D0*RR31*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR32*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f2322 = 0.5D0*(-2.D0*RR23*RR33*(CB2*Lam8 + Lam7*SB2) - 1.D0*RR21*(-2.D0*CB*(Lam4 + Lam5)*RR32*SB + 2.D0*RR31*(CB2*Lam3 +&
  & Lam1*SB2)) - 1.D0*RR22*(-2.D0*CB*(Lam4 + Lam5)*RR31*SB + 2.D0*RR32*(CB2*Lam2 + Lam3*SB2)))
CS1S1S3S3f3111 = 0.5D0*(-2.D0*RR13*RR33*(CB2*Lam7 + Lam8*SB2) - 1.D0*RR32*(2.D0*CB*(Lam4 + Lam5)*RR11*SB + 2.D0*RR12*(CB2*Lam3 + &
  &Lam2*SB2)) - 1.D0*RR31*(2.D0*CB*(Lam4 + Lam5)*RR12*SB + 2.D0*RR11*(CB2*Lam1 + Lam3*SB2)))
CS1S1S3S3f3112 = 0.5D0*(-2.D0*RR13*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR32*(2.D0*RR12*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR11*(CB2 - 1.D0*SB2)) - 1.D0*RR31*(2.D0*RR11*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR12*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f3121 = 0.5D0*(-2.D0*RR13*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR32*(2.D0*RR12*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR11*(CB2 - 1.D0*SB2)) - 1.D0*RR31*(2.D0*RR11*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR12*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f3122 = 0.5D0*(-2.D0*RR13*RR33*(CB2*Lam8 + Lam7*SB2) - 1.D0*RR31*(-2.D0*CB*(Lam4 + Lam5)*RR12*SB + 2.D0*RR11*(CB2*Lam3 +&
  & Lam1*SB2)) - 1.D0*RR32*(-2.D0*CB*(Lam4 + Lam5)*RR11*SB + 2.D0*RR12*(CB2*Lam2 + Lam3*SB2)))
CS1S1S3S3f3211 = 0.5D0*(-2.D0*RR23*RR33*(CB2*Lam7 + Lam8*SB2) - 1.D0*RR32*(2.D0*CB*(Lam4 + Lam5)*RR21*SB + 2.D0*RR22*(CB2*Lam3 + &
  &Lam2*SB2)) - 1.D0*RR31*(2.D0*CB*(Lam4 + Lam5)*RR22*SB + 2.D0*RR21*(CB2*Lam1 + Lam3*SB2)))
CS1S1S3S3f3212 = 0.5D0*(-2.D0*RR23*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR32*(2.D0*RR22*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR21*(CB2 - 1.D0*SB2)) - 1.D0*RR31*(2.D0*RR21*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR22*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f3221 = 0.5D0*(-2.D0*RR23*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR32*(2.D0*RR22*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR21*(CB2 - 1.D0*SB2)) - 1.D0*RR31*(2.D0*RR21*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR22*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f3222 = 0.5D0*(-2.D0*RR23*RR33*(CB2*Lam8 + Lam7*SB2) - 1.D0*RR31*(-2.D0*CB*(Lam4 + Lam5)*RR22*SB + 2.D0*RR21*(CB2*Lam3 +&
  & Lam1*SB2)) - 1.D0*RR32*(-2.D0*CB*(Lam4 + Lam5)*RR21*SB + 2.D0*RR22*(CB2*Lam2 + Lam3*SB2)))
CS1S1S3S3f3311 = 0.5D0*(-1.D0*RR32*(2.D0*CB*(Lam4 + Lam5)*RR31*SB + 2.D0*RR32*(CB2*Lam3 + Lam2*SB2)) - 1.D0*RR31*(2.D0*CB*(Lam4 +&
  & Lam5)*RR32*SB + 2.D0*RR31*(CB2*Lam1 + Lam3*SB2)) - 2.D0*(CB2*Lam7 + Lam8*SB2)*DBLE(RR33**INT(2.D0)))
CS1S1S3S3f3312 = 0.5D0*(-1.D0*RR32*(2.D0*RR32*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (Lam4 + Lam5)*RR31*(CB2 - 1.D0*SB2)) - 1.D0*RR31*(&
  &2.D0*RR31*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR32*(CB2 - 1.D0*SB2)) - 2.D0*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*DBLE(R&
  &R33**INT(2.D0)))
CS1S1S3S3f3321 = 0.5D0*(-1.D0*RR32*(2.D0*RR32*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (Lam4 + Lam5)*RR31*(CB2 - 1.D0*SB2)) - 1.D0*RR31*(&
  &2.D0*RR31*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR32*(CB2 - 1.D0*SB2)) - 2.D0*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*DBLE(R&
  &R33**INT(2.D0)))
CS1S1S3S3f3322 = 0.5D0*(-1.D0*RR31*(-2.D0*CB*(Lam4 + Lam5)*RR32*SB + 2.D0*RR31*(CB2*Lam3 + Lam1*SB2)) - 1.D0*RR32*(-2.D0*CB*(Lam4&
  & + Lam5)*RR31*SB + 2.D0*RR32*(CB2*Lam2 + Lam3*SB2)) - 2.D0*(CB2*Lam8 + Lam7*SB2)*DBLE(RR33**INT(2.D0)))

        ! Calculate the full amplitude for different lambda
        call ltini
            call setdelta(0D0)
            call setmudim(1D0)
            call setlambda(1D0)
            call clearcache
            IRLambda = getlambda()

            ! Reset the UV scale 
            UVDelta = 0D0

            ! Kinematic prefactor together with the symmetry factor of the process
            prefactor = 1D0/1D0 * DSQRT(MA0**4 + MD**4 + MD**4 - 2D0*MA0**2*MD**2 - 2D0*MA0**2*MD**2 &
                        &  - 2D0*MD**2*MD**2 )/(16D0*PI*MA0**3)

            ! Perform the check for gauge dependence for all schemes
            do n = 1, maxNumberSchemes, 1
                ! Reset the checking arrays for all schemes
                isIRDivergentList(n) = .false.
                isIRDivergent = .false.

                ! Info header
                print *, " "
                print *, "==========="
                print *, "IR divergence of scheme ", n
                print *, "==========="

                ! Calculate the amplitude for several different lambda values
                ! Schemes 1, 2, and 9 are without tadpoles
                if ((n == 1) .OR. (n == 2) .OR. (n == 9)) then
                    IRChecks(1) = prefactor*( 2D0*DBLE(A0toDDBarVC()) + &
                                & 2D0*A0toDDBarCT(n) + &
                                & A0toDDBarReal() )
                else
                    IRChecks(1) = prefactor*( 2D0*DBLE(A0toDDBarVC() + &
                                & A0toDDBarTad()) + 2D0*A0toDDBarCT(n) + &
                                & A0toDDBarReal() )
                end if
                do m = 2, 10, 1
                    call clearcache
                    call setlambda(DBLE(10**(m-1)))
                    IRLambda = getlambda()
                    ! Schemes 1, 2, and 9 are without tadpoles
                    if ((n == 1) .OR. (n == 2) .OR. (n == 9)) then
                        IRChecks(m) = ( prefactor*( 2D0*DBLE(A0toDDBarVC()) + &
                                    & 2D0*A0toDDBarCT(n) + A0toDDBarReal() ) &
                                    & - IRChecks(1))/IRChecks(1)
                    else
                        IRChecks(m) = ( prefactor*( 2D0*DBLE(A0toDDBarVC() + &
                                    & A0toDDBarTad()) + 2D0*A0toDDBarCT(n) + &
                                    & A0toDDBarReal() ) - IRChecks(1))/IRChecks(1)
                    end if
                    print *, "Lambda: ", getlambda(), ", Difference from lambda=1: ", IRChecks(m)
                    print *, "Real: ", A0toDDBarReal()
                    if (abs(IRChecks(m)) > IRDivergenceThreshold) then
                        isIRDivergent = .true.
                        isIRDivergentList(n) = .true.
                    end if
                end do
                call clearcache
                call setlambda(1D0)
                IRLambda = getlambda()
                print *, "==========="

                ! Exception handling: if IR divergences are found, inform the user and ask if the program shall be terminated
                if (isIRDivergent) then
                    write (*, '(/,A)') "WARNING: potential IR divergence found! Please check the vertex corrections, &
                        &the counterterm and the real corrections!"
                    do
                        print *, ">>> Do you want to continue with the evaluation of the program? [y/n]"
                        read (*,*) isIRDivergentContinue
                        if (isIRDivergentContinue == 'n') then
                            print *, "Termination requested by user. 2HDMCalc will be terminated now."
                            stop
                        else if (isIRDivergentContinue == 'y') then
                            exit
                        else
                            print *, "Invalid character. Enter y or n."
                        end if
                    end do
                end if
            end do
        call ltexi

        ! Print the summary of the checks for IR divergence
        print *, "==========="
        print *, "Results of the check on IR divergence"
        print *, "==========="
        do n = 1, maxNumberSchemes, 1
            if (isIRDivergentList(n)) then
                print *, "Scheme ", n, ": potentially IR-divergent."
            else
                print *, "Scheme ", n, ": most likely IR-finite."
            end if
        end do
        print *, "==========="

    end if

    ! Perform the numerical evaluation
    if (arguments(4) == 1) then
        print *, "Starting the numerical evaluation ..."

        ! Reset all values
        GaugeXiA = 1D0
        GaugeXiW = 1D0
        GaugeXiZ = 1D0

        ! Calculate all values
        call ltini
            ! Set default values for the loop calculations
            call setlambda(1D0)
            call setdelta(0D0)
            IRLambda = getlambda()

            ! Use this hack to "fill up" the string to the maximum length with whitespace characters so that it can be passed to the subroutine call
            fileName = fileName // ' '
            targetName = targetName // ' '

            ! Get all parameters
            call getParameters(fileName)
            call setmudim(InputScale**2)

            ! Prepare the output file header
            outputFileContent = "MH1,"
            outputFileContent = trim(outputFileContent) // "MH2,"
            outputFileContent = trim(outputFileContent) // "MH3,"
            outputFileContent = trim(outputFileContent) // "MA0,"
            outputFileContent = trim(outputFileContent) // "MHp,"
            outputFileContent = trim(outputFileContent) // "alpha1,"
            outputFileContent = trim(outputFileContent) // "alpha2,"
            outputFileContent = trim(outputFileContent) // "alpha3,"
            outputFileContent = trim(outputFileContent) // "beta,"
            outputFileContent = trim(outputFileContent) // "m122,"
            outputFileContent = trim(outputFileContent) // "2HDMType,"
            outputFileContent = trim(outputFileContent) // "InputScale,"
            outputFileContent = trim(outputFileContent) // "WidthLO,"
            if (debugModeOn) then
                ! outputFileContent = trim(outputFileContent) // "WidthNLOVC,"
                ! outputFileContent = trim(outputFileContent) // "WidthNLOVCwoIR,"
                ! outputFileContent = trim(outputFileContent) // "WidthIRonly,"
                ! outputFileContent = trim(outputFileContent) // "dMW2Usual,"
                ! outputFileContent = trim(outputFileContent) // "dMW2Alter,"
                ! outputFileContent = trim(outputFileContent) // "dMZ2Usual,"
                ! outputFileContent = trim(outputFileContent) // "dMZ2Alter,"
                ! outputFileContent = trim(outputFileContent) // "dMLOSUsual,"
                ! outputFileContent = trim(outputFileContent) // "dMLOSAlter,"
                ! outputFileContent = trim(outputFileContent) // "dMBOSUsual,"
                ! outputFileContent = trim(outputFileContent) // "dMBOSAlter,"
                ! outputFileContent = trim(outputFileContent) // "dZHHHHOS,"
                ! outputFileContent = trim(outputFileContent) // "dZHHh0OSUsual,"
                ! outputFileContent = trim(outputFileContent) // "dZHHh0OSAlter,"
                ! outputFileContent = trim(outputFileContent) // "dZh0HHOSUsual,"
                ! outputFileContent = trim(outputFileContent) // "dZh0HHOSAlter,"
                ! outputFileContent = trim(outputFileContent) // "dZh0h0OS,"
                ! outputFileContent = trim(outputFileContent) // "dZA0A0OS,"
                ! outputFileContent = trim(outputFileContent) // "dZG0A0OSUsual,"
                ! outputFileContent = trim(outputFileContent) // "dZG0A0OSAlter,"
                ! outputFileContent = trim(outputFileContent) // "dZBotBotOSLeft,"
                ! outputFileContent = trim(outputFileContent) // "dZBotBotOSRight,"
                ! outputFileContent = trim(outputFileContent) // "dZTauTauOSLeft,"
                ! outputFileContent = trim(outputFileContent) // "dZTauTauOSRight,"
                ! outputFileContent = trim(outputFileContent) // "dBeta1KanUsual,"
                ! outputFileContent = trim(outputFileContent) // "dBeta1KanAlter,"
                ! outputFileContent = trim(outputFileContent) // "dBeta1PinchPStar,"
                ! outputFileContent = trim(outputFileContent) // "dBeta1PinchOS,"
            end if
            outputFileContent = trim(outputFileContent) // "WidthKanOdd,"
            outputFileContent = trim(outputFileContent) // "WidthKanChar,"
            outputFileContent = trim(outputFileContent) // "WidthTadOdd,"
            outputFileContent = trim(outputFileContent) // "WidthTadChar,"
            outputFileContent = trim(outputFileContent) // "WidthPStarOdd,"
            outputFileContent = trim(outputFileContent) // "WidthPStarChar,"
            outputFileContent = trim(outputFileContent) // "WidthOSPinOdd,"
            outputFileContent = trim(outputFileContent) // "WidthOSPinChar,"
            outputFileContent = trim(outputFileContent) // "WidthUsuProcDep1,"
            outputFileContent = trim(outputFileContent) // "WidthAltProcDep1,"
            outputFileContent = trim(outputFileContent) // "WidthUsuProcDep2,"
            outputFileContent = trim(outputFileContent) // "WidthAltProcDep2,"
            outputFileContent = trim(outputFileContent) // "WidthUsuProcDep3,"
            outputFileContent = trim(outputFileContent) // "WidthAltProcDep3,"
            outputFileContent = trim(outputFileContent) // "DifWidthLO,"
            outputFileContent = trim(outputFileContent) // "DifWidthKanOdd,"
            outputFileContent = trim(outputFileContent) // "DifWidthKanChar,"
            outputFileContent = trim(outputFileContent) // "DifWidthTadOdd,"
            outputFileContent = trim(outputFileContent) // "DifWidthTadChar,"
            outputFileContent = trim(outputFileContent) // "DifWidthPStarOdd,"
            outputFileContent = trim(outputFileContent) // "DifWidthPStarChar,"
            outputFileContent = trim(outputFileContent) // "DifWidthOSPinOdd,"
            outputFileContent = trim(outputFileContent) // "DifWidthOSPinChar,"
            outputFileContent = trim(outputFileContent) // "DifWidthUsuProcDep1,"
            outputFileContent = trim(outputFileContent) // "DifWidthAltProcDep1,"
            outputFileContent = trim(outputFileContent) // "DifWidthUsuProcDep2,"
            outputFileContent = trim(outputFileContent) // "DifWidthAltProcDep2,"
            outputFileContent = trim(outputFileContent) // "DifWidthUsuProcDep3,"
            outputFileContent = trim(outputFileContent) // "DifWidthAltProcDep3\n"

            do point = 1, maxPoint, 1
                ! Set the 2HDM parameters according to the current point in phase-space
                MH2 = MH2List(point)
                MH3 = MH3List(point)
                MA0 = MA0List(point)
                MHp = MHpList(point)
                alpha1 = alpha1List(point)
                alpha2 = alpha2List(point)
                alpha3 = alpha3List(point)
                beta = betaList(point)
                CA1 = CA1List(point)
                CA2 = CA2List(point)
                CA3 = CA3List(point)
                SA1 = SA1List(point)
                SA2 = SA2List(point)
                SA3 = SA3List(point)
                CB = CBList(point)
                SB = SBList(point)
                TB = TBList(point)
                YukS1Lep1 = YukS1Lep1List(point)
                YukS1Lep2 = YukS1Lep2List(point)
                YukS1Lep3 = YukS1Lep3List(point)
                YukS2Lep1 = YukS2Lep1List(point)
                YukS2Lep2 = YukS2Lep2List(point)
                YukS3Lep1 = YukS3Lep1List(point)
                YukS3Lep2 = YukS3Lep2List(point)
                YukS1Quark1 = YukS1Quark1List(point)
                YukS1Quark2 = YukS1Quark2List(point)
                YukS1Quark3 = YukS1Quark3List(point)
                YukS2Quark1 = YukS2Quark1List(point)
                YukS2Quark2 = YukS2Quark2List(point)
                YukS3Quark1 = YukS3Quark1List(point)
                YukS3Quark2 = YukS3Quark2List(point)
                m12squared = m12squaredList(point)
                vS = vSList(point)
                TypeOf2HDM = TypeOf2HDMList(point)

                ! Print out the current point in phase-space (debug mode only)
                if (debugModeOn) then
                    write (*,*) "MW: ", MW
                    write (*,*) "MZ: ", MZ
                    write (*,*) "SW: ", SW
                    write (*,*) "CW: ", CW
                    write (*,*) "EL: ", EL
                    write (*,*) "vev: ", (2D0*MW*SW/EL)
                    write (*,*) "ME: ", ME
                    write (*,*) "MM: ", MM
                    write (*,*) "ML: ", ML
                    write (*,*) "MU: ", MU
                    write (*,*) "MD: ", MD
                    write (*,*) "MS: ", MS
                    write (*,*) "MC: ", MC
                    write (*,*) "MB: ", MB
                    write (*,*) "MT: ", MT
                    write (*,*) "MH1: ", MH1
                    write (*,*) "MH2: ", MH2
                    write (*,*) "MH3: ", MH3
                    write (*,*) "MA0: ", MA0
                    write (*,*) "MHp: ", MHp
                    write (*,*) "alpha1: ", alpha1
                    write (*,*) "alpha2: ", alpha2
                    write (*,*) "alpha3: ", alpha3
                    write (*,*) "beta: ", beta
                    write (*,*) "CA1: ", CA1
                    write (*,*) "CA2: ", CA2
                    write (*,*) "CA3: ", CA3
                    write (*,*) "CB: ", CB
                    write (*,*) "YukS1Lep1: ", YukS1Lep1
                    write (*,*) "YukS1Lep2: ", YukS1Lep2
                    write (*,*) "YukS1Lep3: ", YukS1Lep3
                    write (*,*) "YukS2Lep1: ", YukS2Lep1
                    write (*,*) "YukS2Lep2: ", YukS2Lep2
                    write (*,*) "YukS3Lep1: ", YukS3Lep1
                    write (*,*) "YukS3Lep2: ", YukS3Lep2
                    write (*,*) "YukS1Quark1: ", YukS1Quark1
                    write (*,*) "YukS1Quark2: ", YukS1Quark2
                    write (*,*) "YukS1Quark3: ", YukS1Quark3
                    write (*,*) "YukS2Quark1: ", YukS2Quark1
                    write (*,*) "YukS2Quark2: ", YukS2Quark2
                    write (*,*) "YukS3Quark1: ", YukS3Quark1
                    write (*,*) "YukS3Quark2: ", YukS3Quark2
                    write (*,*) "m12squared: ", m12squared
                    write (*,*) "vS: ", vS
                    write (*,*) "2HDM Type: ", TypeOf2HDM
                    write (*,*) "InputScale: ", InputScale
                end if

                ! Calculate the square of the 2HDM input parameters
                MH22 = MH2**2
                MH32 = MH3**2
                MA02 = MA0**2
                MHp2 = MHp**2
                CA12 = CA1**2
                CA22 = CA2**2
                CA32 = CA3**2
                SA12 = SA1**2
                SA22 = SA2**2
                SA32 = SA3**2
                TB2 = TB**2
                SB2 = SB**2
                CB2 = CB**2

                ! Set the scalar couplings
RR11 = CA1*CA2
RR12 = CA2*SA1
RR13 = SA2
RR21 = -1.D0*CA3*SA1 - 1.D0*CA1*SA2*SA3
RR22 = CA1*CA3 - 1.D0*SA1*SA2*SA3
RR23 = CA2*SA3
RR31 = -1.D0*CA1*CA3*SA2 + SA1*SA3
RR32 = -1.D0*CA3*SA1*SA2 - 1.D0*CA1*SA3
RR33 = CA2*CA3

Lam1 = (0.25D0*EL2*((-1.D0*m12squared*SB)/CB + MH12*DBLE(RR11**INT(2.D0)) + MH22*DBLE(RR21**INT(2.D0)) + MH32*DBLE(RR31**INT(2.D0&
  &))))/(CB2*MW2*SW2)
Lam2 = (0.25D0*EL2*((-1.D0*CB*m12squared)/SB + MH12*DBLE(RR12**INT(2.D0)) + MH22*DBLE(RR22**INT(2.D0)) + MH32*DBLE(RR32**INT(2.D0&
  &))))/(MW2*SB2*SW2)
Lam3 = (0.25D0*EL2*(2.D0*MHp2 - (1.D0*m12squared)/(CB*SB) + (MH12*RR11*RR12 + MH22*RR21*RR22 + MH32*RR31*RR32)/(CB*SB)))/(MW2*SW2&
  &)
Lam4 = (0.25D0*EL2*(MA02 - 2.D0*MHp2 + m12squared/(CB*SB)))/(MW2*SW2)
Lam5 = (0.25D0*EL2*(-1.D0*MA02 + m12squared/(CB*SB)))/(MW2*SW2)
Lam6 = (MH12*DBLE(RR13**INT(2.D0)) + MH22*DBLE(RR23**INT(2.D0)) + MH32*DBLE(RR33**INT(2.D0)))*DBLE(vS**INT(-2.D0))
Lam7 = (0.5D0*EL*(MH12*RR11*RR13 + MH22*RR21*RR23 + MH32*RR31*RR33))/(CB*MW*SW*vS)
Lam8 = (0.5D0*EL*(MH12*RR12*RR13 + MH22*RR22*RR23 + MH32*RR32*RR33))/(MW*SB*SW*vS)

CS1S1S1f111 = -1.D0*RR13*(Lam7*RR11*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR12*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR11*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f112 = -1.D0*RR13*(Lam7*RR11*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR12*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR11*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f113 = -1.D0*RR13*(Lam7*RR11*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR12*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR11*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f121 = -1.D0*RR13*(Lam7*RR21*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR22*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR21*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f122 = -1.D0*RR13*(Lam7*RR21*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR22*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR21*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f123 = -1.D0*RR13*(Lam7*RR21*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR22*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR21*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f131 = -1.D0*RR13*(Lam7*RR31*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR32*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR31*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f132 = -1.D0*RR13*(Lam7*RR31*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR32*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR31*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f133 = -1.D0*RR13*(Lam7*RR31*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR32*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR31*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f211 = -1.D0*RR23*(Lam7*RR11*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR12*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR11*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f212 = -1.D0*RR23*(Lam7*RR11*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR12*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR11*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f213 = -1.D0*RR23*(Lam7*RR11*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR12*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR11*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f221 = -1.D0*RR23*(Lam7*RR21*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR22*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR21*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f222 = -1.D0*RR23*(Lam7*RR21*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR22*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR21*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f223 = -1.D0*RR23*(Lam7*RR21*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR22*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR21*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f231 = -1.D0*RR23*(Lam7*RR31*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR32*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR31*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f232 = -1.D0*RR23*(Lam7*RR31*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR32*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR31*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f233 = -1.D0*RR23*(Lam7*RR31*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR32*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR31*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f311 = -1.D0*RR33*(Lam7*RR11*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR12*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR11*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f312 = -1.D0*RR33*(Lam7*RR11*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR12*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR11*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f313 = -1.D0*RR33*(Lam7*RR11*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR12*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR13*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR12*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR13*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR11*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR11*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR13*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR12*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f321 = -1.D0*RR33*(Lam7*RR21*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR22*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR21*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f322 = -1.D0*RR33*(Lam7*RR21*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR22*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR21*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f323 = -1.D0*RR33*(Lam7*RR21*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR22*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR23*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR22*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR23*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR21*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR21*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR23*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR22*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS1S1S1f331 = -1.D0*RR33*(Lam7*RR31*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + Lam8*RR32*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR11*SW)/EL + (2.D0*Lam8*MW*RR12*SB*SW)/EL + 3.D0*Lam6*RR13*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR13*SW)/EL + RR11*vs) + RR31*((6.D0*CB*Lam1*MW*RR11*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR12*&
  &SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR13*SB*SW)/EL + RR12*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR1&
  &1*SW)/EL + (6.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS1S1S1f332 = -1.D0*RR33*(Lam7*RR31*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + Lam8*RR32*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR21*SW)/EL + (2.D0*Lam8*MW*RR22*SB*SW)/EL + 3.D0*Lam6*RR23*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR23*SW)/EL + RR21*vs) + RR31*((6.D0*CB*Lam1*MW*RR21*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR22*&
  &SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR23*SB*SW)/EL + RR22*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR2&
  &1*SW)/EL + (6.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS1S1S1f333 = -1.D0*RR33*(Lam7*RR31*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + Lam8*RR32*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR33*((&
  &2.D0*CB*Lam7*MW*RR31*SW)/EL + (2.D0*Lam8*MW*RR32*SB*SW)/EL + 3.D0*Lam6*RR33*vs)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR32*((2.D0&
  &*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam7*RR33*((2.D0*CB*MW*RR33*SW)/EL + RR31*vs) + RR31*((6.D0*CB*Lam1*MW*RR31*SW&
  &)/EL + (2.D0*(Lam3 + Lam4 + Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR31*((2.D0*CB*MW*RR32*&
  &SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + Lam8*RR33*((2.D0*MW*RR33*SB*SW)/EL + RR32*vs) + RR32*((2.D0*CB*(Lam3 + Lam4 + Lam5)*MW*RR3&
  &1*SW)/EL + (6.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))

CS2S2S1f111 = -1.D0*CB*(Lam5*SB*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + CB*((2.D0*CB*Lam1*MW*RR11*SW)/EL + (2.D0*(L&
  &am3 + Lam4 - 1.D0*Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*SB*(CB*Lam5*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/&
  &EL) + SB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR11*SW)/EL + (2.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS2S2S1f112 = -1.D0*CB*(Lam5*SB*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + CB*((2.D0*CB*Lam1*MW*RR21*SW)/EL + (2.D0*(L&
  &am3 + Lam4 - 1.D0*Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*SB*(CB*Lam5*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/&
  &EL) + SB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR21*SW)/EL + (2.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS2S2S1f113 = -1.D0*CB*(Lam5*SB*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + CB*((2.D0*CB*Lam1*MW*RR31*SW)/EL + (2.D0*(L&
  &am3 + Lam4 - 1.D0*Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*SB*(CB*Lam5*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/&
  &EL) + SB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR31*SW)/EL + (2.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS2S2S1f121 = -1.D0*CB*(CB*Lam5*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) - 1.D0*SB*((2.D0*CB*Lam1*MW*RR11*SW)/EL + (2.&
  &D0*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*SB*(-1.D0*Lam5*SB*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*R&
  &R11*SB*SW)/EL) + CB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR11*SW)/EL + (2.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS2S2S1f122 = -1.D0*CB*(CB*Lam5*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) - 1.D0*SB*((2.D0*CB*Lam1*MW*RR21*SW)/EL + (2.&
  &D0*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*SB*(-1.D0*Lam5*SB*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*R&
  &R21*SB*SW)/EL) + CB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR21*SW)/EL + (2.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS2S2S1f123 = -1.D0*CB*(CB*Lam5*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) - 1.D0*SB*((2.D0*CB*Lam1*MW*RR31*SW)/EL + (2.&
  &D0*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*SB*(-1.D0*Lam5*SB*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*R&
  &R31*SB*SW)/EL) + CB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR31*SW)/EL + (2.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS2S2S1f211 = SB*(Lam5*SB*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + CB*((2.D0*CB*Lam1*MW*RR11*SW)/EL + (2.D0*(Lam3 + &
  &Lam4 - 1.D0*Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*CB*(CB*Lam5*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) + &
  &SB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR11*SW)/EL + (2.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS2S2S1f212 = SB*(Lam5*SB*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + CB*((2.D0*CB*Lam1*MW*RR21*SW)/EL + (2.D0*(Lam3 + &
  &Lam4 - 1.D0*Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*CB*(CB*Lam5*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) + &
  &SB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR21*SW)/EL + (2.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS2S2S1f213 = SB*(Lam5*SB*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + CB*((2.D0*CB*Lam1*MW*RR31*SW)/EL + (2.D0*(Lam3 + &
  &Lam4 - 1.D0*Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*CB*(CB*Lam5*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) + &
  &SB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR31*SW)/EL + (2.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))
CS2S2S1f221 = SB*(CB*Lam5*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB*SW)/EL) - 1.D0*SB*((2.D0*CB*Lam1*MW*RR11*SW)/EL + (2.D0*(La&
  &m3 + Lam4 - 1.D0*Lam5)*MW*RR12*SB*SW)/EL + Lam7*RR13*vs)) - 1.D0*CB*(-1.D0*Lam5*SB*((2.D0*CB*MW*RR12*SW)/EL + (2.D0*MW*RR11*SB&
  &*SW)/EL) + CB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR11*SW)/EL + (2.D0*Lam2*MW*RR12*SB*SW)/EL + Lam8*RR13*vs))
CS2S2S1f222 = SB*(CB*Lam5*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB*SW)/EL) - 1.D0*SB*((2.D0*CB*Lam1*MW*RR21*SW)/EL + (2.D0*(La&
  &m3 + Lam4 - 1.D0*Lam5)*MW*RR22*SB*SW)/EL + Lam7*RR23*vs)) - 1.D0*CB*(-1.D0*Lam5*SB*((2.D0*CB*MW*RR22*SW)/EL + (2.D0*MW*RR21*SB&
  &*SW)/EL) + CB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR21*SW)/EL + (2.D0*Lam2*MW*RR22*SB*SW)/EL + Lam8*RR23*vs))
CS2S2S1f223 = SB*(CB*Lam5*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB*SW)/EL) - 1.D0*SB*((2.D0*CB*Lam1*MW*RR31*SW)/EL + (2.D0*(La&
  &m3 + Lam4 - 1.D0*Lam5)*MW*RR32*SB*SW)/EL + Lam7*RR33*vs)) - 1.D0*CB*(-1.D0*Lam5*SB*((2.D0*CB*MW*RR32*SW)/EL + (2.D0*MW*RR31*SB&
  &*SW)/EL) + CB*((2.D0*CB*(Lam3 + Lam4 - 1.D0*Lam5)*MW*RR31*SW)/EL + (2.D0*Lam2*MW*RR32*SB*SW)/EL + Lam8*RR33*vs))

CS1S3S3f111 = 0.5D0*(-1.D0*RR12*(CB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + SB*((2.D0*CB2*(Lam4 + La&
  &m5)*MW*SW)/EL + (4.D0*Lam2*MW*SB2*SW)/EL)) - 1.D0*RR11*(SB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) &
  &+ CB*((4.D0*CB2*Lam1*MW*SW)/EL + (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR13*(CB2*Lam7 + Lam8*SB2)*vs)
CS1S3S3f112 = 0.5D0*(-1.D0*RR12*(SB*((4.D0*CB*Lam2*MW*SB*SW)/EL - (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + CB*((2.D0*CB2*(Lam4 + La&
  &m5)*MW*SW)/EL - (4.D0*Lam3*MW*SB2*SW)/EL)) - 1.D0*RR11*(CB*((-4.D0*CB*Lam1*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL)&
  & + SB*((4.D0*CB2*Lam3*MW*SW)/EL - (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR13*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*vs)
CS1S3S3f121 = 0.5D0*(-1.D0*RR12*(-1.D0*SB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + CB*((2.D0*CB2*(Lam&
  &4 + Lam5)*MW*SW)/EL + (4.D0*Lam2*MW*SB2*SW)/EL)) - 1.D0*RR11*(CB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW&
  &)/EL) - 1.D0*SB*((4.D0*CB2*Lam1*MW*SW)/EL + (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR13*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*vs&
  &)
CS1S3S3f122 = 0.5D0*(-1.D0*RR12*(CB*((4.D0*CB*Lam2*MW*SB*SW)/EL - (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) - 1.D0*SB*((2.D0*CB2*(Lam4&
  & + Lam5)*MW*SW)/EL - (4.D0*Lam3*MW*SB2*SW)/EL)) - 1.D0*RR11*(-1.D0*SB*((-4.D0*CB*Lam1*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW&
  &*SB*SW)/EL) + CB*((4.D0*CB2*Lam3*MW*SW)/EL - (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR13*(CB2*Lam8 + Lam7*SB2)*vs)
CS1S3S3f211 = 0.5D0*(-1.D0*RR22*(CB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + SB*((2.D0*CB2*(Lam4 + La&
  &m5)*MW*SW)/EL + (4.D0*Lam2*MW*SB2*SW)/EL)) - 1.D0*RR21*(SB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) &
  &+ CB*((4.D0*CB2*Lam1*MW*SW)/EL + (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR23*(CB2*Lam7 + Lam8*SB2)*vs)
CS1S3S3f212 = 0.5D0*(-1.D0*RR22*(SB*((4.D0*CB*Lam2*MW*SB*SW)/EL - (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + CB*((2.D0*CB2*(Lam4 + La&
  &m5)*MW*SW)/EL - (4.D0*Lam3*MW*SB2*SW)/EL)) - 1.D0*RR21*(CB*((-4.D0*CB*Lam1*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL)&
  & + SB*((4.D0*CB2*Lam3*MW*SW)/EL - (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR23*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*vs)
CS1S3S3f221 = 0.5D0*(-1.D0*RR22*(-1.D0*SB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + CB*((2.D0*CB2*(Lam&
  &4 + Lam5)*MW*SW)/EL + (4.D0*Lam2*MW*SB2*SW)/EL)) - 1.D0*RR21*(CB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW&
  &)/EL) - 1.D0*SB*((4.D0*CB2*Lam1*MW*SW)/EL + (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR23*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*vs&
  &)
CS1S3S3f222 = 0.5D0*(-1.D0*RR22*(CB*((4.D0*CB*Lam2*MW*SB*SW)/EL - (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) - 1.D0*SB*((2.D0*CB2*(Lam4&
  & + Lam5)*MW*SW)/EL - (4.D0*Lam3*MW*SB2*SW)/EL)) - 1.D0*RR21*(-1.D0*SB*((-4.D0*CB*Lam1*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW&
  &*SB*SW)/EL) + CB*((4.D0*CB2*Lam3*MW*SW)/EL - (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR23*(CB2*Lam8 + Lam7*SB2)*vs)
CS1S3S3f311 = 0.5D0*(-1.D0*RR32*(CB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + SB*((2.D0*CB2*(Lam4 + La&
  &m5)*MW*SW)/EL + (4.D0*Lam2*MW*SB2*SW)/EL)) - 1.D0*RR31*(SB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) &
  &+ CB*((4.D0*CB2*Lam1*MW*SW)/EL + (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR33*(CB2*Lam7 + Lam8*SB2)*vs)
CS1S3S3f312 = 0.5D0*(-1.D0*RR32*(SB*((4.D0*CB*Lam2*MW*SB*SW)/EL - (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + CB*((2.D0*CB2*(Lam4 + La&
  &m5)*MW*SW)/EL - (4.D0*Lam3*MW*SB2*SW)/EL)) - 1.D0*RR31*(CB*((-4.D0*CB*Lam1*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL)&
  & + SB*((4.D0*CB2*Lam3*MW*SW)/EL - (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*vs)
CS1S3S3f321 = 0.5D0*(-1.D0*RR32*(-1.D0*SB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) + CB*((2.D0*CB2*(Lam&
  &4 + Lam5)*MW*SW)/EL + (4.D0*Lam2*MW*SB2*SW)/EL)) - 1.D0*RR31*(CB*((4.D0*CB*Lam3*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW&
  &)/EL) - 1.D0*SB*((4.D0*CB2*Lam1*MW*SW)/EL + (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*vs&
  &)
CS1S3S3f322 = 0.5D0*(-1.D0*RR32*(CB*((4.D0*CB*Lam2*MW*SB*SW)/EL - (2.D0*CB*(Lam4 + Lam5)*MW*SB*SW)/EL) - 1.D0*SB*((2.D0*CB2*(Lam4&
  & + Lam5)*MW*SW)/EL - (4.D0*Lam3*MW*SB2*SW)/EL)) - 1.D0*RR31*(-1.D0*SB*((-4.D0*CB*Lam1*MW*SB*SW)/EL + (2.D0*CB*(Lam4 + Lam5)*MW&
  &*SB*SW)/EL) + CB*((4.D0*CB2*Lam3*MW*SW)/EL - (2.D0*(Lam4 + Lam5)*MW*SB2*SW)/EL)) - 2.D0*RR33*(CB2*Lam8 + Lam7*SB2)*vs)

CS1S1S1S1f1111 = -1.D0*RR13*(2.D0*Lam7*RR13*DBLE(RR11**INT(2.D0)) + 2.D0*Lam8*RR13*DBLE(RR12**INT(2.D0)) + RR13*(Lam7*DBLE(RR11**&
  &INT(2.D0)) + Lam8*DBLE(RR12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*DBLE(R&
  &R12**INT(2.D0)) + 2.D0*Lam7*RR11*DBLE(RR13**INT(2.D0)) + RR11*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR1&
  &2**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*DBLE(RR11**INT(2.D0)) + 2.D0*Lam8*RR1&
  &2*DBLE(RR13**INT(2.D0)) + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13*&
  &*INT(2.D0))))
CS1S1S1S1f1112 = -1.D0*RR13*(Lam7*RR11*(RR13*RR21 + RR11*RR23) + Lam8*RR12*(RR13*RR22 + RR12*RR23) + RR13*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR21 + RR11*RR22) + Lam7*RR13*(RR13*RR21 + RR11&
  &*RR23) + RR11*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR13*(RR13*RR22 + RR12*RR23) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f1113 = -1.D0*RR13*(Lam7*RR11*(RR13*RR31 + RR11*RR33) + Lam8*RR12*(RR13*RR32 + RR12*RR33) + RR13*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR31 + RR11*RR32) + Lam7*RR13*(RR13*RR31 + RR11&
  &*RR33) + RR11*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR13*(RR13*RR32 + RR12*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f1121 = -1.D0*RR13*(Lam7*RR11*(RR13*RR21 + RR11*RR23) + Lam8*RR12*(RR13*RR22 + RR12*RR23) + RR13*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR21 + RR11*RR22) + Lam7*RR13*(RR13*RR21 + RR11&
  &*RR23) + RR11*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR13*(RR13*RR22 + RR12*RR23) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f1122 = -1.D0*RR13*(2.D0*Lam7*RR11*RR21*RR23 + 2.D0*Lam8*RR12*RR22*RR23 + RR13*(Lam7*DBLE(RR21**INT(2.D0)) + Lam8*DBLE(R&
  &R22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*RR21*RR22 + 2.D0*Lam7*RR13*RR2&
  &1*RR23 + RR11*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR22**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1&
  &.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR21*RR22 + 2.D0*Lam8*RR13*RR22*RR23 + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS1S1S1S1f1123 = -1.D0*RR13*(Lam7*RR11*(RR23*RR31 + RR21*RR33) + Lam8*RR12*(RR23*RR32 + RR22*RR33) + RR13*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*(RR22*RR31 + RR21*RR32) + Lam7*RR13*(RR23*RR31 + RR21&
  &*RR33) + RR11*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR13*(RR23*RR32 + RR22*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f1131 = -1.D0*RR13*(Lam7*RR11*(RR13*RR31 + RR11*RR33) + Lam8*RR12*(RR13*RR32 + RR12*RR33) + RR13*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR31 + RR11*RR32) + Lam7*RR13*(RR13*RR31 + RR11&
  &*RR33) + RR11*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR13*(RR13*RR32 + RR12*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f1132 = -1.D0*RR13*(Lam7*RR11*(RR23*RR31 + RR21*RR33) + Lam8*RR12*(RR23*RR32 + RR22*RR33) + RR13*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR12*(RR22*RR31 + RR21*RR32) + Lam7*RR13*(RR23*RR31 + RR21&
  &*RR33) + RR11*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR13*(RR23*RR32 + RR22*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f1133 = -1.D0*RR13*(2.D0*Lam7*RR11*RR31*RR33 + 2.D0*Lam8*RR12*RR32*RR33 + RR13*(Lam7*DBLE(RR31**INT(2.D0)) + Lam8*DBLE(R&
  &R32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*RR31*RR32 + 2.D0*Lam7*RR13*RR3&
  &1*RR33 + RR11*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR32**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1&
  &.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR31*RR32 + 2.D0*Lam8*RR13*RR32*RR33 + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS1S1S1S1f1211 = -1.D0*RR13*(2.D0*Lam7*RR11*RR13*RR21 + 2.D0*Lam8*RR12*RR13*RR22 + RR23*(Lam7*DBLE(RR11**INT(2.D0)) + Lam8*DBLE(R&
  &R12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR22 + 2.D0*Lam7*RR11*RR1&
  &3*RR23 + RR21*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR12**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1&
  &.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR21 + 2.D0*Lam8*RR12*RR13*RR23 + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS1S1S1S1f1212 = -1.D0*RR13*(Lam7*RR21*(RR13*RR21 + RR11*RR23) + Lam8*RR22*(RR13*RR22 + RR12*RR23) + RR23*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR21 + RR11*RR22) + Lam7*RR23*(RR13*RR21 + RR11&
  &*RR23) + RR21*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR23*(RR13*RR22 + RR12*RR23) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f1213 = -1.D0*RR13*(Lam7*RR21*(RR13*RR31 + RR11*RR33) + Lam8*RR22*(RR13*RR32 + RR12*RR33) + RR23*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR31 + RR11*RR32) + Lam7*RR23*(RR13*RR31 + RR11&
  &*RR33) + RR21*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR23*(RR13*RR32 + RR12*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f1221 = -1.D0*RR13*(Lam7*RR21*(RR13*RR21 + RR11*RR23) + Lam8*RR22*(RR13*RR22 + RR12*RR23) + RR23*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR21 + RR11*RR22) + Lam7*RR23*(RR13*RR21 + RR11&
  &*RR23) + RR21*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR23*(RR13*RR22 + RR12*RR23) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f1222 = -1.D0*RR13*(2.D0*Lam7*RR23*DBLE(RR21**INT(2.D0)) + 2.D0*Lam8*RR23*DBLE(RR22**INT(2.D0)) + RR23*(Lam7*DBLE(RR21**&
  &INT(2.D0)) + Lam8*DBLE(RR22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*DBLE(R&
  &R22**INT(2.D0)) + 2.D0*Lam7*RR21*DBLE(RR23**INT(2.D0)) + RR21*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR2&
  &2**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR22*DBLE(RR21**INT(2.D0)) + 2.D0*Lam8*RR2&
  &2*DBLE(RR23**INT(2.D0)) + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23*&
  &*INT(2.D0))))
CS1S1S1S1f1223 = -1.D0*RR13*(Lam7*RR21*(RR23*RR31 + RR21*RR33) + Lam8*RR22*(RR23*RR32 + RR22*RR33) + RR23*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*(RR22*RR31 + RR21*RR32) + Lam7*RR23*(RR23*RR31 + RR21&
  &*RR33) + RR21*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR23*(RR23*RR32 + RR22*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f1231 = -1.D0*RR13*(Lam7*RR21*(RR13*RR31 + RR11*RR33) + Lam8*RR22*(RR13*RR32 + RR12*RR33) + RR23*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR31 + RR11*RR32) + Lam7*RR23*(RR13*RR31 + RR11&
  &*RR33) + RR21*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR23*(RR13*RR32 + RR12*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f1232 = -1.D0*RR13*(Lam7*RR21*(RR23*RR31 + RR21*RR33) + Lam8*RR22*(RR23*RR32 + RR22*RR33) + RR23*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR22*(RR22*RR31 + RR21*RR32) + Lam7*RR23*(RR23*RR31 + RR21&
  &*RR33) + RR21*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR23*(RR23*RR32 + RR22*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f1233 = -1.D0*RR13*(2.D0*Lam7*RR21*RR31*RR33 + 2.D0*Lam8*RR22*RR32*RR33 + RR23*(Lam7*DBLE(RR31**INT(2.D0)) + Lam8*DBLE(R&
  &R32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR22*RR31*RR32 + 2.D0*Lam7*RR23*RR3&
  &1*RR33 + RR21*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR32**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1&
  &.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR31*RR32 + 2.D0*Lam8*RR23*RR32*RR33 + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS1S1S1S1f1311 = -1.D0*RR13*(2.D0*Lam7*RR11*RR13*RR31 + 2.D0*Lam8*RR12*RR13*RR32 + RR33*(Lam7*DBLE(RR11**INT(2.D0)) + Lam8*DBLE(R&
  &R12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR32 + 2.D0*Lam7*RR11*RR1&
  &3*RR33 + RR31*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR12**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1&
  &.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR31 + 2.D0*Lam8*RR12*RR13*RR33 + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS1S1S1S1f1312 = -1.D0*RR11*((3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)*RR31 + (Lam3 + Lam4 + Lam5)*&
  &(RR12*RR21 + RR11*RR22)*RR32 + Lam7*(RR13*RR21 + RR11*RR23)*RR33) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*(RR12*RR21 + RR11*RR22)*RR&
  &31 + ((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23)*RR32 + Lam8*(RR13*RR22 + RR12*RR23)*RR33) - 1.D0*&
  &RR13*(Lam7*(RR13*RR21 + RR11*RR23)*RR31 + Lam8*(RR13*RR22 + RR12*RR23)*RR32 + (Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23)*RR33)
CS1S1S1S1f1313 = -1.D0*RR13*(Lam7*RR31*(RR13*RR31 + RR11*RR33) + Lam8*RR32*(RR13*RR32 + RR12*RR33) + RR33*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR32*(RR12*RR31 + RR11*RR32) + Lam7*RR33*(RR13*RR31 + RR11&
  &*RR33) + RR31*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR33*(RR13*RR32 + RR12*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f1321 = -1.D0*RR11*((3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)*RR31 + (Lam3 + Lam4 + Lam5)*&
  &(RR12*RR21 + RR11*RR22)*RR32 + Lam7*(RR13*RR21 + RR11*RR23)*RR33) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*(RR12*RR21 + RR11*RR22)*RR&
  &31 + ((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23)*RR32 + Lam8*(RR13*RR22 + RR12*RR23)*RR33) - 1.D0*&
  &RR13*(Lam7*(RR13*RR21 + RR11*RR23)*RR31 + Lam8*(RR13*RR22 + RR12*RR23)*RR32 + (Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23)*RR33)
CS1S1S1S1f1322 = -1.D0*RR13*(2.D0*Lam7*RR21*RR23*RR31 + 2.D0*Lam8*RR22*RR23*RR32 + RR33*(Lam7*DBLE(RR21**INT(2.D0)) + Lam8*DBLE(R&
  &R22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR22*RR32 + 2.D0*Lam7*RR21*RR2&
  &3*RR33 + RR31*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR22**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1&
  &.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR22*RR31 + 2.D0*Lam8*RR22*RR23*RR33 + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS1S1S1S1f1323 = -1.D0*RR13*(Lam7*RR31*(RR23*RR31 + RR21*RR33) + Lam8*RR32*(RR23*RR32 + RR22*RR33) + RR33*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR32*(RR22*RR31 + RR21*RR32) + Lam7*RR33*(RR23*RR31 + RR21&
  &*RR33) + RR31*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR33*(RR23*RR32 + RR22*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f1331 = -1.D0*RR13*(Lam7*RR31*(RR13*RR31 + RR11*RR33) + Lam8*RR32*(RR13*RR32 + RR12*RR33) + RR33*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR32*(RR12*RR31 + RR11*RR32) + Lam7*RR33*(RR13*RR31 + RR11&
  &*RR33) + RR31*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR33*(RR13*RR32 + RR12*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f1332 = -1.D0*RR13*(Lam7*RR31*(RR23*RR31 + RR21*RR33) + Lam8*RR32*(RR23*RR32 + RR22*RR33) + RR33*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR11*((Lam3 + Lam4 + Lam5)*RR32*(RR22*RR31 + RR21*RR32) + Lam7*RR33*(RR23*RR31 + RR21&
  &*RR33) + RR31*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR12*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR33*(RR23*RR32 + RR22*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f1333 = -1.D0*RR13*(2.D0*Lam7*RR33*DBLE(RR31**INT(2.D0)) + 2.D0*Lam8*RR33*DBLE(RR32**INT(2.D0)) + RR33*(Lam7*DBLE(RR31**&
  &INT(2.D0)) + Lam8*DBLE(RR32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR11*(2.D0*(Lam3 + Lam4 + Lam5)*RR31*DBLE(R&
  &R32**INT(2.D0)) + 2.D0*Lam7*RR31*DBLE(RR33**INT(2.D0)) + RR31*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR3&
  &2**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1.D0*RR12*(2.D0*(Lam3 + Lam4 + Lam5)*RR32*DBLE(RR31**INT(2.D0)) + 2.D0*Lam8*RR3&
  &2*DBLE(RR33**INT(2.D0)) + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33*&
  &*INT(2.D0))))
CS1S1S1S1f2111 = -1.D0*RR23*(2.D0*Lam7*RR13*DBLE(RR11**INT(2.D0)) + 2.D0*Lam8*RR13*DBLE(RR12**INT(2.D0)) + RR13*(Lam7*DBLE(RR11**&
  &INT(2.D0)) + Lam8*DBLE(RR12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*DBLE(R&
  &R12**INT(2.D0)) + 2.D0*Lam7*RR11*DBLE(RR13**INT(2.D0)) + RR11*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR1&
  &2**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*DBLE(RR11**INT(2.D0)) + 2.D0*Lam8*RR1&
  &2*DBLE(RR13**INT(2.D0)) + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13*&
  &*INT(2.D0))))
CS1S1S1S1f2112 = -1.D0*RR23*(Lam7*RR11*(RR13*RR21 + RR11*RR23) + Lam8*RR12*(RR13*RR22 + RR12*RR23) + RR13*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR21 + RR11*RR22) + Lam7*RR13*(RR13*RR21 + RR11&
  &*RR23) + RR11*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR13*(RR13*RR22 + RR12*RR23) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f2113 = -1.D0*RR23*(Lam7*RR11*(RR13*RR31 + RR11*RR33) + Lam8*RR12*(RR13*RR32 + RR12*RR33) + RR13*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR31 + RR11*RR32) + Lam7*RR13*(RR13*RR31 + RR11&
  &*RR33) + RR11*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR13*(RR13*RR32 + RR12*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f2121 = -1.D0*RR23*(Lam7*RR11*(RR13*RR21 + RR11*RR23) + Lam8*RR12*(RR13*RR22 + RR12*RR23) + RR13*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR21 + RR11*RR22) + Lam7*RR13*(RR13*RR21 + RR11&
  &*RR23) + RR11*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR13*(RR13*RR22 + RR12*RR23) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f2122 = -1.D0*RR23*(2.D0*Lam7*RR11*RR21*RR23 + 2.D0*Lam8*RR12*RR22*RR23 + RR13*(Lam7*DBLE(RR21**INT(2.D0)) + Lam8*DBLE(R&
  &R22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*RR21*RR22 + 2.D0*Lam7*RR13*RR2&
  &1*RR23 + RR11*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR22**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1&
  &.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR21*RR22 + 2.D0*Lam8*RR13*RR22*RR23 + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS1S1S1S1f2123 = -1.D0*RR23*(Lam7*RR11*(RR23*RR31 + RR21*RR33) + Lam8*RR12*(RR23*RR32 + RR22*RR33) + RR13*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*(RR22*RR31 + RR21*RR32) + Lam7*RR13*(RR23*RR31 + RR21&
  &*RR33) + RR11*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR13*(RR23*RR32 + RR22*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f2131 = -1.D0*RR23*(Lam7*RR11*(RR13*RR31 + RR11*RR33) + Lam8*RR12*(RR13*RR32 + RR12*RR33) + RR13*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR31 + RR11*RR32) + Lam7*RR13*(RR13*RR31 + RR11&
  &*RR33) + RR11*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR13*(RR13*RR32 + RR12*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f2132 = -1.D0*RR23*(Lam7*RR11*(RR23*RR31 + RR21*RR33) + Lam8*RR12*(RR23*RR32 + RR22*RR33) + RR13*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR12*(RR22*RR31 + RR21*RR32) + Lam7*RR13*(RR23*RR31 + RR21&
  &*RR33) + RR11*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR13*(RR23*RR32 + RR22*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f2133 = -1.D0*RR23*(2.D0*Lam7*RR11*RR31*RR33 + 2.D0*Lam8*RR12*RR32*RR33 + RR13*(Lam7*DBLE(RR31**INT(2.D0)) + Lam8*DBLE(R&
  &R32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*RR31*RR32 + 2.D0*Lam7*RR13*RR3&
  &1*RR33 + RR11*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR32**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1&
  &.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR31*RR32 + 2.D0*Lam8*RR13*RR32*RR33 + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS1S1S1S1f2211 = -1.D0*RR23*(2.D0*Lam7*RR11*RR13*RR21 + 2.D0*Lam8*RR12*RR13*RR22 + RR23*(Lam7*DBLE(RR11**INT(2.D0)) + Lam8*DBLE(R&
  &R12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR22 + 2.D0*Lam7*RR11*RR1&
  &3*RR23 + RR21*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR12**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1&
  &.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR21 + 2.D0*Lam8*RR12*RR13*RR23 + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS1S1S1S1f2212 = -1.D0*RR23*(Lam7*RR21*(RR13*RR21 + RR11*RR23) + Lam8*RR22*(RR13*RR22 + RR12*RR23) + RR23*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR21 + RR11*RR22) + Lam7*RR23*(RR13*RR21 + RR11&
  &*RR23) + RR21*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR23*(RR13*RR22 + RR12*RR23) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f2213 = -1.D0*RR23*(Lam7*RR21*(RR13*RR31 + RR11*RR33) + Lam8*RR22*(RR13*RR32 + RR12*RR33) + RR23*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR31 + RR11*RR32) + Lam7*RR23*(RR13*RR31 + RR11&
  &*RR33) + RR21*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR23*(RR13*RR32 + RR12*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f2221 = -1.D0*RR23*(Lam7*RR21*(RR13*RR21 + RR11*RR23) + Lam8*RR22*(RR13*RR22 + RR12*RR23) + RR23*(Lam7*RR11*RR21 + Lam8*&
  &RR12*RR22 + 3.D0*Lam6*RR13*RR23)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR21 + RR11*RR22) + Lam7*RR23*(RR13*RR21 + RR11&
  &*RR23) + RR21*(3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR21 + RR11*RR22) + Lam8*RR23*(RR13*RR22 + RR12*RR23) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam&
  &8*RR13*RR23))
CS1S1S1S1f2222 = -1.D0*RR23*(2.D0*Lam7*RR23*DBLE(RR21**INT(2.D0)) + 2.D0*Lam8*RR23*DBLE(RR22**INT(2.D0)) + RR23*(Lam7*DBLE(RR21**&
  &INT(2.D0)) + Lam8*DBLE(RR22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*DBLE(R&
  &R22**INT(2.D0)) + 2.D0*Lam7*RR21*DBLE(RR23**INT(2.D0)) + RR21*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR2&
  &2**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR22*DBLE(RR21**INT(2.D0)) + 2.D0*Lam8*RR2&
  &2*DBLE(RR23**INT(2.D0)) + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23*&
  &*INT(2.D0))))
CS1S1S1S1f2223 = -1.D0*RR23*(Lam7*RR21*(RR23*RR31 + RR21*RR33) + Lam8*RR22*(RR23*RR32 + RR22*RR33) + RR23*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*(RR22*RR31 + RR21*RR32) + Lam7*RR23*(RR23*RR31 + RR21&
  &*RR33) + RR21*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR23*(RR23*RR32 + RR22*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f2231 = -1.D0*RR23*(Lam7*RR21*(RR13*RR31 + RR11*RR33) + Lam8*RR22*(RR13*RR32 + RR12*RR33) + RR23*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR31 + RR11*RR32) + Lam7*RR23*(RR13*RR31 + RR11&
  &*RR33) + RR21*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR23*(RR13*RR32 + RR12*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f2232 = -1.D0*RR23*(Lam7*RR21*(RR23*RR31 + RR21*RR33) + Lam8*RR22*(RR23*RR32 + RR22*RR33) + RR23*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR22*(RR22*RR31 + RR21*RR32) + Lam7*RR23*(RR23*RR31 + RR21&
  &*RR33) + RR21*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR23*(RR23*RR32 + RR22*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f2233 = -1.D0*RR23*(2.D0*Lam7*RR21*RR31*RR33 + 2.D0*Lam8*RR22*RR32*RR33 + RR23*(Lam7*DBLE(RR31**INT(2.D0)) + Lam8*DBLE(R&
  &R32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR22*RR31*RR32 + 2.D0*Lam7*RR23*RR3&
  &1*RR33 + RR21*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR32**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1&
  &.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR31*RR32 + 2.D0*Lam8*RR23*RR32*RR33 + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS1S1S1S1f2311 = -1.D0*RR23*(2.D0*Lam7*RR11*RR13*RR31 + 2.D0*Lam8*RR12*RR13*RR32 + RR33*(Lam7*DBLE(RR11**INT(2.D0)) + Lam8*DBLE(R&
  &R12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR32 + 2.D0*Lam7*RR11*RR1&
  &3*RR33 + RR31*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR12**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1&
  &.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR31 + 2.D0*Lam8*RR12*RR13*RR33 + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS1S1S1S1f2312 = -1.D0*RR21*((3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)*RR31 + (Lam3 + Lam4 + Lam5)*&
  &(RR12*RR21 + RR11*RR22)*RR32 + Lam7*(RR13*RR21 + RR11*RR23)*RR33) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*(RR12*RR21 + RR11*RR22)*RR&
  &31 + ((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23)*RR32 + Lam8*(RR13*RR22 + RR12*RR23)*RR33) - 1.D0*&
  &RR23*(Lam7*(RR13*RR21 + RR11*RR23)*RR31 + Lam8*(RR13*RR22 + RR12*RR23)*RR32 + (Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23)*RR33)
CS1S1S1S1f2313 = -1.D0*RR23*(Lam7*RR31*(RR13*RR31 + RR11*RR33) + Lam8*RR32*(RR13*RR32 + RR12*RR33) + RR33*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR32*(RR12*RR31 + RR11*RR32) + Lam7*RR33*(RR13*RR31 + RR11&
  &*RR33) + RR31*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR33*(RR13*RR32 + RR12*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f2321 = -1.D0*RR21*((3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)*RR31 + (Lam3 + Lam4 + Lam5)*&
  &(RR12*RR21 + RR11*RR22)*RR32 + Lam7*(RR13*RR21 + RR11*RR23)*RR33) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*(RR12*RR21 + RR11*RR22)*RR&
  &31 + ((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23)*RR32 + Lam8*(RR13*RR22 + RR12*RR23)*RR33) - 1.D0*&
  &RR23*(Lam7*(RR13*RR21 + RR11*RR23)*RR31 + Lam8*(RR13*RR22 + RR12*RR23)*RR32 + (Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23)*RR33)
CS1S1S1S1f2322 = -1.D0*RR23*(2.D0*Lam7*RR21*RR23*RR31 + 2.D0*Lam8*RR22*RR23*RR32 + RR33*(Lam7*DBLE(RR21**INT(2.D0)) + Lam8*DBLE(R&
  &R22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR22*RR32 + 2.D0*Lam7*RR21*RR2&
  &3*RR33 + RR31*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR22**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1&
  &.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR22*RR31 + 2.D0*Lam8*RR22*RR23*RR33 + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS1S1S1S1f2323 = -1.D0*RR23*(Lam7*RR31*(RR23*RR31 + RR21*RR33) + Lam8*RR32*(RR23*RR32 + RR22*RR33) + RR33*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR32*(RR22*RR31 + RR21*RR32) + Lam7*RR33*(RR23*RR31 + RR21&
  &*RR33) + RR31*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR33*(RR23*RR32 + RR22*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f2331 = -1.D0*RR23*(Lam7*RR31*(RR13*RR31 + RR11*RR33) + Lam8*RR32*(RR13*RR32 + RR12*RR33) + RR33*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR32*(RR12*RR31 + RR11*RR32) + Lam7*RR33*(RR13*RR31 + RR11&
  &*RR33) + RR31*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR33*(RR13*RR32 + RR12*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f2332 = -1.D0*RR23*(Lam7*RR31*(RR23*RR31 + RR21*RR33) + Lam8*RR32*(RR23*RR32 + RR22*RR33) + RR33*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR21*((Lam3 + Lam4 + Lam5)*RR32*(RR22*RR31 + RR21*RR32) + Lam7*RR33*(RR23*RR31 + RR21&
  &*RR33) + RR31*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR22*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR33*(RR23*RR32 + RR22*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f2333 = -1.D0*RR23*(2.D0*Lam7*RR33*DBLE(RR31**INT(2.D0)) + 2.D0*Lam8*RR33*DBLE(RR32**INT(2.D0)) + RR33*(Lam7*DBLE(RR31**&
  &INT(2.D0)) + Lam8*DBLE(RR32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR21*(2.D0*(Lam3 + Lam4 + Lam5)*RR31*DBLE(R&
  &R32**INT(2.D0)) + 2.D0*Lam7*RR31*DBLE(RR33**INT(2.D0)) + RR31*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR3&
  &2**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1.D0*RR22*(2.D0*(Lam3 + Lam4 + Lam5)*RR32*DBLE(RR31**INT(2.D0)) + 2.D0*Lam8*RR3&
  &2*DBLE(RR33**INT(2.D0)) + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33*&
  &*INT(2.D0))))
CS1S1S1S1f3111 = -1.D0*RR33*(2.D0*Lam7*RR13*DBLE(RR11**INT(2.D0)) + 2.D0*Lam8*RR13*DBLE(RR12**INT(2.D0)) + RR13*(Lam7*DBLE(RR11**&
  &INT(2.D0)) + Lam8*DBLE(RR12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*DBLE(R&
  &R12**INT(2.D0)) + 2.D0*Lam7*RR11*DBLE(RR13**INT(2.D0)) + RR11*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR1&
  &2**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*DBLE(RR11**INT(2.D0)) + 2.D0*Lam8*RR1&
  &2*DBLE(RR13**INT(2.D0)) + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13*&
  &*INT(2.D0))))
CS1S1S1S1f3112 = -1.D0*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR21 + RR11*RR22) + Lam7*RR13*(RR13*RR21 + RR11*RR23) + RR11*(3.D0*Lam1*R&
  &R11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23))*RR31 - 1.D0*((Lam3 + Lam4 + Lam5)*RR11*(RR12*RR21 + RR11*RR22) + &
  &Lam8*RR13*(RR13*RR22 + RR12*RR23) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23))*RR32 - 1.D0*&
  &(Lam7*RR11*(RR13*RR21 + RR11*RR23) + Lam8*RR12*(RR13*RR22 + RR12*RR23) + RR13*(Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23))*RR33
CS1S1S1S1f3113 = -1.D0*RR33*(Lam7*RR11*(RR13*RR31 + RR11*RR33) + Lam8*RR12*(RR13*RR32 + RR12*RR33) + RR13*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR31 + RR11*RR32) + Lam7*RR13*(RR13*RR31 + RR11&
  &*RR33) + RR11*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR13*(RR13*RR32 + RR12*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f3121 = -1.D0*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR21 + RR11*RR22) + Lam7*RR13*(RR13*RR21 + RR11*RR23) + RR11*(3.D0*Lam1*R&
  &R11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23))*RR31 - 1.D0*((Lam3 + Lam4 + Lam5)*RR11*(RR12*RR21 + RR11*RR22) + &
  &Lam8*RR13*(RR13*RR22 + RR12*RR23) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23))*RR32 - 1.D0*&
  &(Lam7*RR11*(RR13*RR21 + RR11*RR23) + Lam8*RR12*(RR13*RR22 + RR12*RR23) + RR13*(Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23))*RR33
CS1S1S1S1f3122 = -1.D0*RR33*(2.D0*Lam7*RR11*RR21*RR23 + 2.D0*Lam8*RR12*RR22*RR23 + RR13*(Lam7*DBLE(RR21**INT(2.D0)) + Lam8*DBLE(R&
  &R22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*RR21*RR22 + 2.D0*Lam7*RR13*RR2&
  &1*RR23 + RR11*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR22**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1&
  &.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR21*RR22 + 2.D0*Lam8*RR13*RR22*RR23 + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS1S1S1S1f3123 = -1.D0*RR33*(Lam7*RR11*(RR23*RR31 + RR21*RR33) + Lam8*RR12*(RR23*RR32 + RR22*RR33) + RR13*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR12*(RR22*RR31 + RR21*RR32) + Lam7*RR13*(RR23*RR31 + RR21&
  &*RR33) + RR11*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR13*(RR23*RR32 + RR22*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f3131 = -1.D0*RR33*(Lam7*RR11*(RR13*RR31 + RR11*RR33) + Lam8*RR12*(RR13*RR32 + RR12*RR33) + RR13*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR12*(RR12*RR31 + RR11*RR32) + Lam7*RR13*(RR13*RR31 + RR11&
  &*RR33) + RR11*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR13*(RR13*RR32 + RR12*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f3132 = -1.D0*RR33*(Lam7*RR11*(RR23*RR31 + RR21*RR33) + Lam8*RR12*(RR23*RR32 + RR22*RR33) + RR13*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR12*(RR22*RR31 + RR21*RR32) + Lam7*RR13*(RR23*RR31 + RR21&
  &*RR33) + RR11*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR11*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR13*(RR23*RR32 + RR22*RR33) + RR12*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f3133 = -1.D0*RR33*(2.D0*Lam7*RR11*RR31*RR33 + 2.D0*Lam8*RR12*RR32*RR33 + RR13*(Lam7*DBLE(RR31**INT(2.D0)) + Lam8*DBLE(R&
  &R32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR12*RR31*RR32 + 2.D0*Lam7*RR13*RR3&
  &1*RR33 + RR11*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR32**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1&
  &.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR31*RR32 + 2.D0*Lam8*RR13*RR32*RR33 + RR12*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS1S1S1S1f3211 = -1.D0*RR33*(2.D0*Lam7*RR11*RR13*RR21 + 2.D0*Lam8*RR12*RR13*RR22 + RR23*(Lam7*DBLE(RR11**INT(2.D0)) + Lam8*DBLE(R&
  &R12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR22 + 2.D0*Lam7*RR11*RR1&
  &3*RR23 + RR21*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR12**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1&
  &.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR21 + 2.D0*Lam8*RR12*RR13*RR23 + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS1S1S1S1f3212 = -1.D0*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR21 + RR11*RR22) + Lam7*RR23*(RR13*RR21 + RR11*RR23) + RR21*(3.D0*Lam1*R&
  &R11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23))*RR31 - 1.D0*((Lam3 + Lam4 + Lam5)*RR21*(RR12*RR21 + RR11*RR22) + &
  &Lam8*RR23*(RR13*RR22 + RR12*RR23) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23))*RR32 - 1.D0*&
  &(Lam7*RR21*(RR13*RR21 + RR11*RR23) + Lam8*RR22*(RR13*RR22 + RR12*RR23) + RR23*(Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23))*RR33
CS1S1S1S1f3213 = -1.D0*RR33*(Lam7*RR21*(RR13*RR31 + RR11*RR33) + Lam8*RR22*(RR13*RR32 + RR12*RR33) + RR23*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR31 + RR11*RR32) + Lam7*RR23*(RR13*RR31 + RR11&
  &*RR33) + RR21*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR23*(RR13*RR32 + RR12*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f3221 = -1.D0*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR21 + RR11*RR22) + Lam7*RR23*(RR13*RR21 + RR11*RR23) + RR21*(3.D0*Lam1*R&
  &R11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23))*RR31 - 1.D0*((Lam3 + Lam4 + Lam5)*RR21*(RR12*RR21 + RR11*RR22) + &
  &Lam8*RR23*(RR13*RR22 + RR12*RR23) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23))*RR32 - 1.D0*&
  &(Lam7*RR21*(RR13*RR21 + RR11*RR23) + Lam8*RR22*(RR13*RR22 + RR12*RR23) + RR23*(Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23))*RR33
CS1S1S1S1f3222 = -1.D0*RR33*(2.D0*Lam7*RR23*DBLE(RR21**INT(2.D0)) + 2.D0*Lam8*RR23*DBLE(RR22**INT(2.D0)) + RR23*(Lam7*DBLE(RR21**&
  &INT(2.D0)) + Lam8*DBLE(RR22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*DBLE(R&
  &R22**INT(2.D0)) + 2.D0*Lam7*RR21*DBLE(RR23**INT(2.D0)) + RR21*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR2&
  &2**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR22*DBLE(RR21**INT(2.D0)) + 2.D0*Lam8*RR2&
  &2*DBLE(RR23**INT(2.D0)) + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23*&
  &*INT(2.D0))))
CS1S1S1S1f3223 = -1.D0*RR33*(Lam7*RR21*(RR23*RR31 + RR21*RR33) + Lam8*RR22*(RR23*RR32 + RR22*RR33) + RR23*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR22*(RR22*RR31 + RR21*RR32) + Lam7*RR23*(RR23*RR31 + RR21&
  &*RR33) + RR21*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR23*(RR23*RR32 + RR22*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f3231 = -1.D0*RR33*(Lam7*RR21*(RR13*RR31 + RR11*RR33) + Lam8*RR22*(RR13*RR32 + RR12*RR33) + RR23*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR22*(RR12*RR31 + RR11*RR32) + Lam7*RR23*(RR13*RR31 + RR11&
  &*RR33) + RR21*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR23*(RR13*RR32 + RR12*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f3232 = -1.D0*RR33*(Lam7*RR21*(RR23*RR31 + RR21*RR33) + Lam8*RR22*(RR23*RR32 + RR22*RR33) + RR23*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR22*(RR22*RR31 + RR21*RR32) + Lam7*RR23*(RR23*RR31 + RR21&
  &*RR33) + RR21*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR21*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR23*(RR23*RR32 + RR22*RR33) + RR22*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f3233 = -1.D0*RR33*(2.D0*Lam7*RR21*RR31*RR33 + 2.D0*Lam8*RR22*RR32*RR33 + RR23*(Lam7*DBLE(RR31**INT(2.D0)) + Lam8*DBLE(R&
  &R32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR22*RR31*RR32 + 2.D0*Lam7*RR23*RR3&
  &1*RR33 + RR21*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR32**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1&
  &.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR31*RR32 + 2.D0*Lam8*RR23*RR32*RR33 + RR22*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS1S1S1S1f3311 = -1.D0*RR33*(2.D0*Lam7*RR11*RR13*RR31 + 2.D0*Lam8*RR12*RR13*RR32 + RR33*(Lam7*DBLE(RR11**INT(2.D0)) + Lam8*DBLE(R&
  &R12**INT(2.D0)) + 3.D0*Lam6*DBLE(RR13**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR32 + 2.D0*Lam7*RR11*RR1&
  &3*RR33 + RR31*(3.D0*Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR12**INT(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1&
  &.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR11*RR12*RR31 + 2.D0*Lam8*RR12*RR13*RR33 + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR11**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS1S1S1S1f3312 = -1.D0*RR31*((3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)*RR31 + (Lam3 + Lam4 + Lam5)*&
  &(RR12*RR21 + RR11*RR22)*RR32 + Lam7*(RR13*RR21 + RR11*RR23)*RR33) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*(RR12*RR21 + RR11*RR22)*RR&
  &31 + ((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23)*RR32 + Lam8*(RR13*RR22 + RR12*RR23)*RR33) - 1.D0*&
  &RR33*(Lam7*(RR13*RR21 + RR11*RR23)*RR31 + Lam8*(RR13*RR22 + RR12*RR23)*RR32 + (Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23)*RR33)
CS1S1S1S1f3313 = -1.D0*RR33*(Lam7*RR31*(RR13*RR31 + RR11*RR33) + Lam8*RR32*(RR13*RR32 + RR12*RR33) + RR33*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR32*(RR12*RR31 + RR11*RR32) + Lam7*RR33*(RR13*RR31 + RR11&
  &*RR33) + RR31*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR33*(RR13*RR32 + RR12*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f3321 = -1.D0*RR31*((3.D0*Lam1*RR11*RR21 + (Lam3 + Lam4 + Lam5)*RR12*RR22 + Lam7*RR13*RR23)*RR31 + (Lam3 + Lam4 + Lam5)*&
  &(RR12*RR21 + RR11*RR22)*RR32 + Lam7*(RR13*RR21 + RR11*RR23)*RR33) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*(RR12*RR21 + RR11*RR22)*RR&
  &31 + ((Lam3 + Lam4 + Lam5)*RR11*RR21 + 3.D0*Lam2*RR12*RR22 + Lam8*RR13*RR23)*RR32 + Lam8*(RR13*RR22 + RR12*RR23)*RR33) - 1.D0*&
  &RR33*(Lam7*(RR13*RR21 + RR11*RR23)*RR31 + Lam8*(RR13*RR22 + RR12*RR23)*RR32 + (Lam7*RR11*RR21 + Lam8*RR12*RR22 + 3.D0*Lam6*RR1&
  &3*RR23)*RR33)
CS1S1S1S1f3322 = -1.D0*RR33*(2.D0*Lam7*RR21*RR23*RR31 + 2.D0*Lam8*RR22*RR23*RR32 + RR33*(Lam7*DBLE(RR21**INT(2.D0)) + Lam8*DBLE(R&
  &R22**INT(2.D0)) + 3.D0*Lam6*DBLE(RR23**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR22*RR32 + 2.D0*Lam7*RR21*RR2&
  &3*RR33 + RR31*(3.D0*Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR22**INT(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1&
  &.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR21*RR22*RR31 + 2.D0*Lam8*RR22*RR23*RR33 + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR21**INT(2.D0&
  &)) + 3.D0*Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS1S1S1S1f3323 = -1.D0*RR33*(Lam7*RR31*(RR23*RR31 + RR21*RR33) + Lam8*RR32*(RR23*RR32 + RR22*RR33) + RR33*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR32*(RR22*RR31 + RR21*RR32) + Lam7*RR33*(RR23*RR31 + RR21&
  &*RR33) + RR31*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR33*(RR23*RR32 + RR22*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f3331 = -1.D0*RR33*(Lam7*RR31*(RR13*RR31 + RR11*RR33) + Lam8*RR32*(RR13*RR32 + RR12*RR33) + RR33*(Lam7*RR11*RR31 + Lam8*&
  &RR12*RR32 + 3.D0*Lam6*RR13*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR32*(RR12*RR31 + RR11*RR32) + Lam7*RR33*(RR13*RR31 + RR11&
  &*RR33) + RR31*(3.D0*Lam1*RR11*RR31 + (Lam3 + Lam4 + Lam5)*RR12*RR32 + Lam7*RR13*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR12*RR31 + RR11*RR32) + Lam8*RR33*(RR13*RR32 + RR12*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR11*RR31 + 3.D0*Lam2*RR12*RR32 + Lam&
  &8*RR13*RR33))
CS1S1S1S1f3332 = -1.D0*RR33*(Lam7*RR31*(RR23*RR31 + RR21*RR33) + Lam8*RR32*(RR23*RR32 + RR22*RR33) + RR33*(Lam7*RR21*RR31 + Lam8*&
  &RR22*RR32 + 3.D0*Lam6*RR23*RR33)) - 1.D0*RR31*((Lam3 + Lam4 + Lam5)*RR32*(RR22*RR31 + RR21*RR32) + Lam7*RR33*(RR23*RR31 + RR21&
  &*RR33) + RR31*(3.D0*Lam1*RR21*RR31 + (Lam3 + Lam4 + Lam5)*RR22*RR32 + Lam7*RR23*RR33)) - 1.D0*RR32*((Lam3 + Lam4 + Lam5)*RR31*&
  &(RR22*RR31 + RR21*RR32) + Lam8*RR33*(RR23*RR32 + RR22*RR33) + RR32*((Lam3 + Lam4 + Lam5)*RR21*RR31 + 3.D0*Lam2*RR22*RR32 + Lam&
  &8*RR23*RR33))
CS1S1S1S1f3333 = -1.D0*RR33*(2.D0*Lam7*RR33*DBLE(RR31**INT(2.D0)) + 2.D0*Lam8*RR33*DBLE(RR32**INT(2.D0)) + RR33*(Lam7*DBLE(RR31**&
  &INT(2.D0)) + Lam8*DBLE(RR32**INT(2.D0)) + 3.D0*Lam6*DBLE(RR33**INT(2.D0)))) - 1.D0*RR31*(2.D0*(Lam3 + Lam4 + Lam5)*RR31*DBLE(R&
  &R32**INT(2.D0)) + 2.D0*Lam7*RR31*DBLE(RR33**INT(2.D0)) + RR31*(3.D0*Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 + Lam5)*DBLE(RR3&
  &2**INT(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1.D0*RR32*(2.D0*(Lam3 + Lam4 + Lam5)*RR32*DBLE(RR31**INT(2.D0)) + 2.D0*Lam8*RR3&
  &2*DBLE(RR33**INT(2.D0)) + RR32*((Lam3 + Lam4 + Lam5)*DBLE(RR31**INT(2.D0)) + 3.D0*Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33*&
  &*INT(2.D0))))

CS2S2S2S2f1111 = -1.D0*SB*(2.D0*CB2*(Lam3 + Lam4 + Lam5)*SB + SB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam2*SB2)) - 1.D0*CB*(2.D0*CB*(&
  &Lam3 + Lam4 + Lam5)*SB2 + CB*(3.D0*CB2*Lam1 + (Lam3 + Lam4 + Lam5)*SB2))
CS2S2S2S2f1112 = -1.D0*SB*(SB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1&
  &.D0*CB*(CB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + (Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f1121 = -1.D0*SB*(SB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1&
  &.D0*CB*(CB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + (Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f1122 = -1.D0*CB*(-2.D0*CB*(Lam3 + Lam4 + Lam5)*SB2 + CB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam1*SB2)) - 1.D0*SB*(-2.D0*CB&
  &2*(Lam3 + Lam4 + Lam5)*SB + SB*(3.D0*CB2*Lam2 + (Lam3 + Lam4 + Lam5)*SB2))
CS2S2S2S2f1211 = -1.D0*SB*(-2.D0*CB*(Lam3 + Lam4 + Lam5)*SB2 + CB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam2*SB2)) - 1.D0*CB*(2.D0*CB2&
  &*(Lam3 + Lam4 + Lam5)*SB - 1.D0*SB*(3.D0*CB2*Lam1 + (Lam3 + Lam4 + Lam5)*SB2))
CS2S2S2S2f1212 = -1.D0*CB*(-1.D0*SB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) -&
  & 1.D0*SB*(CB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) - 1.D0*(Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f1221 = -1.D0*CB*(-1.D0*SB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) -&
  & 1.D0*SB*(CB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) - 1.D0*(Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f1222 = -1.D0*CB*(-2.D0*CB2*(Lam3 + Lam4 + Lam5)*SB - 1.D0*SB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam1*SB2)) - 1.D0*SB*(2.D&
  &0*CB*(Lam3 + Lam4 + Lam5)*SB2 + CB*(3.D0*CB2*Lam2 + (Lam3 + Lam4 + Lam5)*SB2))
CS2S2S2S2f2111 = -1.D0*CB*(2.D0*CB2*(Lam3 + Lam4 + Lam5)*SB + SB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam2*SB2)) + SB*(2.D0*CB*(Lam3 &
  &+ Lam4 + Lam5)*SB2 + CB*(3.D0*CB2*Lam1 + (Lam3 + Lam4 + Lam5)*SB2))
CS2S2S2S2f2112 = -1.D0*CB*(SB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) + S&
  &B*(CB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + (Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f2121 = -1.D0*CB*(SB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) + S&
  &B*(CB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + (Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f2122 = SB*(-2.D0*CB*(Lam3 + Lam4 + Lam5)*SB2 + CB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam1*SB2)) - 1.D0*CB*(-2.D0*CB2*(Lam&
  &3 + Lam4 + Lam5)*SB + SB*(3.D0*CB2*Lam2 + (Lam3 + Lam4 + Lam5)*SB2))
CS2S2S2S2f2211 = -1.D0*CB*(-2.D0*CB*(Lam3 + Lam4 + Lam5)*SB2 + CB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam2*SB2)) + SB*(2.D0*CB2*(Lam&
  &3 + Lam4 + Lam5)*SB - 1.D0*SB*(3.D0*CB2*Lam1 + (Lam3 + Lam4 + Lam5)*SB2))
CS2S2S2S2f2212 = SB*(-1.D0*SB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*&
  &CB*(CB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) - 1.D0*(Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f2221 = SB*(-1.D0*SB*(-3.D0*CB*Lam1*SB + CB*(Lam3 + Lam4 + Lam5)*SB) + CB*(Lam3 + Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*&
  &CB*(CB*(3.D0*CB*Lam2*SB - 1.D0*CB*(Lam3 + Lam4 + Lam5)*SB) - 1.D0*(Lam3 + Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2))
CS2S2S2S2f2222 = SB*(-2.D0*CB2*(Lam3 + Lam4 + Lam5)*SB - 1.D0*SB*(CB2*(Lam3 + Lam4 + Lam5) + 3.D0*Lam1*SB2)) - 1.D0*CB*(2.D0*CB*(&
  &Lam3 + Lam4 + Lam5)*SB2 + CB*(3.D0*CB2*Lam2 + (Lam3 + Lam4 + Lam5)*SB2))

CS3S3S3S3f1111 = -1.D0*SB*(2.D0*CB2*(Lam3 + Lam4)*SB + 2.D0*SB*(CB2*Lam5 + Lam2*SB2)) - 1.D0*CB*(2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D0&
  &*CB*(CB2*Lam1 + Lam5*SB2))
CS3S3S3S3f1112 = -1.D0*SB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*(-1.D0*&
  &CB*Lam1*SB + CB*Lam5*SB) + (Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f1121 = -1.D0*SB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*(-1.D0*&
  &CB*Lam1*SB + CB*Lam5*SB) + (Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f1122 = -1.D0*CB*(-2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D0*CB*(CB2*Lam5 + Lam1*SB2)) - 1.D0*SB*(-2.D0*CB2*(Lam3 + Lam4)*SB + 2.&
  &D0*SB*(CB2*Lam2 + Lam5*SB2))
CS3S3S3S3f1211 = -1.D0*SB*(-2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D0*CB*(CB2*Lam5 + Lam2*SB2)) - 1.D0*CB*(2.D0*CB2*(Lam3 + Lam4)*SB - 2.D&
  &0*SB*(CB2*Lam1 + Lam5*SB2))
CS3S3S3S3f1212 = -1.D0*CB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) - 1.D0*SB*(2.D0*CB*(CB*L&
  &am2*SB - 1.D0*CB*Lam5*SB) - 1.D0*(Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f1221 = -1.D0*CB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) - 1.D0*SB*(2.D0*CB*(CB*L&
  &am2*SB - 1.D0*CB*Lam5*SB) - 1.D0*(Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f1222 = -1.D0*CB*(-2.D0*CB2*(Lam3 + Lam4)*SB - 2.D0*SB*(CB2*Lam5 + Lam1*SB2)) - 1.D0*SB*(2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D&
  &0*CB*(CB2*Lam2 + Lam5*SB2))
CS3S3S3S3f2111 = -1.D0*CB*(2.D0*CB2*(Lam3 + Lam4)*SB + 2.D0*SB*(CB2*Lam5 + Lam2*SB2)) + SB*(2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D0*CB*(&
  &CB2*Lam1 + Lam5*SB2))
CS3S3S3S3f2112 = -1.D0*CB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) + SB*(2.D0*CB*(-1.D0*CB*La&
  &m1*SB + CB*Lam5*SB) + (Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f2121 = -1.D0*CB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) + SB*(2.D0*CB*(-1.D0*CB*La&
  &m1*SB + CB*Lam5*SB) + (Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f2122 = SB*(-2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D0*CB*(CB2*Lam5 + Lam1*SB2)) - 1.D0*CB*(-2.D0*CB2*(Lam3 + Lam4)*SB + 2.D0*SB*&
  &(CB2*Lam2 + Lam5*SB2))
CS3S3S3S3f2211 = -1.D0*CB*(-2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D0*CB*(CB2*Lam5 + Lam2*SB2)) + SB*(2.D0*CB2*(Lam3 + Lam4)*SB - 2.D0*SB*&
  &(CB2*Lam1 + Lam5*SB2))
CS3S3S3S3f2212 = SB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*(CB*Lam2*SB&
  & - 1.D0*CB*Lam5*SB) - 1.D0*(Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f2221 = SB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam5*SB) + CB*(Lam3 + Lam4)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*(CB*Lam2*SB&
  & - 1.D0*CB*Lam5*SB) - 1.D0*(Lam3 + Lam4)*SB*(CB2 - 1.D0*SB2))
CS3S3S3S3f2222 = SB*(-2.D0*CB2*(Lam3 + Lam4)*SB - 2.D0*SB*(CB2*Lam5 + Lam1*SB2)) - 1.D0*CB*(2.D0*CB*(Lam3 + Lam4)*SB2 + 2.D0*CB*(&
  &CB2*Lam2 + Lam5*SB2))

CS2S2S1S1f1111 = -1.D0*CB*(2.D0*Lam5*RR11*RR12*SB + CB*(Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR12**INT(2.D&
  &0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1.D0*SB*(2.D0*CB*Lam5*RR11*RR12 + SB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR11**INT(2.D0)) + L&
  &am2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS2S2S1S1f1112 = -1.D0*CB*(CB*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR13*RR23) + Lam5*(RR12*RR21 + RR11*RR&
  &22)*SB) - 1.D0*SB*(CB*Lam5*(RR12*RR21 + RR11*RR22) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23)*S&
  &B)
CS2S2S1S1f1113 = -1.D0*CB*(CB*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR13*RR33) + Lam5*(RR12*RR31 + RR11*RR&
  &32)*SB) - 1.D0*SB*(CB*Lam5*(RR12*RR31 + RR11*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33)*S&
  &B)
CS2S2S1S1f1121 = -1.D0*CB*(CB*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR13*RR23) + Lam5*(RR12*RR21 + RR11*RR&
  &22)*SB) - 1.D0*SB*(CB*Lam5*(RR12*RR21 + RR11*RR22) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23)*S&
  &B)
CS2S2S1S1f1122 = -1.D0*CB*(2.D0*Lam5*RR21*RR22*SB + CB*(Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR22**INT(2.D&
  &0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1.D0*SB*(2.D0*CB*Lam5*RR21*RR22 + SB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR21**INT(2.D0)) + L&
  &am2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS2S2S1S1f1123 = -1.D0*CB*(CB*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR23*RR33) + Lam5*(RR22*RR31 + RR21*RR&
  &32)*SB) - 1.D0*SB*(CB*Lam5*(RR22*RR31 + RR21*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33)*S&
  &B)
CS2S2S1S1f1131 = -1.D0*CB*(CB*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR13*RR33) + Lam5*(RR12*RR31 + RR11*RR&
  &32)*SB) - 1.D0*SB*(CB*Lam5*(RR12*RR31 + RR11*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33)*S&
  &B)
CS2S2S1S1f1132 = -1.D0*CB*(CB*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR23*RR33) + Lam5*(RR22*RR31 + RR21*RR&
  &32)*SB) - 1.D0*SB*(CB*Lam5*(RR22*RR31 + RR21*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33)*S&
  &B)
CS2S2S1S1f1133 = -1.D0*CB*(2.D0*Lam5*RR31*RR32*SB + CB*(Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR32**INT(2.D&
  &0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1.D0*SB*(2.D0*CB*Lam5*RR31*RR32 + SB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR31**INT(2.D0)) + L&
  &am2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS2S2S1S1f1211 = -1.D0*CB*(2.D0*CB*Lam5*RR11*RR12 - 1.D0*SB*(Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR12**IN&
  &T(2.D0)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1.D0*SB*(-2.D0*Lam5*RR11*RR12*SB + CB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR11**INT(2.D0&
  &)) + Lam2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS2S2S1S1f1212 = -1.D0*SB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23) - 1.D0*Lam5*(RR12*RR21 + RR&
  &11*RR22)*SB) - 1.D0*CB*(CB*Lam5*(RR12*RR21 + RR11*RR22) - 1.D0*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR&
  &13*RR23)*SB)
CS2S2S1S1f1213 = -1.D0*SB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33) - 1.D0*Lam5*(RR12*RR31 + RR&
  &11*RR32)*SB) - 1.D0*CB*(CB*Lam5*(RR12*RR31 + RR11*RR32) - 1.D0*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR&
  &13*RR33)*SB)
CS2S2S1S1f1221 = -1.D0*SB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23) - 1.D0*Lam5*(RR12*RR21 + RR&
  &11*RR22)*SB) - 1.D0*CB*(CB*Lam5*(RR12*RR21 + RR11*RR22) - 1.D0*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR&
  &13*RR23)*SB)
CS2S2S1S1f1222 = -1.D0*CB*(2.D0*CB*Lam5*RR21*RR22 - 1.D0*SB*(Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR22**IN&
  &T(2.D0)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1.D0*SB*(-2.D0*Lam5*RR21*RR22*SB + CB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR21**INT(2.D0&
  &)) + Lam2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS2S2S1S1f1223 = -1.D0*SB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33) - 1.D0*Lam5*(RR22*RR31 + RR&
  &21*RR32)*SB) - 1.D0*CB*(CB*Lam5*(RR22*RR31 + RR21*RR32) - 1.D0*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR&
  &23*RR33)*SB)
CS2S2S1S1f1231 = -1.D0*SB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33) - 1.D0*Lam5*(RR12*RR31 + RR&
  &11*RR32)*SB) - 1.D0*CB*(CB*Lam5*(RR12*RR31 + RR11*RR32) - 1.D0*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR&
  &13*RR33)*SB)
CS2S2S1S1f1232 = -1.D0*SB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33) - 1.D0*Lam5*(RR22*RR31 + RR&
  &21*RR32)*SB) - 1.D0*CB*(CB*Lam5*(RR22*RR31 + RR21*RR32) - 1.D0*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR&
  &23*RR33)*SB)
CS2S2S1S1f1233 = -1.D0*CB*(2.D0*CB*Lam5*RR31*RR32 - 1.D0*SB*(Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR32**IN&
  &T(2.D0)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1.D0*SB*(-2.D0*Lam5*RR31*RR32*SB + CB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR31**INT(2.D0&
  &)) + Lam2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS2S2S1S1f2111 = SB*(2.D0*Lam5*RR11*RR12*SB + CB*(Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR12**INT(2.D0)) + &
  &Lam7*DBLE(RR13**INT(2.D0)))) - 1.D0*CB*(2.D0*CB*Lam5*RR11*RR12 + SB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR11**INT(2.D0)) + Lam2*DB&
  &LE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS2S2S1S1f2112 = SB*(CB*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR13*RR23) + Lam5*(RR12*RR21 + RR11*RR22)*SB&
  &) - 1.D0*CB*(CB*Lam5*(RR12*RR21 + RR11*RR22) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23)*SB)
CS2S2S1S1f2113 = SB*(CB*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR13*RR33) + Lam5*(RR12*RR31 + RR11*RR32)*SB&
  &) - 1.D0*CB*(CB*Lam5*(RR12*RR31 + RR11*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33)*SB)
CS2S2S1S1f2121 = SB*(CB*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR13*RR23) + Lam5*(RR12*RR21 + RR11*RR22)*SB&
  &) - 1.D0*CB*(CB*Lam5*(RR12*RR21 + RR11*RR22) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23)*SB)
CS2S2S1S1f2122 = SB*(2.D0*Lam5*RR21*RR22*SB + CB*(Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR22**INT(2.D0)) + &
  &Lam7*DBLE(RR23**INT(2.D0)))) - 1.D0*CB*(2.D0*CB*Lam5*RR21*RR22 + SB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR21**INT(2.D0)) + Lam2*DB&
  &LE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS2S2S1S1f2123 = SB*(CB*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR23*RR33) + Lam5*(RR22*RR31 + RR21*RR32)*SB&
  &) - 1.D0*CB*(CB*Lam5*(RR22*RR31 + RR21*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33)*SB)
CS2S2S1S1f2131 = SB*(CB*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR13*RR33) + Lam5*(RR12*RR31 + RR11*RR32)*SB&
  &) - 1.D0*CB*(CB*Lam5*(RR12*RR31 + RR11*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33)*SB)
CS2S2S1S1f2132 = SB*(CB*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR23*RR33) + Lam5*(RR22*RR31 + RR21*RR32)*SB&
  &) - 1.D0*CB*(CB*Lam5*(RR22*RR31 + RR21*RR32) + ((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33)*SB)
CS2S2S1S1f2133 = SB*(2.D0*Lam5*RR31*RR32*SB + CB*(Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR32**INT(2.D0)) + &
  &Lam7*DBLE(RR33**INT(2.D0)))) - 1.D0*CB*(2.D0*CB*Lam5*RR31*RR32 + SB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR31**INT(2.D0)) + Lam2*DB&
  &LE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))
CS2S2S1S1f2211 = SB*(2.D0*CB*Lam5*RR11*RR12 - 1.D0*SB*(Lam1*DBLE(RR11**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR12**INT(2.D0&
  &)) + Lam7*DBLE(RR13**INT(2.D0)))) - 1.D0*CB*(-2.D0*Lam5*RR11*RR12*SB + CB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR11**INT(2.D0)) + L&
  &am2*DBLE(RR12**INT(2.D0)) + Lam8*DBLE(RR13**INT(2.D0))))
CS2S2S1S1f2212 = -1.D0*CB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23) - 1.D0*Lam5*(RR12*RR21 + RR&
  &11*RR22)*SB) + SB*(CB*Lam5*(RR12*RR21 + RR11*RR22) - 1.D0*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR13*RR&
  &23)*SB)
CS2S2S1S1f2213 = -1.D0*CB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33) - 1.D0*Lam5*(RR12*RR31 + RR&
  &11*RR32)*SB) + SB*(CB*Lam5*(RR12*RR31 + RR11*RR32) - 1.D0*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR13*RR&
  &33)*SB)
CS2S2S1S1f2221 = -1.D0*CB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR21 + Lam2*RR12*RR22 + Lam8*RR13*RR23) - 1.D0*Lam5*(RR12*RR21 + RR&
  &11*RR22)*SB) + SB*(CB*Lam5*(RR12*RR21 + RR11*RR22) - 1.D0*(Lam1*RR11*RR21 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR22 + Lam7*RR13*RR&
  &23)*SB)
CS2S2S1S1f2222 = SB*(2.D0*CB*Lam5*RR21*RR22 - 1.D0*SB*(Lam1*DBLE(RR21**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR22**INT(2.D0&
  &)) + Lam7*DBLE(RR23**INT(2.D0)))) - 1.D0*CB*(-2.D0*Lam5*RR21*RR22*SB + CB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR21**INT(2.D0)) + L&
  &am2*DBLE(RR22**INT(2.D0)) + Lam8*DBLE(RR23**INT(2.D0))))
CS2S2S1S1f2223 = -1.D0*CB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33) - 1.D0*Lam5*(RR22*RR31 + RR&
  &21*RR32)*SB) + SB*(CB*Lam5*(RR22*RR31 + RR21*RR32) - 1.D0*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR23*RR&
  &33)*SB)
CS2S2S1S1f2231 = -1.D0*CB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR11*RR31 + Lam2*RR12*RR32 + Lam8*RR13*RR33) - 1.D0*Lam5*(RR12*RR31 + RR&
  &11*RR32)*SB) + SB*(CB*Lam5*(RR12*RR31 + RR11*RR32) - 1.D0*(Lam1*RR11*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR12*RR32 + Lam7*RR13*RR&
  &33)*SB)
CS2S2S1S1f2232 = -1.D0*CB*(CB*((Lam3 + Lam4 - 1.D0*Lam5)*RR21*RR31 + Lam2*RR22*RR32 + Lam8*RR23*RR33) - 1.D0*Lam5*(RR22*RR31 + RR&
  &21*RR32)*SB) + SB*(CB*Lam5*(RR22*RR31 + RR21*RR32) - 1.D0*(Lam1*RR21*RR31 + (Lam3 + Lam4 - 1.D0*Lam5)*RR22*RR32 + Lam7*RR23*RR&
  &33)*SB)
CS2S2S1S1f2233 = SB*(2.D0*CB*Lam5*RR31*RR32 - 1.D0*SB*(Lam1*DBLE(RR31**INT(2.D0)) + (Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR32**INT(2.D0&
  &)) + Lam7*DBLE(RR33**INT(2.D0)))) - 1.D0*CB*(-2.D0*Lam5*RR31*RR32*SB + CB*((Lam3 + Lam4 - 1.D0*Lam5)*DBLE(RR31**INT(2.D0)) + L&
  &am2*DBLE(RR32**INT(2.D0)) + Lam8*DBLE(RR33**INT(2.D0))))

CS2S2S3S3f1111 = 0.5D0*(-1.D0*SB*(2.D0*CB2*(Lam4 + Lam5)*SB + 2.D0*SB*(CB2*Lam3 + Lam2*SB2)) - 1.D0*CB*(2.D0*CB*(Lam4 + Lam5)*SB2&
  & + 2.D0*CB*(CB2*Lam1 + Lam3*SB2)))
CS2S2S3S3f1112 = 0.5D0*(-1.D0*SB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*&
  &(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f1121 = 0.5D0*(-1.D0*SB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*&
  &(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f1122 = 0.5D0*(-1.D0*CB*(-2.D0*CB*(Lam4 + Lam5)*SB2 + 2.D0*CB*(CB2*Lam3 + Lam1*SB2)) - 1.D0*SB*(-2.D0*CB2*(Lam4 + Lam5)*&
  &SB + 2.D0*SB*(CB2*Lam2 + Lam3*SB2)))
CS2S2S3S3f1211 = 0.5D0*(-1.D0*SB*(-2.D0*CB*(Lam4 + Lam5)*SB2 + 2.D0*CB*(CB2*Lam3 + Lam2*SB2)) - 1.D0*CB*(2.D0*CB2*(Lam4 + Lam5)*S&
  &B - 2.D0*SB*(CB2*Lam1 + Lam3*SB2)))
CS2S2S3S3f1212 = 0.5D0*(-1.D0*CB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*SB*(2.D0*C&
  &B*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) - 1.D0*(Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f1221 = 0.5D0*(-1.D0*CB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*SB*(2.D0*C&
  &B*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) - 1.D0*(Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f1222 = 0.5D0*(-1.D0*CB*(-2.D0*CB2*(Lam4 + Lam5)*SB - 2.D0*SB*(CB2*Lam3 + Lam1*SB2)) - 1.D0*SB*(2.D0*CB*(Lam4 + Lam5)*SB&
  &2 + 2.D0*CB*(CB2*Lam2 + Lam3*SB2)))
CS2S2S3S3f2111 = 0.5D0*(-1.D0*CB*(2.D0*CB2*(Lam4 + Lam5)*SB + 2.D0*SB*(CB2*Lam3 + Lam2*SB2)) + SB*(2.D0*CB*(Lam4 + Lam5)*SB2 + 2.&
  &D0*CB*(CB2*Lam1 + Lam3*SB2)))
CS2S2S3S3f2112 = 0.5D0*(-1.D0*CB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) + SB*(2.D0*CB*(-1.D&
  &0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f2121 = 0.5D0*(-1.D0*CB*(2.D0*SB*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) + SB*(2.D0*CB*(-1.D&
  &0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f2122 = 0.5D0*(SB*(-2.D0*CB*(Lam4 + Lam5)*SB2 + 2.D0*CB*(CB2*Lam3 + Lam1*SB2)) - 1.D0*CB*(-2.D0*CB2*(Lam4 + Lam5)*SB + 2&
  &.D0*SB*(CB2*Lam2 + Lam3*SB2)))
CS2S2S3S3f2211 = 0.5D0*(-1.D0*CB*(-2.D0*CB*(Lam4 + Lam5)*SB2 + 2.D0*CB*(CB2*Lam3 + Lam2*SB2)) + SB*(2.D0*CB2*(Lam4 + Lam5)*SB - 2&
  &.D0*SB*(CB2*Lam1 + Lam3*SB2)))
CS2S2S3S3f2212 = 0.5D0*(SB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*(CB*&
  &Lam2*SB - 1.D0*CB*Lam3*SB) - 1.D0*(Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f2221 = 0.5D0*(SB*(-2.D0*SB*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + CB*(Lam4 + Lam5)*(CB2 - 1.D0*SB2)) - 1.D0*CB*(2.D0*CB*(CB*&
  &Lam2*SB - 1.D0*CB*Lam3*SB) - 1.D0*(Lam4 + Lam5)*SB*(CB2 - 1.D0*SB2)))
CS2S2S3S3f2222 = 0.5D0*(SB*(-2.D0*CB2*(Lam4 + Lam5)*SB - 2.D0*SB*(CB2*Lam3 + Lam1*SB2)) - 1.D0*CB*(2.D0*CB*(Lam4 + Lam5)*SB2 + 2.&
  &D0*CB*(CB2*Lam2 + Lam3*SB2)))

CS1S1S3S3f1111 = 0.5D0*(-1.D0*RR12*(2.D0*CB*(Lam4 + Lam5)*RR11*SB + 2.D0*RR12*(CB2*Lam3 + Lam2*SB2)) - 1.D0*RR11*(2.D0*CB*(Lam4 +&
  & Lam5)*RR12*SB + 2.D0*RR11*(CB2*Lam1 + Lam3*SB2)) - 2.D0*(CB2*Lam7 + Lam8*SB2)*DBLE(RR13**INT(2.D0)))
CS1S1S3S3f1112 = 0.5D0*(-1.D0*RR12*(2.D0*RR12*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (Lam4 + Lam5)*RR11*(CB2 - 1.D0*SB2)) - 1.D0*RR11*(&
  &2.D0*RR11*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR12*(CB2 - 1.D0*SB2)) - 2.D0*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*DBLE(R&
  &R13**INT(2.D0)))
CS1S1S3S3f1121 = 0.5D0*(-1.D0*RR12*(2.D0*RR12*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (Lam4 + Lam5)*RR11*(CB2 - 1.D0*SB2)) - 1.D0*RR11*(&
  &2.D0*RR11*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR12*(CB2 - 1.D0*SB2)) - 2.D0*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*DBLE(R&
  &R13**INT(2.D0)))
CS1S1S3S3f1122 = 0.5D0*(-1.D0*RR11*(-2.D0*CB*(Lam4 + Lam5)*RR12*SB + 2.D0*RR11*(CB2*Lam3 + Lam1*SB2)) - 1.D0*RR12*(-2.D0*CB*(Lam4&
  & + Lam5)*RR11*SB + 2.D0*RR12*(CB2*Lam2 + Lam3*SB2)) - 2.D0*(CB2*Lam8 + Lam7*SB2)*DBLE(RR13**INT(2.D0)))
CS1S1S3S3f1211 = 0.5D0*(-2.D0*RR13*RR23*(CB2*Lam7 + Lam8*SB2) - 1.D0*RR12*(2.D0*CB*(Lam4 + Lam5)*RR21*SB + 2.D0*RR22*(CB2*Lam3 + &
  &Lam2*SB2)) - 1.D0*RR11*(2.D0*CB*(Lam4 + Lam5)*RR22*SB + 2.D0*RR21*(CB2*Lam1 + Lam3*SB2)))
CS1S1S3S3f1212 = 0.5D0*(-2.D0*RR13*RR23*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR12*(2.D0*RR22*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR21*(CB2 - 1.D0*SB2)) - 1.D0*RR11*(2.D0*RR21*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR22*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f1221 = 0.5D0*(-2.D0*RR13*RR23*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR12*(2.D0*RR22*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR21*(CB2 - 1.D0*SB2)) - 1.D0*RR11*(2.D0*RR21*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR22*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f1222 = 0.5D0*(-2.D0*RR13*RR23*(CB2*Lam8 + Lam7*SB2) - 1.D0*RR11*(-2.D0*CB*(Lam4 + Lam5)*RR22*SB + 2.D0*RR21*(CB2*Lam3 +&
  & Lam1*SB2)) - 1.D0*RR12*(-2.D0*CB*(Lam4 + Lam5)*RR21*SB + 2.D0*RR22*(CB2*Lam2 + Lam3*SB2)))
CS1S1S3S3f1311 = 0.5D0*(-2.D0*RR13*RR33*(CB2*Lam7 + Lam8*SB2) - 1.D0*RR12*(2.D0*CB*(Lam4 + Lam5)*RR31*SB + 2.D0*RR32*(CB2*Lam3 + &
  &Lam2*SB2)) - 1.D0*RR11*(2.D0*CB*(Lam4 + Lam5)*RR32*SB + 2.D0*RR31*(CB2*Lam1 + Lam3*SB2)))
CS1S1S3S3f1312 = 0.5D0*(-2.D0*RR13*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR12*(2.D0*RR32*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR31*(CB2 - 1.D0*SB2)) - 1.D0*RR11*(2.D0*RR31*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR32*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f1321 = 0.5D0*(-2.D0*RR13*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR12*(2.D0*RR32*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR31*(CB2 - 1.D0*SB2)) - 1.D0*RR11*(2.D0*RR31*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR32*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f1322 = 0.5D0*(-2.D0*RR13*RR33*(CB2*Lam8 + Lam7*SB2) - 1.D0*RR11*(-2.D0*CB*(Lam4 + Lam5)*RR32*SB + 2.D0*RR31*(CB2*Lam3 +&
  & Lam1*SB2)) - 1.D0*RR12*(-2.D0*CB*(Lam4 + Lam5)*RR31*SB + 2.D0*RR32*(CB2*Lam2 + Lam3*SB2)))
CS1S1S3S3f2111 = 0.5D0*(-2.D0*RR13*RR23*(CB2*Lam7 + Lam8*SB2) - 1.D0*RR22*(2.D0*CB*(Lam4 + Lam5)*RR11*SB + 2.D0*RR12*(CB2*Lam3 + &
  &Lam2*SB2)) - 1.D0*RR21*(2.D0*CB*(Lam4 + Lam5)*RR12*SB + 2.D0*RR11*(CB2*Lam1 + Lam3*SB2)))
CS1S1S3S3f2112 = 0.5D0*(-2.D0*RR13*RR23*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR22*(2.D0*RR12*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR11*(CB2 - 1.D0*SB2)) - 1.D0*RR21*(2.D0*RR11*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR12*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f2121 = 0.5D0*(-2.D0*RR13*RR23*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR22*(2.D0*RR12*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR11*(CB2 - 1.D0*SB2)) - 1.D0*RR21*(2.D0*RR11*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR12*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f2122 = 0.5D0*(-2.D0*RR13*RR23*(CB2*Lam8 + Lam7*SB2) - 1.D0*RR21*(-2.D0*CB*(Lam4 + Lam5)*RR12*SB + 2.D0*RR11*(CB2*Lam3 +&
  & Lam1*SB2)) - 1.D0*RR22*(-2.D0*CB*(Lam4 + Lam5)*RR11*SB + 2.D0*RR12*(CB2*Lam2 + Lam3*SB2)))
CS1S1S3S3f2211 = 0.5D0*(-1.D0*RR22*(2.D0*CB*(Lam4 + Lam5)*RR21*SB + 2.D0*RR22*(CB2*Lam3 + Lam2*SB2)) - 1.D0*RR21*(2.D0*CB*(Lam4 +&
  & Lam5)*RR22*SB + 2.D0*RR21*(CB2*Lam1 + Lam3*SB2)) - 2.D0*(CB2*Lam7 + Lam8*SB2)*DBLE(RR23**INT(2.D0)))
CS1S1S3S3f2212 = 0.5D0*(-1.D0*RR22*(2.D0*RR22*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (Lam4 + Lam5)*RR21*(CB2 - 1.D0*SB2)) - 1.D0*RR21*(&
  &2.D0*RR21*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR22*(CB2 - 1.D0*SB2)) - 2.D0*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*DBLE(R&
  &R23**INT(2.D0)))
CS1S1S3S3f2221 = 0.5D0*(-1.D0*RR22*(2.D0*RR22*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (Lam4 + Lam5)*RR21*(CB2 - 1.D0*SB2)) - 1.D0*RR21*(&
  &2.D0*RR21*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR22*(CB2 - 1.D0*SB2)) - 2.D0*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*DBLE(R&
  &R23**INT(2.D0)))
CS1S1S3S3f2222 = 0.5D0*(-1.D0*RR21*(-2.D0*CB*(Lam4 + Lam5)*RR22*SB + 2.D0*RR21*(CB2*Lam3 + Lam1*SB2)) - 1.D0*RR22*(-2.D0*CB*(Lam4&
  & + Lam5)*RR21*SB + 2.D0*RR22*(CB2*Lam2 + Lam3*SB2)) - 2.D0*(CB2*Lam8 + Lam7*SB2)*DBLE(RR23**INT(2.D0)))
CS1S1S3S3f2311 = 0.5D0*(-2.D0*RR23*RR33*(CB2*Lam7 + Lam8*SB2) - 1.D0*RR22*(2.D0*CB*(Lam4 + Lam5)*RR31*SB + 2.D0*RR32*(CB2*Lam3 + &
  &Lam2*SB2)) - 1.D0*RR21*(2.D0*CB*(Lam4 + Lam5)*RR32*SB + 2.D0*RR31*(CB2*Lam1 + Lam3*SB2)))
CS1S1S3S3f2312 = 0.5D0*(-2.D0*RR23*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR22*(2.D0*RR32*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR31*(CB2 - 1.D0*SB2)) - 1.D0*RR21*(2.D0*RR31*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR32*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f2321 = 0.5D0*(-2.D0*RR23*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR22*(2.D0*RR32*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR31*(CB2 - 1.D0*SB2)) - 1.D0*RR21*(2.D0*RR31*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR32*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f2322 = 0.5D0*(-2.D0*RR23*RR33*(CB2*Lam8 + Lam7*SB2) - 1.D0*RR21*(-2.D0*CB*(Lam4 + Lam5)*RR32*SB + 2.D0*RR31*(CB2*Lam3 +&
  & Lam1*SB2)) - 1.D0*RR22*(-2.D0*CB*(Lam4 + Lam5)*RR31*SB + 2.D0*RR32*(CB2*Lam2 + Lam3*SB2)))
CS1S1S3S3f3111 = 0.5D0*(-2.D0*RR13*RR33*(CB2*Lam7 + Lam8*SB2) - 1.D0*RR32*(2.D0*CB*(Lam4 + Lam5)*RR11*SB + 2.D0*RR12*(CB2*Lam3 + &
  &Lam2*SB2)) - 1.D0*RR31*(2.D0*CB*(Lam4 + Lam5)*RR12*SB + 2.D0*RR11*(CB2*Lam1 + Lam3*SB2)))
CS1S1S3S3f3112 = 0.5D0*(-2.D0*RR13*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR32*(2.D0*RR12*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR11*(CB2 - 1.D0*SB2)) - 1.D0*RR31*(2.D0*RR11*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR12*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f3121 = 0.5D0*(-2.D0*RR13*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR32*(2.D0*RR12*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR11*(CB2 - 1.D0*SB2)) - 1.D0*RR31*(2.D0*RR11*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR12*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f3122 = 0.5D0*(-2.D0*RR13*RR33*(CB2*Lam8 + Lam7*SB2) - 1.D0*RR31*(-2.D0*CB*(Lam4 + Lam5)*RR12*SB + 2.D0*RR11*(CB2*Lam3 +&
  & Lam1*SB2)) - 1.D0*RR32*(-2.D0*CB*(Lam4 + Lam5)*RR11*SB + 2.D0*RR12*(CB2*Lam2 + Lam3*SB2)))
CS1S1S3S3f3211 = 0.5D0*(-2.D0*RR23*RR33*(CB2*Lam7 + Lam8*SB2) - 1.D0*RR32*(2.D0*CB*(Lam4 + Lam5)*RR21*SB + 2.D0*RR22*(CB2*Lam3 + &
  &Lam2*SB2)) - 1.D0*RR31*(2.D0*CB*(Lam4 + Lam5)*RR22*SB + 2.D0*RR21*(CB2*Lam1 + Lam3*SB2)))
CS1S1S3S3f3212 = 0.5D0*(-2.D0*RR23*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR32*(2.D0*RR22*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR21*(CB2 - 1.D0*SB2)) - 1.D0*RR31*(2.D0*RR21*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR22*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f3221 = 0.5D0*(-2.D0*RR23*RR33*(-1.D0*CB*Lam7*SB + CB*Lam8*SB) - 1.D0*RR32*(2.D0*RR22*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (&
  &Lam4 + Lam5)*RR21*(CB2 - 1.D0*SB2)) - 1.D0*RR31*(2.D0*RR21*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR22*(CB2 - 1.D0*SB&
  &2)))
CS1S1S3S3f3222 = 0.5D0*(-2.D0*RR23*RR33*(CB2*Lam8 + Lam7*SB2) - 1.D0*RR31*(-2.D0*CB*(Lam4 + Lam5)*RR22*SB + 2.D0*RR21*(CB2*Lam3 +&
  & Lam1*SB2)) - 1.D0*RR32*(-2.D0*CB*(Lam4 + Lam5)*RR21*SB + 2.D0*RR22*(CB2*Lam2 + Lam3*SB2)))
CS1S1S3S3f3311 = 0.5D0*(-1.D0*RR32*(2.D0*CB*(Lam4 + Lam5)*RR31*SB + 2.D0*RR32*(CB2*Lam3 + Lam2*SB2)) - 1.D0*RR31*(2.D0*CB*(Lam4 +&
  & Lam5)*RR32*SB + 2.D0*RR31*(CB2*Lam1 + Lam3*SB2)) - 2.D0*(CB2*Lam7 + Lam8*SB2)*DBLE(RR33**INT(2.D0)))
CS1S1S3S3f3312 = 0.5D0*(-1.D0*RR32*(2.D0*RR32*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (Lam4 + Lam5)*RR31*(CB2 - 1.D0*SB2)) - 1.D0*RR31*(&
  &2.D0*RR31*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR32*(CB2 - 1.D0*SB2)) - 2.D0*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*DBLE(R&
  &R33**INT(2.D0)))
CS1S1S3S3f3321 = 0.5D0*(-1.D0*RR32*(2.D0*RR32*(CB*Lam2*SB - 1.D0*CB*Lam3*SB) + (Lam4 + Lam5)*RR31*(CB2 - 1.D0*SB2)) - 1.D0*RR31*(&
  &2.D0*RR31*(-1.D0*CB*Lam1*SB + CB*Lam3*SB) + (Lam4 + Lam5)*RR32*(CB2 - 1.D0*SB2)) - 2.D0*(-1.D0*CB*Lam7*SB + CB*Lam8*SB)*DBLE(R&
  &R33**INT(2.D0)))
CS1S1S3S3f3322 = 0.5D0*(-1.D0*RR31*(-2.D0*CB*(Lam4 + Lam5)*RR32*SB + 2.D0*RR31*(CB2*Lam3 + Lam1*SB2)) - 1.D0*RR32*(-2.D0*CB*(Lam4&
  & + Lam5)*RR31*SB + 2.D0*RR32*(CB2*Lam2 + Lam3*SB2)) - 2.D0*(CB2*Lam8 + Lam7*SB2)*DBLE(RR33**INT(2.D0)))

                ! Kinematic prefactor together with the symmetry factor of the process
                prefactor = 1D0/1D0 * DSQRT(MA0**4 + MD**4 + MD**4 - 2D0*MA0**2*MD**2 - 2D0*MA0**2*MD**2 &
                            &  - 2D0*MD**2*MD**2 )/(16D0*PI*MA0**3)

                ! Get the full tree-level decay width
                call clearcache
                treeLevelWidth = prefactor*A0toDDBarTree()

                ! Calculate the NLO ingredients
                treeLevelTemp = A0toDDBarTree()
                vertexCorrectionsTemp = A0toDDBarVC()
                vertexTadpolesTemp = A0toDDBarTad()
                realCorrectionsTemp = A0toDDBarReal()

                ! Calculate the NLO width w/o counterterm contributions (debug mode only)
                if (debugModeOn) then
                    call clearcache
                    NLOVCwidth = prefactor*( 2D0*DBLE(vertexCorrectionsTemp) + realCorrectionsTemp )
                    call clearcache
                    NLOVCwoIRwidth = prefactor*( 2D0*DBLE(vertexCorrectionsTemp) )
                    NLOIRonlywidth = prefactor*( realCorrectionsTemp )
                end if

                ! Get the full NLO decay width for all schemes
                do m = 1, maxNumberSchemes, 1
                    call clearcache

                    ! Schemes 1, 2, and 9 are without tadpoles
                    if ((m == 1) .OR. (m == 2) .OR. (m == 9)) then
                        fullamplitude(m) = treeLevelTemp + 2D0*DBLE(vertexCorrectionsTemp) + &
                                            & 2D0*A0toDDBarCT(m)
                        call clearcache
                        NLOWidth(m) = prefactor*( fullamplitude(m) + realCorrectionsTemp )
                    else
                        fullamplitude(m) = treeLevelTemp + 2D0*DBLE(vertexCorrectionsTemp + vertexTadpolesTemp) + &
                                            & 2D0*A0toDDBarCT(m)
                        call clearcache
                        NLOWidth(m) = prefactor*( fullamplitude(m) + realCorrectionsTemp )
                    end if
                end do

                ! Write the results to the output string
                ! Format: we use 18 characters in total for double precision. 3 are reserved for the exponent, 1 for the sign of the exponent,
                ! 1 for the symbol E denoting the exponent, 1 for the dot, 1 for a possible negative sign, 1 for the digit before the comma
                ! and 10 for the digits after the comma
                write( tempVal, '(ES18.10E3)' ) MH1
                outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                write( tempVal, '(ES18.10E3)' ) MH2
                outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                write( tempVal, '(ES18.10E3)' ) MH3
                outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                write( tempVal, '(ES18.10E3)' ) MA0
                outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                write( tempVal, '(ES18.10E3)' ) MHp
                outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                write( tempVal, '(ES18.10E3)' ) alpha1
                outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                write( tempVal, '(ES18.10E3)' ) alpha2
                outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                write( tempVal, '(ES18.10E3)' ) alpha3
                outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                write( tempVal, '(ES18.10E3)' ) beta
                outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                write( tempVal, '(ES18.10E3)' ) m12squared
                outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                write( tempVal, '(I1)' ) TypeOf2HDM
                outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                write( tempVal, '(ES18.10E3)' ) InputScale
                outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                write( tempVal, '(ES18.10E3)' ) treeLevelWidth
                outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                if (debugModeOn) then
                    ! write( tempVal, '(ES18.10E3)' ) NLOVCwidth
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) NLOVCwoIRwidth
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) NLOIRonlywidth
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dMW2Usual()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dMW2Alter()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dMZ2Usual()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dMZ2Alter()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dMLOSUsual()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dMLOSAlter()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dMBOSUsual()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dMBOSAlter()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dZHHHHOS()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dZHHh0OSUsual()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dZHHh0OSAlter()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dZh0HHOSUsual()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dZh0HHOSAlter()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dZh0h0OS()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dZA0A0OS()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dZG0A0OSUsual()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dZG0A0OSAlter()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dZTauTauOSLeft()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dZTauTauOSRight()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dZBBOSLeft()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dZBBOSRight()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dBeta1KanUsual()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dBeta1KanAlter()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dBeta1PinchPStar()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    ! write( tempVal, '(ES18.10E3)' ) dBeta1PinchOS()
                    ! outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                end if
                do m = 1, maxNumberSchemes, 1
                    write( tempVal, '(ES18.10E3)') NLOWidth(m)
                    outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                end do
                write( tempVal, '(ES18.10E3)' ) ((treeLevelWidth - treeLevelWidth)*100D0/treeLevelWidth)
                outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                do m = 1, maxNumberSchemes, 1
                    write( tempVal, '(ES18.10E3)') ((NLOWidth(m) - treeLevelWidth)*100D0/treeLevelWidth)
                    if (m == maxNumberSchemes) then
                        outputFileContent = trim(outputFileContent) // (trim(tempVal) // "\n")
                    else
                        outputFileContent = trim(outputFileContent) // (trim(tempVal) // ",")
                    end if
                end do
            end do

            ! Write the results to the output file
            open(unit=44, file=trim(pathToOutputFiles)//trim(targetName), status='new', &
            &action='write', iostat=statWrite)
                if ( statWrite == 0) then
                    write(44,*) trim(outputFileContent)
                else
                   write(*,*) 'ERROR: could not create output file for writing!'
                end if
           close(unit=44)
        call ltexi
    end if

end program decayWidth
