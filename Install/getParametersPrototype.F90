subroutine getParameters(filePath)
    use constants
    implicit none
    character(6) dump
    double precision dump2
    double precision MH2Temp, MH3Temp, MA0Temp, MHpTemp, alpha1Temp, alpha2Temp, alpha3Temp, TBTemp, m12squaredTemp, vSTemp
    integer TypeOf2HDMTemp    
    integer statOpen, statRead
    integer :: currentLine = 1
    logical fileExistsSM, fileExists2HDM
    character isContinue
    character(50), intent(in) :: filePath
    character(300), parameter :: pathToInputFiles = 'PLACEHOLDERRESULTPATH'

    ! Check if the SM input file exists
    inquire(file=trim(pathToInputFiles)//'ParametersSM.txt', exist=fileExistsSM)
    if (.not. fileExistsSM) then
        do
            print *, "ERROR: Could not find ParametersSM.txt!"
            print *, ">>> Do you want to continue with the evaluation of the program? [y/n]"
            read (*,*) isContinue
            if (isContinue == 'n') then
                print *, "Termination requested by user. 2HDMCalc will be terminated now."
                stop
            else if (isContinue == 'y') then
                exit
            else
                print *, "Invalid character. Enter y or n."
            end if
        end do
    end if

    ! Read the SM input parameters
    open(unit=42, file=trim(pathToInputFiles)//'ParametersSM.txt', iostat=statOpen)
        if (statOpen == 0) then
            read(42, *) dump, MH1
            read(42, *) dump, MW
            read(42, *) dump, MZ
            read(42, *) dump, ME
            read(42, *) dump, MM
            read(42, *) dump, ML
            read(42, *) dump, MU
            read(42, *) dump, MC
            read(42, *) dump, MT
            read(42, *) dump, MD
            read(42, *) dump, MS
            read(42, *) dump, MB
            read(42, *) dump, EL
            read(42, *) dump, SW
            read(42, *) dump, CW
            read(42, *) dump, CKM11
            read(42, *) dump, CKM12
            read(42, *) dump, CKM13
            read(42, *) dump, CKM21
            read(42, *) dump, CKM22
            read(42, *) dump, CKM23
            read(42, *) dump, CKM31
            read(42, *) dump, CKM32
            read(42, *) dump, CKM33
            read(42, *) dump, CKMC11
            read(42, *) dump, CKMC12
            read(42, *) dump, CKMC13
            read(42, *) dump, CKMC21
            read(42, *) dump, CKMC22
            read(42, *) dump, CKMC23
            read(42, *) dump, CKMC31
            read(42, *) dump, CKMC32
            read(42, *) dump, CKMC33
            read(42, *) dump, InputScale
        else
            do
                print *, "ERROR: Generic error when reading file ParametersSM.txt!"
                print *, ">>> Do you want to continue with the evaluation of the program? [y/n]"
                read (*,*) isContinue
                if (isContinue == 'n') then
                    print *, "Termination requested by user. 2HDMCalc will be terminated now."
                    stop
                else if (isContinue == 'y') then
                    exit
                else
                    print *, "Invalid character. Enter y or n."
                end if
            end do
        end if
    close(42)

    ! Generate the square of the SM input parameters
    MH12 = MH1**2
    MW2 = MW**2
    MZ2 = MZ**2
    ME2 = ME**2
    MM2 = MM**2
    ML2 = ML**2
    MU2 = MU**2
    MC2 = MC**2
    MT2 = MT**2
    MD2 = MD**2
    MS2 = MS**2
    MB2 = MB**2
    EL2 = EL**2
    SW2 = SW**2
    CW2 = CW**2

    ! Check if the 2HDM input file exists
    inquire(file=trim(pathToInputFiles)//trim(filePath), exist=fileExists2HDM)
    if (.not. fileExists2HDM) then
        do
            print *, "ERROR: Could not find the 2HDM input parameter file!"
            print *, ">>> Do you want to continue with the evaluation of the program? [y/n]"
            read (*,*) isContinue
            if (isContinue == 'n') then
                print *, "Termination requested by user. 2HDMCalc will be terminated now."
                stop
            else if (isContinue == 'y') then
                exit
            else
                print *, "Invalid character. Enter y or n."
            end if
        end do
    end if

    ! Read the 2HDM input parameters
    open(unit=43, file=trim(pathToInputFiles)//trim(filePath), iostat=statOpen)
        if (statOpen == 0) then
            do
                ! Read the parameters from the current line
                ! read(43, *, iostat=statRead) MHHTemp, dump2, MA0Temp, MHpTemp, alphaTemp, TBTemp, dump2, dump2, dump2, &
                !     & dump2, dump2, dump2, dump2, m12squaredTemp
                read(43, *, iostat=statRead) MH2Temp, MH3Temp, MA0Temp, MHpTemp, alpha1Temp, alpha2Temp, alpha3Temp,&
                & TBTemp, m12squaredTemp, vSTemp, TypeOf2HDMTemp

                ! Stop reading the file at end-of-file
                if (statRead /= 0) exit
                !if (currentLine == 1000) exit

                !print *, "Reading line ", currentLine, "..."

                ! If end-of-file is not reached yet, copy the values into the lists
                MH2List(currentLine) = MH2Temp
                MH3List(currentLine) = MH3Temp
                MA0List(currentLine) = MA0Temp
                MHpList(currentLine) = MHpTemp
                alpha1List(currentLine) = alpha1Temp
                alpha2List(currentLine) = alpha2Temp
                alpha3List(currentLine) = alpha3Temp
                TBList(currentLine) = TBTemp
                m12squaredList(currentLine) = m12squaredTemp
                vSList(currentLine) = vSTemp
                TypeOf2HDMList(currentLine) = TypeOf2HDMTemp

                ! Calculate the current beta with the given tan(beta)
                betaList(currentLine) = datan(TBList(currentLine))

                ! Calculate the rest of the 2HDM parameters
                CA1List(currentLine) = dcos(alpha1List(currentLine))
                CA2List(currentLine) = dcos(alpha2List(currentLine))
                CA3List(currentLine) = dcos(alpha3List(currentLine))
                SA1List(currentLine) = dsin(alpha1List(currentLine))
                SA2List(currentLine) = dsin(alpha2List(currentLine))
                SA3List(currentLine) = dsin(alpha3List(currentLine))
                CBList(currentLine) = dcos(betaList(currentLine))
                SBList(currentLine) = dsin(betaList(currentLine))
                if (TypeOf2HDMList(currentLine) == 1) then
                    YukS1Lep1List(currentLine) = CA2List(currentLine)*SA1List(currentLine)/SBList(currentLine)
                    YukS1Lep2List(currentLine) = ( CA1List(currentLine)*CA3List(currentLine) - SA1List(currentLine)*&
                                                    &SA2List(currentLine)*SA3List(currentLine) )/SBList(currentLine)
                    YukS1Lep3List(currentLine) = -( CA1List(currentLine)*SA3List(currentLine) + CA3List(currentLine)*&
                                                    &SA1List(currentLine)*SA2List(currentLine) )/SBList(currentLine)
                    YukS2Lep1List(currentLine) = 1D0
                    YukS2Lep2List(currentLine) = CBList(currentLine)/SBList(currentLine)
                    YukS3Lep1List(currentLine) = 1D0
                    YukS3Lep2List(currentLine) = CBList(currentLine)/SBList(currentLine)
                    YukS1Quark1List(currentLine) = CA2List(currentLine)*SA1List(currentLine)/SBList(currentLine)
                    YukS1Quark2List(currentLine) = ( CA1List(currentLine)*CA3List(currentLine) - SA1List(currentLine)*&
                                                    &SA2List(currentLine)*SA3List(currentLine) )/SBList(currentLine)
                    YukS1Quark3List(currentLine) = -( CA1List(currentLine)*SA3List(currentLine) + CA3List(currentLine)*&
                                                    &SA1List(currentLine)*SA2List(currentLine) )/SBList(currentLine)
                    YukS2Quark1List(currentLine) = 1D0
                    YukS2Quark2List(currentLine) = CBList(currentLine)/SBList(currentLine)
                    YukS3Quark1List(currentLine) = 1D0
                    YukS3Quark2List(currentLine) = CBList(currentLine)/SBList(currentLine)
                else if (TypeOf2HDMList(currentLine) == 2) then
                    YukS1Lep1List(currentLine) = CA1List(currentLine)*CA2List(currentLine)/CBList(currentLine)
                    YukS1Lep2List(currentLine) = -( CA3List(currentLine)*SA1List(currentLine) + CA1List(currentLine)*&
                                                    &SA2List(currentLine)*SA3List(currentLine) )/CBList(currentLine)
                    YukS1Lep3List(currentLine) = ( SA1List(currentLine)*SA3List(currentLine) - CA1List(currentLine)*&
                                                    &CA3List(currentLine)*SA2List(currentLine) )/CBList(currentLine)
                    YukS2Lep1List(currentLine) = 1D0
                    YukS2Lep2List(currentLine) = -SBList(currentLine)/CBList(currentLine)
                    YukS3Lep1List(currentLine) = 1D0
                    YukS3Lep2List(currentLine) = -SBList(currentLine)/CBList(currentLine)
                    YukS1Quark1List(currentLine) = CA1List(currentLine)*CA2List(currentLine)/CBList(currentLine)
                    YukS1Quark2List(currentLine) = -( CA3List(currentLine)*SA1List(currentLine) + CA1List(currentLine)*&
                                                    &SA2List(currentLine)*SA3List(currentLine) )/CBList(currentLine)
                    YukS1Quark3List(currentLine) = ( SA1List(currentLine)*SA3List(currentLine) - CA1List(currentLine)*&
                                                    &CA3List(currentLine)*SA2List(currentLine) )/CBList(currentLine)
                    YukS2Quark1List(currentLine) = 1D0
                    YukS2Quark2List(currentLine) = -SBList(currentLine)/CBList(currentLine)
                    YukS3Quark1List(currentLine) = 1D0
                    YukS3Quark2List(currentLine) = -SBList(currentLine)/CBList(currentLine)
                else if (TypeOf2HDMList(currentLine) == 3) then
                    YukS1Lep1List(currentLine) = CA1List(currentLine)*CA2List(currentLine)/CBList(currentLine)
                    YukS1Lep2List(currentLine) = -( CA3List(currentLine)*SA1List(currentLine) + CA1List(currentLine)*&
                                                    &SA2List(currentLine)*SA3List(currentLine) )/CBList(currentLine)
                    YukS1Lep3List(currentLine) = ( SA1List(currentLine)*SA3List(currentLine) - CA1List(currentLine)*&
                                                    &CA3List(currentLine)*SA2List(currentLine) )/CBList(currentLine)
                    YukS2Lep1List(currentLine) = 1D0
                    YukS2Lep2List(currentLine) = -SBList(currentLine)/CBList(currentLine)
                    YukS3Lep1List(currentLine) = 1D0
                    YukS3Lep2List(currentLine) = -SBList(currentLine)/CBList(currentLine)
                    YukS1Quark1List(currentLine) = CA2List(currentLine)*SA1List(currentLine)/SBList(currentLine)
                    YukS1Quark2List(currentLine) = ( CA1List(currentLine)*CA3List(currentLine) - SA1List(currentLine)*&
                                                    &SA2List(currentLine)*SA3List(currentLine) )/SBList(currentLine)
                    YukS1Quark3List(currentLine) = -( CA1List(currentLine)*SA3List(currentLine) + CA3List(currentLine)*&
                                                    &SA1List(currentLine)*SA2List(currentLine) )/SBList(currentLine)
                    YukS2Quark1List(currentLine) = 1D0
                    YukS2Quark2List(currentLine) = CBList(currentLine)/SBList(currentLine)
                    YukS3Quark1List(currentLine) = 1D0
                    YukS3Quark2List(currentLine) = CBList(currentLine)/SBList(currentLine)
                else
                    YukS1Lep1List(currentLine) = CA2List(currentLine)*SA1List(currentLine)/SBList(currentLine)
                    YukS1Lep2List(currentLine) = ( CA1List(currentLine)*CA3List(currentLine) - SA1List(currentLine)*&
                                                    &SA2List(currentLine)*SA3List(currentLine) )/SBList(currentLine)
                    YukS1Lep3List(currentLine) = -( CA1List(currentLine)*SA3List(currentLine) + CA3List(currentLine)*&
                                                    &SA1List(currentLine)*SA2List(currentLine) )/SBList(currentLine)
                    YukS2Lep1List(currentLine) = 1D0
                    YukS2Lep2List(currentLine) = CBList(currentLine)/SBList(currentLine)
                    YukS3Lep1List(currentLine) = 1D0
                    YukS3Lep2List(currentLine) = CBList(currentLine)/SBList(currentLine)
                    YukS1Quark1List(currentLine) = CA1List(currentLine)*CA2List(currentLine)/CBList(currentLine)
                    YukS1Quark2List(currentLine) = -( CA3List(currentLine)*SA1List(currentLine) + CA1List(currentLine)*&
                                                    &SA2List(currentLine)*SA3List(currentLine) )/CBList(currentLine)
                    YukS1Quark3List(currentLine) = ( SA1List(currentLine)*SA3List(currentLine) - CA1List(currentLine)*&
                                                    &CA3List(currentLine)*SA2List(currentLine) )/CBList(currentLine)
                    YukS2Quark1List(currentLine) = 1D0
                    YukS2Quark2List(currentLine) = -SBList(currentLine)/CBList(currentLine)
                    YukS3Quark1List(currentLine) = 1D0
                    YukS3Quark2List(currentLine) = -SBList(currentLine)/CBList(currentLine)
                end if

                ! Set the maximum number of points contained in the file to the current line number
                maxPoint = currentLine

                ! Increment the line counter
                currentLine = currentLine + 1
            end do
            numberOfPoints = currentLine - 1
        else
            do
                print *, "ERROR: Generic error when reading the 2HDM input parameter file!"
                print *, ">>> Do you want to continue with the evaluation of the program? [y/n]"
                read (*,*) isContinue
                if (isContinue == 'n') then
                    print *, "Termination requested by user. 2HDMCalc will be terminated now."
                    stop
                else if (isContinue == 'y') then
                    exit
                else
                    print *, "Invalid character. Enter y or n."
                end if
            end do
        end if
    close(43)

end subroutine getParameters
