# N2HDMCalc

A program for the Semi-Automated Calculation of Electroweak One-Loop Corrections to Higgs Decays in the Next-to-Two-Higgs-Doublet Model

## Program information

**Program** N2HDMCalc 1.0.2

**Date** 04.11.2019

**Author** [Marcel Krause](mailto:marcel.krause@alumni.kit.edu)

**Changelog** For a documentation about the changes made in N2HDMCalc, check the [Changelog.md](Changelog.md) file.

**Copyright** Copyright (C) 2019, Marcel Krause

**License** GNU General Public License (GNU GPL-3.0-or-later). N2HDMCalc is released under GNU General Public License (GNU GPL-3.0-or-later). This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You have received a copy ([LICENSE.md](LICENSE.md)) of the GNU General Public License along with this program.

**Contact** For feedback, complaints and bug reports, please send an e-mail to <marcel.krause@alumni.kit.edu>.

### System requirements

Supported operating systems:
- Windows 7 and Windows 10 (with [Cygwin](https://www.cygwin.com/ "Cygwin") installed)
- Linux
- macOS

The following components must be installed on your system to compile and run N2HDMCalc:
- GNU gcc (tested with versions 6.4.0 and 7.3.1)
- GNU g++
- GNU gfortran
- Python 2 or 3 (tested with versions 2.7.14 and 3.5.0)

If you install N2HDMCalc on Windows, make sure that you install [Cygwin](https://www.cygwin.com/ "Cygwin") first, including the following components:
- GNU gcc
- GNU g++
- GNU gfortran

### Installing

For an easy installation, we recommend using the automatic installer of N2HDMCalc, which guides you through the installation. Download the latest version of N2HDMCalc from https://github.com/marcel-krause/N2HDMCalc. Open a shell, navigate to the N2HDMCalc root folder and execute the following command:
```
python N2HDMCalc.py
```
The main N2HDMCalc script asks you whether you want to start the installation routine and subsequently guides you through the installation.

### Using N2HDMCalc

To start N2HDMCalc, execute the following command in the main folder of N2HDMCalc:
```
python N2HDMCalc.py
```
The program ask you what process you want to calculate. For an overview over the implemented particles, you can use the
```
help
```
command. 