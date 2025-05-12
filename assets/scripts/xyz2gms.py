#!/usr/bin/env python3
"""
This script prepares a set of
input files for each xyz structure provided as an
input. An input template is required. If none is found a standard
one is prepared.

If -g (--grep) option is used the data files should be the log files from which
the energies are extracted.

Usage: xyz2gzm.py [options] *.xyz or *.log file(s)

Options:
  -h, --help       show this help
  -g, --grep       extract energies from logs
  -e, --eds        prepare files for EDS analysis
  -m, --mrsf       prepare files for MRSF calculations
"""

#     Copyright (C) 2011, Robert W. Gora (robert.gora@pwr.wroc.pl)
#  
#     This program is free software; you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation; either version 2 of the License, or (at your
#     option) any later version.
#  
#     This program is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#     Public License for more details.
#  
#     You should have received a copy of the GNU General Public License along
#     with this program; if not, write to the Free Software Foundation, Inc.,
#     675 Mass Ave, Cambridge, MA 02139, USA.

__author__  = "Robert Gora (robert.gora@pwr.edu.pl)"
__version__ = list(filter(str.isdigit, "$Revision: 1 $"))

# Import necessary modules
import os, sys, getopt, re

from string import Template
from numpy import *

# Regular expressions
reflags = re.DOTALL

#----------------------------------------------------------------------------
# Usage
#----------------------------------------------------------------------------
def Usage():
    """Print usage information and exit."""
    print(__doc__)

    sys.exit()

#----------------------------------------------------------------------------
# Main
#----------------------------------------------------------------------------
def Main(argv):
    '''Parse commandline and loop throught the logs'''

    data = {}

    # Set up defaults
    eds = False
    mrsf = False
    grep = False

    # Parse commandline
    try:
        opts, args = getopt.getopt(argv, "hgem",
                                        ["help",
                                         "grep",
                                         "eds",
                                         "mrsf"
                                         ])
    except getopt.GetoptError as error:
        print(error)
        Usage()
    if not argv:
        Usage()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            Usage()
        elif opt in ("-g", "--grep"):
            grep=True
        elif opt in ("-e", "--eds"):
            eds=True
        elif opt in ("-m", "--mrsf"):
            mrsf=True

    if not (eds or grep or mrsf):
        print("Choose at least one type of job")
        sys.exit()

    # Parse each data file (with xyz coords) or dir (with the results)
    data_files = args

    if grep:
        # Exctract energies
        if eds:
            for f in data_files:
                l=EDS_LOGS(f)
            EDS_LOGS.PrintEnergies(l)
    else:
        # Prepare inputs
        for f in data_files:
            if eds:
                EDS_INPUTS(f)
            if mrsf:
                MRSF_INPUTS(f)

#----------------------------------------------------------------------------
# Common input template routines
#----------------------------------------------------------------------------
class INPUT_TEMPLATE(Template):
    delimiter = '@'

class INPUTS:
    """Common input routines"""
    def __init__(self, data):
        self.data = data
        self.ReadTemplate()

    def ReadTemplate(self):
        """Read or punch standard template"""
        try:
            self.tmpl = open(self.pkg+'.tmpl','r').read()
            self.tmpl = INPUT_TEMPLATE(self.tmpl)
            self.WriteInputs()
        except IOError:
            print("There's no " + self.pkg + " template. I'm punching one - please check")
            open(self.pkg+'.tmpl','w').write(self.tmpl)
            sys.exit()
        except AttributeError:
            pass

    def WriteInputs(self):
        pass

class MRSF_INPUTS(INPUTS):
    """Gamess US input routines"""

    def __init__(self, data):
        # template name
        self.pkg = "mrsf"
        # template content
        self.tmpl="""\
 $contrl scftyp=rohf runtyp=conical dfttyp=bhhlyp icharg=0
         tddft=mrsf maxit=400 mult=3 nprint=7 units=angs
         icut=20 itol=30 qmttol=1.0e-05 ispher=1 $end
 $conicl ixroot(1)=1,2 $end
 $tddft nstate=3 iroot=1 $end
 $scf dirscf=.t. fdiff=.f. diis=.t. soscf=.f.
      conv=1d-7 ethrsh=2.0 swdiis=0.005 cuhf=.t.
      shift=.t. damp=.t. npreo(1)=1,-1,2,1 $end
 $punch molden=.true. $end
 $dft swoff=0d-0 switch=0d-0 nrad0=96 nleb0=302
      nrad=99 nleb=590 $end
 $tddft nrad=99 nleb=590 $end
 $guess guess=huckel $end
 $basis gbasis=ccd extfil=.f. $end
 $system timlim=999999100 mwords=100 $end
@data
"""
        INPUTS.__init__(self, data)

    def WriteInputs(self):

        # initialize periodic table
        p=Periodic(0)
        
        # read xyz file
        try:
            xyz=open(self.data,'r').readlines()
            xyz=xyz[2:int(xyz[0])+2]
        except ValueError:
            print("Problem with *.xyz file?")
            sys.exit(1)
        
        for i in range(len(xyz)):
            xyz[i]=xyz[i].split()
            xyz[i].insert(1, str( Atomn(xyz[i][0], p) ))
            xyz[i]='%-5s %5s %15s %15s %15s\n' % tuple([xyz[i][j] for j in range(5)])
        
        # Prepare $data block
        xyz.insert(0, ' $data\n%s\nc1 0\n' % self.data)
        xyz.append(' $end')
        xyz=''.join(xyz)

        # Set filename
        filename = self.data.replace('.xyz','_mrsf')
        filename = filename.replace(' ','_')
        #
        finput = self.tmpl.substitute(fname=filename, data=xyz)
        #
        open(filename+'.inp','w').write(finput)

class EDS_INPUTS(INPUTS):
    """Gamess US input routines"""

    def __init__(self, data):
        # template name
        self.pkg = "eds"
        # template content
        self.tmpl="""\
 $system mwords=10 memddi=0 parall=.t. $end
 $contrl scftyp=rhf runtyp=eds icharg=0 mult=1 units=angs
         maxit=100 nprint=8 exetyp=run mplevl=2 ispher=-1
         icut=20 itol=30 aimpac=.f. cctyp=none $end
 $eds    mch(1)=0,0 mmul(1)=1,1 mnr(1)=1 ffeds=.t. $end
 $mp2    mp2prp=.t. code=ddi cutoff=1d-20 $end
 $ccinp  iconv=14 maxcc=100 $end
 $trans  cuttrf=1d-15 $end
 $scf    dirscf=.t. fdiff=.f. diis=.t. soscf=.f.
         conv=1d-11 swdiis=0.0001 npreo(1)=1,-1,2,1 $end
 $basis  gbasis=sto ngauss=3 $end
@data
"""
        INPUTS.__init__(self, data)

    def WriteInputs(self):

        # initialize periodic table
        p=Periodic(0)
        
        # read xyz file
        try:
            xyz=open(self.data,'r').readlines()
            xyz=xyz[2:int(xyz[0])+2]
        except ValueError:
            print("Problem with *.xyz file?")
            sys.exit(1)
        
        for i in range(len(xyz)):
            xyz[i]=xyz[i].split()
            xyz[i].insert(1, str( Atomn(xyz[i][0], p) ))
            xyz[i]='%-5s %5s %15s %15s %15s\n' % tuple([xyz[i][j] for j in range(5)])
        
        # Prepare $data block
        xyz.insert(0, ' $data\n%s\nc1 0\n' % self.data)
        xyz.append(' $end')
        xyz=''.join(xyz)

        # Set filename
        filename = self.data.replace('.xyz','_eds')
        filename = filename.replace(' ','_')
        #
        finput = self.tmpl.substitute(fname=filename, data=xyz)
        #
        open(filename+'.inp','w').write(finput)


#----------------------------------------------------------------------------
# Utilities
#----------------------------------------------------------------------------

def Run(program, args):
    '''Run shell command'''
    os.system(program + ' ' + args)

def SkipLines(open_file,n):
    '''Read n lines from file f.'''

    for i in range(n):
        line = open_file.readline()
        if line == '' :
            break

    return line

def FindLine(open_file,pattern):
    '''Read lines until pattern matches.'''

    while 1:
        line = open_file.readline()
        if line.find(pattern) !=-1 :
            break
        if line == '' :
            line = -1
            break

    return line

def ReFindLine(open_file,pattern):
    '''Read lines until pattern matches.'''

    while 1:
        line = open_file.readline()
        if re.search(pattern,line) !=None : break
        if line == '' :
            line = -1
            break

    return line

def Periodic(mendeleiev):
    '''Returns the mendeleiev table as a python list of tuples. Each cell
    contains either None or a tuple (symbol, atomic number), or a list of pairs
    for the cells * and **. Requires: "import re". Source: Gribouillis at
    www.daniweb.com - 2008 '''

    # L is a consecutive list of tuples ('Symbol', atomic number)
    L = [ (e,i+1) for (i,e) in enumerate( re.compile ("[A-Z][a-z]*").findall('''
    HHeLiBeBCNOFNeNaMgAlSiPSClArKCaScTiVCrMnFeCoNiCuZnGaGeAsSeBrKr
    RbSrYZrNbMoTcRuRhPdAgCdInSnSbTeIXeCsBaLaCePrNdPmSmEuGdTbDyHoEr
    TmYbLuHfTaWReOsIrPtAuHgTlPbBiPoAtRnFrRaAcThPaUNpPuAmCmBkCfEsFm
    MdNoLrRfDbSgBhHsMtDsRgUubUutUuqUupUuhUusUuo'''))]

    # The following fills the void with nones and returns the list of lists
    mendeleiev = 0

    if mendeleiev:
        for i,j in ( (88,103), (56,71) ):
            L[i] = L[i:j]
            L[i+1:] = L[j:]
        for i,j in ( (12,10 ), (4,10), (1,16) ):
            L[i:i]=[None]*j 

        return [ L[18*i:18*(i+1)] for i in range (7) ]

    # Return a plain list of tuples
    else:
        return L

def Atomn(s,ptable):
    '''Returns the atomic number based on atomic symbol string
    ptable is a list of consecutive (symbol, atomic number) tuples.'''

    for n,a in enumerate(ptable):
        if a[0].lower().find(s.strip().lower()) !=-1 :
            return float(n+1)

def Atoms(n,ptable):
    '''Returns the atomic symbol based on atomic number
    ptable is a list of consecutive (symbol, atomic number) tuples.'''

    return ptable[int(n)-1][0]

#----------------------------------------------------------------------------
# Main routine
#----------------------------------------------------------------------------
if __name__ == "__main__":
    Main(sys.argv[1:])

