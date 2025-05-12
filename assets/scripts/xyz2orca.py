#!/usr/bin/env python3
"""
This script generates orca input files from xyz file provided based on the
available template file. If there's no template the default template shall be
generated. 

If -g (--grep) option is used the data files should be the log files from which
the energies are extracted.

Usage: xyz2orca.py [options] data file(s)

Options:
  -h, --help       show this help
  -d, --dft        prepare nevpt2 inputs
  -n, --nevpt2     prepare nevpt2 inputs
  -g, --grep       extract energies from logs
  -l, --nm         extract wavelengths in nm
  -x, --xyz        copy xyz files for HTML output
  -w, --www        output data in HTML
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

__author__  = "Robert Gora (robert.gora@pwr.wroc.pl)"
__version__ = list(filter(str.isdigit, "$Revision: 1 $"))

# Import necessary modules
import os, sys, getopt, re

from string import Template
from io import StringIO
import numpy as np

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
    nevpt2 = False
    dft    = False
    grep   = False
    nm     = False
    xyz    = False
    www    = False

    # Parse commandline
    try:
        opts, args = getopt.getopt(argv, "hndglxw",
                                        ["help",
                                         "nevpt2",
                                         "dft",
                                         "grep",
                                         "nm",
                                         "xyz",
                                         "www"
                                         ])
    except getopt.GetoptError as error:
        print(error)
        Usage()
    if not argv:
        Usage()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            Usage()
        elif opt in ("-n", "--nevpt2"):
            nevpt2=True
        elif opt in ("-d", "--dft"):
            dft=True
        elif opt in ("-g", "--grep"):
            grep=True
        elif opt in ("-l", "--nm"):
            nm=True
        elif opt in ("-x", "--xyz"):
            xyz=True
        elif opt in ("-w", "--www"):
            www=True

    if not (nevpt2 or dft):
        print("Choose at least one type of job")
        sys.exit()

    # Parse each data file (with xyz coords)
    data_files = args

    if grep:
        # Exctract energies
        if nevpt2:
            for f in data_files:
                l=NEVPT2_LOGS(f,nm)
            NEVPT2_LOGS.PrintEnergies(l)
        if dft:
            for f in data_files:
                l=DFT_LOGS(f,nm)
                if xyz: CopyXYZ(f)
            if www:
                DFT_LOGS.PrintHTML(l)
            else:
                DFT_LOGS.PrintEnergies(l)
    else:
        # Prepare inputs
        for f in data_files:
            if nevpt2:
                NEVPT2_INPUTS(f)
            if dft:
                DFT_INPUTS(f)

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

class DFT_INPUTS(INPUTS):
    """Molpro input routines"""

    def __init__(self, data):
        # template name
        self.pkg = "dft"
        # template content
        self.tmpl="""\
! ROKS PBE0 D3BJ def2-SVP def2/J UseSym CPCM(THF) TightSCF RIJCOSX 

%output
        print[p_mos] true
        print[p_basis] 5
end

* xyzfile 0 1 @xyz

"""
        INPUTS.__init__(self, data)

    def WriteInputs(self):

        # read XYZ input file
        try:
            # read atomic coordinates
            xyz=open(self.data,'r').readlines()
            xyz=xyz[2:int(xyz[0])+2]
            xyz=''.join(xyz).rstrip()
            #
        except ValueError:
            print("Problem with *.xyz file?")
            sys.exit(1)

        # Set filename
        filename = self.data.replace('.xyz','_dft')
        filename = filename.replace(' ','_')
        #
        finput = self.tmpl.substitute(fname=filename, data=xyz, xyz=self.data)
        #
        open(filename+'.inp','w').write(finput)

class NEVPT2_INPUTS(INPUTS):
    """Molpro input routines"""

    def __init__(self, data):
        # template name
        self.pkg = "nevpt2"
        # template content
        self.tmpl="""\
# DFT guess orbitals
! ROKS PBE def2-SV(P) def2-SV(P)/JK SmallPrint UNO UseSym CPCM(THF)
! TightSCF RIJCOSX 

%output
        print[p_mos] true
        print[p_basis] 5
end

* xyzfile 0 1 @xyz

# Run CASSCF and NEVPT2 calculations
$new_job
! def2-svp def2/J RIJCOSX def2-svp/C TightSCF UseSym RIJCOSX moread
#%moinp "dft_@fname.gbw"

#%scf
#    rotate
#        {15, 12, 90} # interchange orbs 15 and 12
#    end
#end

%output
        print[p_mos] true
        print[p_basis] 5
end

%casscf
    cistep csfci
    trafostep ri       # RI approximation for CASSCF and NEVPT2
    PTMethod SC_NEVPT2
    PTSettings
        QDType QD_VanVleck
    end
    nel ?              # act el
    norb ?             # act orbs
    mult 1,3
    irrep 0,1          # irreps for s and t
    nroots 4,2         # 4S and 2T
    weights [0] = 1,1,1,1
    weights [1] = 1,1
    switchstep diis    # diis, nr, soscf
    switchconv 0.08    # gradient tol.
    maxiter 200        # max macro iter.
    shiftup 2.0        # level shifts
    shiftdn 2.0
    printwf false      # print the wave function
    #freezeie 0.4      # keep cas until
                       # int-ext rot. contribute <40%
    #freezeactive 0.03 # keep almost doubly occ orbs unless
                       # their contribution is less than 3%
    #freezegrad 0.2    # keep frozencore orbitals as long as
                       # their contribution is less than 2%
end

* xyzfile 0 1 @xyz

"""
        INPUTS.__init__(self, data)

    def WriteInputs(self):

        # read XYZ input file
        try:
            # read atomic coordinates
            xyz=open(self.data,'r').readlines()
            xyz=xyz[2:int(xyz[0])+2]
            xyz=''.join(xyz).rstrip()
            #
        except ValueError:
            print("Problem with *.xyz file?")
            sys.exit(1)

        # Set filename
        filename = self.data.replace('.xyz','_nevpt2')
        filename = filename.replace(' ','_')
        #
        finput = self.tmpl.substitute(fname=filename, data=xyz, xyz=self.data)
        #
        open(filename+'.inp','w').write(finput)

#----------------------------------------------------------------------------
# Parser routines
#----------------------------------------------------------------------------
class LOGS:
    """A common parser routine"""

    # common variables
    re_number=re.compile(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?")
    eLabels = []
    Energies = {}
    Properties = {}
    lLen = 0  # filename
    eLen = 15 # energies
    oLen = 0  # excitation type

    def __init__(self, logfile, nm):
        # initialize
        self.logfile = logfile
        self.nm = nm
        self.E = {}
        self.P = {}
        self.Sep = ","
        self.au2ev = 27.211396132
        self.ev2nm = 1239.84193
        # parse files
        self.ParseFile()

    def ParseFile(self):
        """Read properties from this file."""
        pass

    def PrintEnergies(self):
        """Print properties from this file."""
        pass

class NEVPT2_LOGS(LOGS):
    """A nevpt2 parser routine"""

    def ParseFile(self):
        """Read properties from this file."""
        self.log = open(self.logfile)
        # read MCSCF energies
        line = FindLine(self.log,'CAS-SCF STATES FOR BLOCK')
        if line !=-1 :
            if line.split()[6] == '1':
                self.Slabel = 'S'
            else:
                self.Slabel = 'T'
            self.ReadMultiEnergies()
            LOGS.Energies[self.logfile] = self.E
            LOGS.Properties[self.logfile] = self.P
            LOGS.lLen = max(len(self.logfile), LOGS.lLen)
        else:
            print("Couldn't find CASSCF energies")

        # read RS2C energies
        line = FindLine(self.log,'1PROGRAM * RS2C')
        if line !=-1 :
            self.ReadRS2CEnergies()
            LOGS.Energies[self.logfile] = self.E
            LOGS.lLen = max(len(self.logfile), LOGS.lLen)
        else:
            print("Couldn't find RS2C excitation energies")

    def ReadMultiEnergies(self):
        '''Read excitation energies from current log'''
        #
        SkipLines(self.log,2)
        while 1:
            line = self.log.readline()
            if line == '':
                break
            if line.find('--------------') !=-1 :
                break

            match = re.search(r"ROOT\s+(\d+):\s+E=", line)
            if match:
                data = line.split()
                label = self.Slabel + data[1].replace(':','')
                en  = float(data[3])
                #
                self.E[label]=[en]
                #
                if not label in LOGS.eLabels:
                    LOGS.eLabels.append(label)

        # Return energies
        return 0

    def ReadRS2CEnergies(self):
        '''Read excitation energies from current log'''
        #
        while 1:
            line = self.log.readline()
            if line == '':
                break
            if line.find('*******') !=-1 :
                break

            if re.search('!RSPT2 STATE .+ Energy',line) !=None :
                data = line.split()
                label = data[2]
                en  = float(data[-1])
                #
                self.E[label].append(en)

        # Return energies
        return 0

    def PrintEnergies(self):
        """Print properties from this batch."""
        # Prepare data for printout
        lLen = LOGS.lLen
        eLen = LOGS.eLen

        # print total energies header
        log = '# Total energies from %s calculations:\n' % 'MCSCF'
        # print labels
        line = '%'+str(lLen)+'s'+self.Sep
        log += line % '#'.ljust(lLen)
        line = '%'+str(eLen)+'s'+self.Sep
        for l in LOGS.eLabels:
            log += line % l.center(eLen)
        log += '\n'
        # print energies
        for f, p in sorted(LOGS.Energies.items()):
            # print filename ...
            line = '%'+str(lLen)+'s'+self.Sep
            log += line % f.rjust(lLen)
            # ... and energies for all states
            line = '%'+str(eLen)+'.7'+'f'+self.Sep
            for s, e in sorted(p.items()):
                log += line % e[0]
            log += '\n'

        # print total energies header
        log += '\n'
        log += '# Total energies from %s calculations:\n' % 'NEVPT2'
        # print labels
        line = '%'+str(lLen)+'s'+self.Sep
        log += line % '#'.ljust(lLen)
        line = '%'+str(eLen)+'s'+self.Sep
        for l in LOGS.eLabels:
            log += line % l.center(eLen)
        log += '\n'
        # print energies
        line0 = '%'+str(eLen)+'s'+self.Sep
        for f, p in sorted(LOGS.Energies.items()):
            # print filename ...
            line = '%'+str(lLen)+'s'+self.Sep
            log += line % f.rjust(lLen)
            # ... and energies for all states
            line = '%'+str(eLen)+'.7'+'f'+self.Sep
            for s, e in sorted(p.items()):
                if len(e)==2:
                    log += line % e[1]
                else:
                    log += line0 % '-'

            log += '\n'

        print(log)

class DFT_LOGS(LOGS):
    """A TDDFT parser routine"""

    def ParseFile(self):
        """Read properties from this file."""
        try:
            self.log = open(self.logfile)
            while 1:
                line = self.log.readline()
                if line == '':
                    break
                # read SCF energies
                if line.find('TOTAL SCF ENERGY') !=-1 :
                    self.ReadSCFEnergies()
                # read TDDFT energies
                if line.find('ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC')!=-1 :
                    self.ReadTDDFTEnergies()
                if self.E:
                    LOGS.Energies[self.logfile] = self.E
                    LOGS.Properties[self.logfile] = self.P
                    LOGS.lLen = max(len(self.logfile), LOGS.lLen)
        except FileNotFoundError:
            pass

    def ReadSCFEnergies(self):
        '''Read DFT energies from current log'''
        #
        while 1:
            line = self.log.readline()
            if line == '':
                break
            if line.find('TIMINGS') !=-1 :
                break
            if re.search('Total Energy       :',line) !=None :
                data = line.split()
                label = 'E_S0'
                en  = float(data[3])
                #
                self.E[label]=en
                #
                if not label in LOGS.eLabels:
                    LOGS.eLabels.append(label)
                #
                break

    def ReadTDDFTEnergies(self):
        '''Read TDDFT energies from current log'''
        #
        while 1:
            line = self.log.readline()
            if line == '':
                break
            if line.find('ABSORPTION SPECTRUM VIA TRANSITION VELOCITY') !=-1 :
                break

            if re.search(r"(\d+-\d+[A-Z])\s+->\s+(\d+-\d+[A-Z])",line):
                data = line.split()
                if re.search(r"(\d+-1[A-Z])",line):
                    elabel = 'dE_S'+data[2]
                    tlabel = 'E_S'+data[2]
                if re.search(r"(\d+-3[A-Z])",line):
                    elabel = 'dE_T'+data[2]
                    tlabel = 'E_T'+data[2]
                en   = float(data[3])
                lamb = float(data[5])
                fosc = float(data[6])
                #
                self.E[tlabel] = self.E['E_S0']+en/self.au2ev
                if self.nm:
                    self.E[elabel] = lamb
                else:
                    self.E[elabel] = en
                self.P[elabel] = fosc
                #
                if not elabel in LOGS.eLabels:
                    LOGS.eLabels.append(elabel)
                if not tlabel in LOGS.eLabels:
                    LOGS.eLabels.append(tlabel)

        # Return energies
        return 0

    def PrintEnergies(self):
        """Print properties from this batch."""
        # Prepare data for printout
        lLen = LOGS.lLen
        eLen = LOGS.eLen

        # print total energies header
        log = '# Total energies from %s calculations:\n' % 'TDDFT'
        # print labels
        line = '%'+str(lLen)+'s'+self.Sep
        lint = '%'+str(eLen)+'s'+self.Sep
        log += line % '#'.ljust(lLen)
        line = '%'+str(eLen)+'s'+self.Sep
        for l in LOGS.eLabels:
            log += line % l.rjust(eLen)
            log += line % 'f_osc'.rjust(eLen)
        log += '\n'
        # print energies
        for f, p in sorted(LOGS.Energies.items()):
            # print filename ...
            line = '%'+str(lLen)+'s'+self.Sep
            log += line % f.rjust(lLen)
            # ... and energies for all states
            line = '%'+str(eLen)+'.6'+'f'+self.Sep
            for l in LOGS.eLabels:
                try:
                    log += line % p[l]
                except KeyError:
                    log += lint % '-'.rjust(eLen)
                try:
                    log += line % LOGS.Properties[f][l]
                except KeyError:
                    log += lint % '-'.rjust(eLen)
            log += '\n'
            #line = '%'+str(lLen)+'s'+self.Sep
            #log += line % 'f_osc'.rjust(lLen)
            #line = '%'+str(eLen)+'.6'+'f'+self.Sep
            #for l in LOGS.eLabels:
            #    try:
            #        log += line % LOGS.Properties[f][l]
            #    except KeyError:
            #        log += lint % '-'.rjust(eLen)
            #log += '\n'

        print(log)

    def PrintHTML(self):
        """Print properties from this batch."""
        # Prepare data for printout
        lLen = LOGS.lLen
        eLen = LOGS.eLen

        # print total energies header
        log = '<!-- Total energies from %s calculations: -->\n' % 'TDDFT'
        log +='''<html>
<head>
<style>
#molecules {
  font-family: Arial, Helvetica, sans-serif;
  border-collapse: collapse;
  width: 100%;
}

#molecules td, #molecules th {
  border: 1px solid #ddd;
  padding: 8px;
}

#molecules tr:nth-child(even){background-color: #ffffff;}

#molecules tr:hover {background-color: #f2f2f2;}

#molecules th {
  padding-top: 12px;
  padding-bottom: 12px;
  text-align: left;
  background-color: #AAAAAA;
  color: white;
}
</style>
</head>
<body>
<table id="molecules">
'''
        # print labels
        log += '<tr>\n'
        line = '<td>%'+str(lLen)+'s</td>'
        lint = '<td>%'+str(eLen)+'s</td>'
        log += line % 'Molecule, E/au or dE/eV / f_osc'.ljust(lLen)
        line = '<td>%'+str(eLen)+'s</td>'
        for l in LOGS.eLabels:
            log += line % l.rjust(eLen)
        log += '</tr>\n'
        # print energies
        for f, p in sorted(LOGS.Energies.items()):
            # print filename ...
            line = '''<td rowspan="2">
<!-- data-surface='opacity:.7;color:white' -->
<div style="width: 400px; height: 400px; position: relative;"
  class="viewer_3Dmoljs"
  data-href="%s"
  data-type="xyz"
  data-style='stick'
  data-ui="true"
  data-backgroundcolor='0xffffff'
  data-backgroundalpha=0.1
></div>
%s
</td>
'''
            xyz_file = f.replace('.log','.xyz')
            xyz_file = xyz_file.replace(' ','_')
            xyz_file = xyz_file.replace('/','_')
            log += line % (xyz_file,xyz_file)
            # ... and energies for all states
            line = '<td>%'+str(eLen)+'.6'+'f</td>'
            for l in LOGS.eLabels:
                try:
                    log += line % p[l]
                except KeyError:
                    log += lint % '-'.rjust(eLen)
            log += '\n</tr>\n'
            log += '\n<tr>\n'
            line = '<td>%'+str(eLen)+'.6'+'f</td>'
            for l in LOGS.eLabels:
                try:
                    log += line % LOGS.Properties[f][l]
                except KeyError:
                    log += lint % '-'.rjust(eLen)
            log += '\n</tr>\n'

        log += '''</table>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
</body>
</html>'''

        print(log)
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

