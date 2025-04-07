#!/usr/bin/env python3
"""
This script generates CIOpt templates from xyz file provided based on the
available template file. If there's no template the default template shall be
generated. If -g option is given the data_file should be the CIOpt log file.

Usage: xyz2ciopt.py [options] data_file

Options:
  -h, --help       show this help
  -r, --rs2        prepare molpro rs2  templates
  -c, --rs2c       prepare molpro rs2c templates
  -x, --xms        prepare molpro xms-rs2 templates
  -e, --eomcc      prepare molpro eom-ccsd templates
  -a, --adc2       prepare turbomole adc2 templates (requires adcmp2.sh in ~/bin dir)
                   (in case of CC2 change last line to:&%08(f20.10)00230)
  -l, --log        in case of adc2 read data from ricc2.log instead of gradient files
  -g, --grep       extract trajectory in xyz format from CIOpt.log
"""

#     Copyright (C) 2016, Robert W. Gora (robert.gora@pwr.wroc.pl)
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
import os, sys, getopt, re, subprocess

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
    '''Main driver routine'''

    data = {}

    # Set up defaults
    rs2 = False
    rs2c = False
    xms = False
    eomcc = False
    adc2 = False
    log = False
    grep = False

    # Parse commandline
    try:
        opts, args = getopt.getopt(argv, "hrcxealg",
                                        ["help",
                                         "rs2",
                                         "rs2c",
                                         "xms",
                                         "eomcc",
                                         "adc2",
                                         "log",
                                         "grep"
                                         ])
    except getopt.GetoptError as error:
        print(error)
        Usage()
    if not argv:
        Usage()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            Usage()
        elif opt in ("-r", "--rs2"):
            rs2=True
        elif opt in ("-c", "--rs2c"):
            rs2c=True
        elif opt in ("-x", "--xms"):
            xms=True
        elif opt in ("-e", "--eomcc"):
            eomcc=True
        elif opt in ("-a", "--adc2"):
            adc2=True
        elif opt in ("-l", "--log"):
            log=True
        elif opt in ("-g", "--grep"):
            grep=True

    if not (rs2 or rs2c or xms or eomcc or adc2 or grep):
        print("Choose at least one type of job")
        sys.exit()

    # Parse each data file (with xyz coords)
    data_files = args

    if grep:
        # Exctract trajectory
        for f in data_files:
            GetXYZ(f)
    else:
        # Prepare templates for CIOpt
        for f in data_files:
            if rs2:
                RS2_TEMPLATES(f)
            if rs2c:
                RS2C_TEMPLATES(f)
            if xms:
                XMS_TEMPLATES(f)
            if adc2 and log:
                ADC2_LOG_TEMPLATES(f)
            if adc2 and not log:
                ADC2_GR_TEMPLATES(f)

#----------------------------------------------------------------------------
# CIOpt templates
#----------------------------------------------------------------------------

class INPUT_TEMPLATE(Template):
    delimiter = '@'

class TEMPLATES:
    """Common input routines"""
    def __init__(self, data):
        self.data = data

    def ReadTemplate(self):
        """Read or punch standard template"""
        try:
            tmpl = open(self.pkg+'.tmpl','r').read()
            tmpl = INPUT_TEMPLATE(tmpl)
        except IOError:
            print("There's no " + self.pkg + " template. I'm punching one - please check")
            open(self.pkg+'.tmpl','w').write(self.tmpl)
            sys.exit()
        except AttributeError:
            pass

        return tmpl

    def WriteTemplates(self):
        pass

class RS2_TEMPLATES(TEMPLATES):
    """Molpro input templates"""

    def __init__(self, data):
        # init TEMPLATES
        TEMPLATES.__init__(self, data)
        # template content
        self.pkg = "rs2_e"
        self.tmpl="""\
***,rs2 for @fname;
memory,50,m
file,2,@fname.wfu
punch,@fname.dat

basis=avdz;                 ! aug-cc-pVDZ
gprint,basis;               ! print basis information
gprint,orbital,orbitals=10; ! Print +10 orbitals in SCF and MCSCF 

gthresh,twoint=1.d-12;  ! neglect 2e integrals (default 1.d-12)
gthresh,energy=1.d-6;   ! for energy (default 1.d-6)
gthresh,gradient=1.d-2; ! for orbital gradient in MCSCF (default 1.d-2)
gthresh,orbital=1.d-5;  ! for orbital optimization in the SCF (default 1.d-5)
gthresh,civec=1.d-5;    ! for CI coeffs in MCSCF and ref vec in CI (default 1.-d.5) 
gthresh,coeff=1.d-4;    ! for coefficients in CI and CCSD (default 1.d-4) 
gthresh,printci=0.05;   ! for printing CI coefficients (default 0.05) 
gthresh,symtol=1.d-6;   ! for finding symmetry equivalent atoms (default 1.d-6) 

angstrom
!symmetry, ?
nosym
orient, noorient

geometry={
@data}

{hf;orbital,2100.2}
put,molden,@fname.molden;orb,2100.2

{multi
occ,?,?;closed,?,?;     ! cas(?,?)
start,2100.2;           ! Start with RHF orbitals from above
!rotate,14.1,16.1       ! Rotate orbitals (angle 0 means interchange)
save,ref=4000.2;        ! Save configuration weights for CI in record 4000.2
wf,?,1,0;state,4;       ! SA 4 singlets from symm 1 
wf,?,2,0;state,2;       ! SA 2 singlets from symm 2
natorb,,ci;             ! Print natural orbitals and CI coefficients
canonical,,ci;          ! Print canonical orbitals and CI coefficients
maxit, 40;
} 

! ss-caspt2 with level shift and altered convergence thresholds
{rs2, shift=0.3, thrvar=1.0d-6, thrden=1.0d-6;
wf,?,1,0; state,2,1,5     ! for 2 singlets from symm 1
natorb,,ci;               ! print natural orbitals and CI coefficients
maxiter,200,300           ! set maximum macro and micro iterations
}

"""
        self.ReadTemplate()

        self.pkg = "rs2_g"
        self.tmpl="""\
***,rs2 for @fname;
memory,50,m
file,2,@fname.wfu
punch,@fname.dat

basis=avdz;                 ! aug-cc-pVDZ
gprint,basis;               ! print basis information
gprint,orbital,orbitals=10; ! Print +10 orbitals in SCF and MCSCF 

gthresh,twoint=1.d-12;  ! neglect 2e integrals (default 1.d-12)
gthresh,energy=1.d-6;   ! for energy (default 1.d-6)
gthresh,gradient=1.d-2; ! for orbital gradient in MCSCF (default 1.d-2)
gthresh,orbital=1.d-5;  ! for orbital optimization in the SCF (default 1.d-5)
gthresh,civec=1.d-5;    ! for CI coeffs in MCSCF and ref vec in CI (default 1.-d.5) 
gthresh,coeff=1.d-4;    ! for coefficients in CI and CCSD (default 1.d-4) 
gthresh,printci=0.05;   ! for printing CI coefficients (default 0.05) 
gthresh,symtol=1.d-6;   ! for finding symmetry equivalent atoms (default 1.d-6) 

angstrom
!symmetry, ?
nosym
orient, noorient

geometry={
@data}

{hf;orbital,2100.2}
put,molden,@fname.molden;orb,2100.2

{multi
occ,?,?;closed,?,?;     ! cas(?,?)
start,2100.2;           ! Start with RHF orbitals from above
!rotate,14.1,16.1       ! Rotate orbitals (angle 0 means interchange)
save,ref=4000.2;        ! Save configuration weights for CI in record 4000.2
wf,?,1,0;state,4;       ! SA 4 singlets from symm 1 
wf,?,2,0;state,2;       ! SA 2 singlets from symm 2
natorb,,ci;             ! Print natural orbitals and CI coefficients
canonical,,ci;          ! Print canonical orbitals and CI coefficients
maxit, 40;
} 

! ss-caspt2 with level shift and altered convergence thresholds
{rs2, shift=0.3, ignoreshift, thrvar=1.0d-6, thrden=1.0d-6, root=%#ISTATE;
wf,?,1,0; state,1,5       ! for 5th singlet from symm 1
natorb,,ci;               ! print natural orbitals and CI coefficients
maxiter,200,300           ! set maximum macro and micro iterations
}

forces

! ss-caspt2 with level shift and altered convergence thresholds
{rs2, shift=0.3, ignoreshift, thrvar=1.0d-6, thrden=1.0d-6, root=%#JSTATE;
wf,?,1,0; state,1,1       ! for 1st singlet from symm 1
natorb,,ci;               ! print natural orbitals and CI coefficients
maxiter,200,300           ! set maximum macro and micro iterations
}

forces

"""
        self.ReadTemplate()

        self.tmpl_read="""\
^002 Heff Half Sum symmetrized
^001 STATE       ENERGY     MIXING COEFFICIENTS
&%07(f13.7)00108
&%07(f13.7)00208
"""

        self.tmpl_control="""\
 &control
 nopt=3
 natoms=@natoms
 nstates=2
 istate=2
 jstate=1
 nefunc=7
 dlambdagap=3.5
 alpha=0.02
 tol=1.0d-06
 gtol=5.0d-03
 cigap=0.001
 znoncart=.false.
 zangrad=.true.
 zmultigrad=.true.
 crunstr='/opt/molpro/2012.1.25/bin/molpro tmp.com'
/
@data
"""
        self.tmpl_control=INPUT_TEMPLATE(self.tmpl_control)

        self.tmpl_readg="""\
^001 RSPT2 GRADIENT FOR STATE
@003
"""

        self.tmpl_readg2="""\
^002 RSPT2 GRADIENT FOR STATE
@003
"""
        self.pkg = "rs2_e"
        self.tmpl_write=self.ReadTemplate()
        self.pkg = "rs2_g"
        self.tmpl_writeg=self.ReadTemplate()
        self.WriteTemplates()

    def Coords(self):
        var=StringIO(self.xyz)
        var=np.genfromtxt(var, dtype=str)
        for i in range(len(var)):
            var[i][0]+=',, '
            for j in range(len(var[i])-1):
                var[i][j+1]='%%'+'%03d' % (i*3+j+1)
                if j<len(var[i])-2:
                    var[i][j+1]+=', '
        s=StringIO()
        np.savetxt(s,var,fmt='%s')
        s.seek(0)
        var=s.read()
        return var

    def Grads(self):
        var=StringIO(self.xyz)
        var=np.genfromtxt(var, dtype=str)
        var=var[:,1:]
        fmt=''
        for i in range(len(var)):
            fmt+='&%07(f15.9)' + '%03d10'   % (i*3+1)
            fmt+='%07(f15.9)'  + '%03d30'   % (i*3+2)
            fmt+='%07(f15.9)'  + '%03d49\n' % (i*3+3)
        return fmt

    def WriteTemplates(self):

        # read XYZ input file
        try:
            # read atomic coordinates
            xyz=open(self.data,'r').readlines()
            xyz=xyz[2:int(xyz[0])+2]
            self.xyz=''.join(xyz)
            #
        except ValueError:
            print("Problem with *.xyz file?")
            sys.exit(1)

        # Set filename
        filename = self.data.replace('.xyz','')
        # Write templates
        finput = self.tmpl_write.substitute(fname=filename, data=self.Coords())
        open('template.write','w').write(finput)
        #
        finput = self.tmpl_writeg.substitute(fname=filename, data=self.Coords())
        open('template.writeg','w').write(finput)
        #
        open('template.read','w').write(self.tmpl_read)
        #
        finput = self.tmpl_control.substitute(natoms=len(xyz), data=self.xyz)
        open('Control.dat','w').write(finput)
        #
        self.tmpl_readg=self.tmpl_readg+self.Grads()
        open('template.readg','w').write(self.tmpl_readg)
        #
        self.tmpl_readg2=self.tmpl_readg2+self.Grads()
        open('template.readg2','w').write(self.tmpl_readg2)

class RS2C_TEMPLATES(TEMPLATES):
    """Molpro input templates"""

    def __init__(self, data):
        # init TEMPLATES
        TEMPLATES.__init__(self, data)
        # template content
        self.pkg = "rs2_e"
        self.tmpl="""\
***,rs2 for @fname;
memory,50,m
file,2,@fname.wfu
punch,@fname.dat

basis=avdz;                 ! aug-cc-pVDZ
gprint,basis;               ! print basis information
gprint,orbital,orbitals=10; ! Print +10 orbitals in SCF and MCSCF 

gthresh,twoint=1.d-12;  ! neglect 2e integrals (default 1.d-12)
gthresh,energy=1.d-6;   ! for energy (default 1.d-6)
gthresh,gradient=1.d-2; ! for orbital gradient in MCSCF (default 1.d-2)
gthresh,orbital=1.d-5;  ! for orbital optimization in the SCF (default 1.d-5)
gthresh,civec=1.d-5;    ! for CI coeffs in MCSCF and ref vec in CI (default 1.-d.5) 
gthresh,coeff=1.d-4;    ! for coefficients in CI and CCSD (default 1.d-4) 
gthresh,printci=0.05;   ! for printing CI coefficients (default 0.05) 
gthresh,symtol=1.d-6;   ! for finding symmetry equivalent atoms (default 1.d-6) 

angstrom
!symmetry, ?
nosym
orient, noorient

geometry={
@data}

{hf;orbital,2100.2}
put,molden,@fname.molden;orb,2100.2

{multi
occ,?,?;closed,?,?;     ! cas(?,?)
start,2100.2;           ! Start with RHF orbitals from above
!rotate,14.1,16.1       ! Rotate orbitals (angle 0 means interchange)
save,ref=4000.2;        ! Save configuration weights for CI in record 4000.2
wf,?,1,0;state,4;       ! SA 4 singlets from symm 1 
wf,?,2,0;state,2;       ! SA 2 singlets from symm 2
natorb,,ci;             ! Print natural orbitals and CI coefficients
canonical,,ci;          ! Print canonical orbitals and CI coefficients
maxit, 40;
} 

! ss-caspt2 with level shift and altered convergence thresholds
{rs2c, shift=0.3, ignore, thrvar=1.0d-6, thrden=1.0d-6;
wf,?,1,0; state,2,1,5     ! for 2 singlets from symm 1
natorb,,ci;               ! print natural orbitals and CI coefficients
maxiter,200,300           ! set maximum macro and micro iterations
}

"""
        self.ReadTemplate()

        self.pkg = "rs2_g"
        self.tmpl="""\
***,rs2 for @fname;
memory,50,m
file,2,@fname.wfu
punch,@fname.dat

basis=avdz;                 ! aug-cc-pVDZ
gprint,basis;               ! print basis information
gprint,orbital,orbitals=10; ! Print +10 orbitals in SCF and MCSCF 

gthresh,twoint=1.d-12;  ! neglect 2e integrals (default 1.d-12)
gthresh,energy=1.d-6;   ! for energy (default 1.d-6)
gthresh,gradient=1.d-2; ! for orbital gradient in MCSCF (default 1.d-2)
gthresh,orbital=1.d-5;  ! for orbital optimization in the SCF (default 1.d-5)
gthresh,civec=1.d-5;    ! for CI coeffs in MCSCF and ref vec in CI (default 1.-d.5) 
gthresh,coeff=1.d-4;    ! for coefficients in CI and CCSD (default 1.d-4) 
gthresh,printci=0.05;   ! for printing CI coefficients (default 0.05) 
gthresh,symtol=1.d-6;   ! for finding symmetry equivalent atoms (default 1.d-6) 

angstrom
!symmetry, ?
nosym
orient, noorient

geometry={
@data}

{hf;orbital,2100.2}
put,molden,@fname.molden;orb,2100.2

{multi
occ,?,?;closed,?,?;     ! cas(?,?)
start,2100.2;           ! Start with RHF orbitals from above
!rotate,14.1,16.1       ! Rotate orbitals (angle 0 means interchange)
save,ref=4000.2;        ! Save configuration weights for CI in record 4000.2
wf,?,1,0;state,4;       ! SA 4 singlets from symm 1 
wf,?,2,0;state,2;       ! SA 2 singlets from symm 2
natorb,,ci;             ! Print natural orbitals and CI coefficients
canonical,,ci;          ! Print canonical orbitals and CI coefficients
maxit, 40;
} 

! ss-caspt2 with level shift and altered convergence thresholds
{rs2c, shift=0.3, ignore, thrvar=1.0d-6, thrden=1.0d-6, root=%#ISTATE;
wf,?,1,0; state,1,5       ! for 5th singlet from symm 1
!natorb,,ci;              ! print natural orbitals and CI coefficients
maxiter,200,300           ! set maximum macro and micro iterations
}

forces,numerical

! ss-caspt2 with level shift and altered convergence thresholds
{rs2c, shift=0.3, ignore, thrvar=1.0d-6, thrden=1.0d-6, root=%#JSTATE;
wf,?,1,0; state,1,1       ! for 1st singlet from symm 1
natorb,,ci;               ! print natural orbitals and CI coefficients
maxiter,200,300           ! set maximum macro and micro iterations
}

forces,numerical

"""
        self.ReadTemplate()

        self.tmpl_read="""\
^001 Correlation energy
&%07(f14.8)00136
^001 Correlation energy
&%07(f14.8)00236
"""

        self.tmpl_control="""\
 &control
 nopt=3
 natoms=@natoms
 nstates=2
 istate=2
 jstate=1
 nefunc=7
 dlambdagap=3.5
 alpha=0.02
 tol=1.0d-06
 gtol=5.0d-03
 cigap=0.001
 znoncart=.false.
 zangrad=.true.
 zmultigrad=.true.
 crunstr='/opt/molpro/2012.1.25/bin/molpro tmp.com'
/
@data
"""
        self.tmpl_control=INPUT_TEMPLATE(self.tmpl_control)

        self.tmpl_readg="""\
^001 Numerical gradient for RS2C
@004
"""

        self.tmpl_readg2="""\
^002 Numerical gradient for RS2C
@004
"""
        self.pkg = "rs2_e"
        self.tmpl_write=self.ReadTemplate()
        self.pkg = "rs2_g"
        self.tmpl_writeg=self.ReadTemplate()
        self.WriteTemplates()

    def Coords(self):
        var=StringIO(self.xyz)
        var=np.genfromtxt(var, dtype=str)
        for i in range(len(var)):
            var[i][0]+=',, '
            for j in range(len(var[i])-1):
                var[i][j+1]='%%'+'%03d' % (i*3+j+1)
                if j<len(var[i])-2:
                    var[i][j+1]+=', '
        s=StringIO()
        np.savetxt(s,var,fmt='%s')
        s.seek(0)
        var=s.read()
        return var

    def Grads(self):
        var=StringIO(self.xyz)
        var=np.genfromtxt(var, dtype=str)
        var=var[:,1:]
        fmt=''
        for i in range(len(var)):
            fmt+='&%07(f15.9)' + '%03d10'   % (i*3+1)
            fmt+='%07(f15.9)'  + '%03d30'   % (i*3+2)
            fmt+='%07(f15.9)'  + '%03d49\n' % (i*3+3)
        return fmt

    def WriteTemplates(self):

        # read XYZ input file
        try:
            # read atomic coordinates
            xyz=open(self.data,'r').readlines()
            xyz=xyz[2:int(xyz[0])+2]
            self.xyz=''.join(xyz)
            #
        except ValueError:
            print("Problem with *.xyz file?")
            sys.exit(1)

        # Set filename
        filename = self.data.replace('.xyz','')
        # Write templates
        finput = self.tmpl_write.substitute(fname=filename, data=self.Coords())
        open('template.write','w').write(finput)
        #
        finput = self.tmpl_writeg.substitute(fname=filename, data=self.Coords())
        open('template.writeg','w').write(finput)
        #
        open('template.read','w').write(self.tmpl_read)
        #
        finput = self.tmpl_control.substitute(natoms=len(xyz), data=self.xyz)
        open('Control.dat','w').write(finput)
        #
        self.tmpl_readg=self.tmpl_readg+self.Grads()
        open('template.readg','w').write(self.tmpl_readg)
        #
        self.tmpl_readg2=self.tmpl_readg2+self.Grads()
        open('template.readg2','w').write(self.tmpl_readg2)

class XMS_TEMPLATES(TEMPLATES):
    """Molpro input templates"""

    def __init__(self, data):
        # init TEMPLATES
        TEMPLATES.__init__(self, data)
        # template content
        self.pkg = "rs2_e"
        self.tmpl="""\
***,rs2 for @fname;
memory,50,m
file,2,@fname.wfu
punch,@fname.dat

basis=avdz;                 ! aug-cc-pVDZ
gprint,basis;               ! print basis information
gprint,orbital,orbitals=10; ! Print +10 orbitals in SCF and MCSCF 

gthresh,twoint=1.d-12;  ! neglect 2e integrals (default 1.d-12)
gthresh,energy=1.d-6;   ! for energy (default 1.d-6)
gthresh,gradient=1.d-2; ! for orbital gradient in MCSCF (default 1.d-2)
gthresh,orbital=1.d-5;  ! for orbital optimization in the SCF (default 1.d-5)
gthresh,civec=1.d-5;    ! for CI coeffs in MCSCF and ref vec in CI (default 1.-d.5) 
gthresh,coeff=1.d-4;    ! for coefficients in CI and CCSD (default 1.d-4) 
gthresh,printci=0.05;   ! for printing CI coefficients (default 0.05) 
gthresh,symtol=1.d-6;   ! for finding symmetry equivalent atoms (default 1.d-6) 

angstrom
!symmetry, ?
nosym
orient, noorient

geometry={
@data}

{hf;orbital,2100.2}
put,molden,@fname.molden;orb,2100.2

{multi
occ,?,?;closed,?,?;     ! cas(?,?)
start,2100.2;           ! Start with RHF orbitals from above
!rotate,14.1,16.1       ! Rotate orbitals (angle 0 means interchange)
save,ref=4000.2;        ! Save configuration weights for CI in record 4000.2
wf,?,1,0;state,4;       ! SA 4 singlets from symm 1 
wf,?,2,0;state,2;       ! SA 2 singlets from symm 2
natorb,,ci;             ! Print natural orbitals and CI coefficients
canonical,,ci;          ! Print canonical orbitals and CI coefficients
maxit, 40;
} 

! ms-caspt2 with level shift and altered convergence thresholds
{rs2, shift=0.3, xms=0, mix=2, thrvar=1.0d-6, thrden=1.0d-6;
wf,?,1,0; state,2         ! for 2 singlets from symm 1
natorb,,ci;               ! print natural orbitals and CI coefficients
maxiter,200,300           ! set maximum macro and micro iterations
}

"""
        self.ReadTemplate()

        self.pkg = "rs2_g"
        self.tmpl="""\
***,rs2 for @fname;
memory,50,m
file,2,@fname.wfu
punch,@fname.dat

basis=avdz;                 ! aug-cc-pVDZ
gprint,basis;               ! print basis information
gprint,orbital,orbitals=10; ! Print +10 orbitals in SCF and MCSCF 

gthresh,twoint=1.d-12;  ! neglect 2e integrals (default 1.d-12)
gthresh,energy=1.d-6;   ! for energy (default 1.d-6)
gthresh,gradient=1.d-2; ! for orbital gradient in MCSCF (default 1.d-2)
gthresh,orbital=1.d-5;  ! for orbital optimization in the SCF (default 1.d-5)
gthresh,civec=1.d-5;    ! for CI coeffs in MCSCF and ref vec in CI (default 1.-d.5) 
gthresh,coeff=1.d-4;    ! for coefficients in CI and CCSD (default 1.d-4) 
gthresh,printci=0.05;   ! for printing CI coefficients (default 0.05) 
gthresh,symtol=1.d-6;   ! for finding symmetry equivalent atoms (default 1.d-6) 

angstrom
!symmetry, ?
nosym
orient, noorient

geometry={
@data}

{hf;orbital,2100.2}
put,molden,@fname.molden;orb,2100.2

{multi
occ,?,?;closed,?,?;     ! cas(?,?)
start,2100.2;           ! Start with RHF orbitals from above
!rotate,14.1,16.1       ! Rotate orbitals (angle 0 means interchange)
save,ref=4000.2;        ! Save configuration weights for CI in record 4000.2
wf,?,1,0;state,4;       ! SA 4 singlets from symm 1 
wf,?,2,0;state,2;       ! SA 2 singlets from symm 2
natorb,,ci;             ! Print natural orbitals and CI coefficients
canonical,,ci;          ! Print canonical orbitals and CI coefficients
maxit, 40;
} 

! ms-caspt2 with level shift and altered convergence thresholds
{rs2, shift=0.3, xms=0, mix=2, thrvar=1.0d-6, thrden=1.0d-6, root=%#ISTATE;
wf,?,1,0; state,2         ! for 2 singlets from symm 1
natorb,,ci;               ! print natural orbitals and CI coefficients
maxiter,200,300           ! set maximum macro and micro iterations
}

forces

! ms-caspt2 with level shift and altered convergence thresholds
{rs2, shift=0.3, xms=0, mix=2, thrvar=1.0d-6, thrden=1.0d-6, root=%#JSTATE;
wf,?,1,0; state,2         ! for 2 singlets from symm 1
natorb,,ci;               ! print natural orbitals and CI coefficients
maxiter,200,300           ! set maximum macro and micro iterations
}

forces

"""
        self.ReadTemplate()

        self.tmpl_read="""\
^002 Heff Half Sum symmetrized
^001 STATE       ENERGY     MIXING COEFFICIENTS
&%07(f13.7)00108
&%07(f13.7)00208
"""

        self.tmpl_control="""\
 &control
 nopt=3
 natoms=@natoms
 nstates=2
 istate=2
 nefunc=7
 dlambdagap=3.5
 alpha=0.02
 tol=1.0d-06
 gtol=5.0d-03
 cigap=0.001
 znoncart=.false.
 zangrad=.true.
 zmultigrad=.true.
 crunstr='/opt/molpro/2012.1.25/bin/molpro tmp.com'
/
@data
"""
        self.tmpl_control=INPUT_TEMPLATE(self.tmpl_control)

        self.tmpl_readg="""\
^001 RSPT2 GRADIENT FOR STATE
@003
"""

        self.tmpl_readg2="""\
^002 RSPT2 GRADIENT FOR STATE
@003
"""
        self.pkg = "rs2_e"
        self.tmpl_write=self.ReadTemplate()
        self.pkg = "rs2_g"
        self.tmpl_writeg=self.ReadTemplate()
        self.WriteTemplates()

    def Coords(self):
        var=StringIO(self.xyz)
        var=np.genfromtxt(var, dtype=str)
        for i in range(len(var)):
            var[i][0]+=',, '
            for j in range(len(var[i])-1):
                var[i][j+1]='%%'+'%03d' % (i*3+j+1)
                if j<len(var[i])-2:
                    var[i][j+1]+=', '
        s=StringIO()
        np.savetxt(s,var,fmt='%s')
        s.seek(0)
        var=s.read()
        return var

    def Grads(self):
        var=StringIO(self.xyz)
        var=np.genfromtxt(var, dtype=str)
        var=var[:,1:]
        fmt=''
        for i in range(len(var)):
            fmt+='&%07(f15.9)' + '%03d10'   % (i*3+1)
            fmt+='%07(f15.9)'  + '%03d30'   % (i*3+2)
            fmt+='%07(f15.9)'  + '%03d49\n' % (i*3+3)
        return fmt

    def WriteTemplates(self):

        # read XYZ input file
        try:
            # read atomic coordinates
            xyz=open(self.data,'r').readlines()
            xyz=xyz[2:int(xyz[0])+2]
            self.xyz=''.join(xyz)
            #
        except ValueError:
            print("Problem with *.xyz file?")
            sys.exit(1)

        # Set filename
        filename = self.data.replace('.xyz','')
        # Write templates
        finput = self.tmpl_write.substitute(fname=filename, data=self.Coords())
        open('template.write','w').write(finput)
        #
        finput = self.tmpl_writeg.substitute(fname=filename, data=self.Coords())
        open('template.writeg','w').write(finput)
        #
        open('template.read','w').write(self.tmpl_read)
        #
        finput = self.tmpl_control.substitute(natoms=len(xyz), data=self.xyz)
        open('Control.dat','w').write(finput)
        #
        self.tmpl_readg=self.tmpl_readg+self.Grads()
        open('template.readg','w').write(self.tmpl_readg)
        #
        self.tmpl_readg2=self.tmpl_readg2+self.Grads()
        open('template.readg2','w').write(self.tmpl_readg2)

class ADC2_LOG_TEMPLATES(TEMPLATES):
    """Molpro input templates"""

    def __init__(self, data):
        # init TEMPLATES
        TEMPLATES.__init__(self, data)
        # template content
        self.pkg = "adc2"

        self.tmpl_write="""@natoms\n\n@data"""
        self.tmpl_write=INPUT_TEMPLATE(self.tmpl_write)

        self.tmpl_read="""\
^001     *   MP2 correlation energy
@001
&%08(F16.10)00153
^001 *<<<<<<<<<<<<<<<  EXCITED STATE PROPERTIES
^001  |    number, symmetry, multiplicity:    1 a    1
@003
&%08(F16.10)00240
"""

        self.tmpl_control="""\
 &control
 nopt=3
 natoms=@natoms
 nstates=2
 istate=2
 nefunc=7
 dlambdagap=3.5
 alpha=0.02
 tol=1.0d-06
 gtol=5.0d-03
 cigap=0.001
 znoncart=.false.
 zangrad=.true.
 zmultigrad=.true.
 cinpdeck='last.xyz'
 coutfile='ricc2.log'
 crunstr='./adcmp2.sh -v smp -n 1 -p 4 -m 800mb -r last.xyz '
/
@data
"""
        self.tmpl_control=INPUT_TEMPLATE(self.tmpl_control)

        self.tmpl_readg="""\
^001           total gradient of ADC(2) energy of excited state
^001           cartesian gradient of the energy (hartree/bohr)
@003
"""

        self.tmpl_readg2="""\
^001                   total gradient of MP2 energy
^001           cartesian gradient of the energy (hartree/bohr)
@003
"""
        self.WriteTemplates()

    def Coords(self):
        var=StringIO(self.xyz)
        var=np.genfromtxt(var, dtype=str)
        for i in range(len(var)):
            for j in range(len(var[i])-1):
                var[i][j+1]='%%'+'%03d' % (i*3+j+1)
        s=StringIO()
        np.savetxt(s,var,fmt='%s')
        s.seek(0)
        var=s.read()
        return var

    def Grads(self):
        var=StringIO(self.xyz)
        var=np.genfromtxt(var, dtype=str)
        var=var[:,1:]
        fmt=''
        # these come in 5 column rows with x,y,z elements 
        # organized colummnwise for each atom
        l=len(var)/5
        if l>0:
            # get gradients from 5 el. rows
            for i in range(l):
                for j in range(3):
                    fmt+='&%07(d15.7)' + '%03d06'   % (i*5+j+1)
                    fmt+='%07(d15.7)'  + '%03d21'   % (i*5+j+4)
                    fmt+='%07(d15.7)'  + '%03d36'   % (i*5+j+7)
                    fmt+='%07(d15.7)'  + '%03d51'   % (i*5+j+10)
                    fmt+='%07(d15.7)'  + '%03d66\n' % (i*5+j+13)
                fmt+='@002\n'
            # get the remaining gradients (if any)
            m=len(var)-l*5
            if m>0:
                for i in range(m):
                    for j in range(3):
                        if i==0:
                            fmt+='&%07(d15.7)' + '%03d06' % (l*15+i*m+j+1)
                        if i==1:
                            fmt+='%07(d15.7)'  + '%03d21' % (l*15+i*m+j+4)
                        if i==2:
                            fmt+='%07(d15.7)'  + '%03d36' % (l*15+i*m+j+7)
                        if i==3:
                            fmt+='%07(d15.7)'  + '%03d51' % (l*15+i*m+j+10)
                        if i==m-1:
                            fmt+='\n'
        else:
            # get the gradients
            m=len(var)
            for i in range(m):
                for j in range(3):
                    if i==0:
                        fmt+='&%07(d15.7)' + '%03d06' % (i*m+j+1)
                    if i==1:
                        fmt+='%07(d15.7)'  + '%03d21' % (i*m+j+4)
                    if i==2:
                        fmt+='%07(d15.7)'  + '%03d36' % (i*m+j+7)
                    if i==3:
                        fmt+='%07(d15.7)'  + '%03d51' % (i*m+j+10)
                    if i==m-1:
                        fmt+='\n'

        return fmt

    def WriteTemplates(self):
        """ write templates """

        # read XYZ input file
        try:
            # read atomic coordinates
            xyz=open(self.data,'r').readlines()
            xyz=xyz[2:int(xyz[0])+2]
            self.xyz=''.join(xyz)
        except ValueError:
            print("Problem with *.xyz file?")
            sys.exit(1)

        # Set filename
        filename = self.data.replace('.xyz','')

        # Write templates
        finput = self.tmpl_write.substitute(natoms=len(xyz), data=self.Coords())
        open('template.write','w').write(finput)
        open('template.writeg','w').write(finput)
        #
        open('template.read','w').write(self.tmpl_read)
        #
        finput = self.tmpl_control.substitute(natoms=len(xyz), data=self.xyz)
        open('Control.dat','w').write(finput)
        #
        self.tmpl_readg=self.tmpl_readg+self.Grads()
        open('template.readg','w').write(self.tmpl_readg)
        #
        self.tmpl_readg2=self.tmpl_readg2+self.Grads()
        open('template.readg2','w').write(self.tmpl_readg2)

        # Copy adcmp2.sh to $CWD
        try:
                subprocess.check_output('cp ~/bin/adcmp2.sh .', shell=True)
        except:
                print('There is no ~/bin/adcmp2.sh?')

class ADC2_GR_TEMPLATES(TEMPLATES):
    """Molpro input templates"""

    def __init__(self, data):
        # init TEMPLATES
        TEMPLATES.__init__(self, data)
        # template content
        self.pkg = "adc2"

        self.tmpl_write="""@natoms\n\n@data"""
        self.tmpl_write=INPUT_TEMPLATE(self.tmpl_write)

        self.tmpl_read="""\
^001$grad          cartesian gradients
&%08(f20.10)00130
^001$grad          cartesian gradients
&%08(f20.10)00233
"""

        self.tmpl_control="""\
 &control
 nopt=3
 natoms=@natoms
 nstates=2
 istate=2
 nefunc=7
 dlambdagap=3.5
 alpha=0.02
 tol=1.0d-06
 gtol=5.0d-03
 cigap=0.001
 znoncart=.false.
 zangrad=.true.
 zmultigrad=.true.
 cinpdeck='last.xyz'
 coutfile='gradients'
 crunstr='./adcmp2.sh -v smp -n 1 -p 4 -m 800mb -r last.xyz '
/
@data
"""
        self.tmpl_control=INPUT_TEMPLATE(self.tmpl_control)

        self.tmpl_readg="""\
^002$grad          cartesian gradients
"""

        self.tmpl_readg2="""\
^001$grad          cartesian gradients
"""
        self.WriteTemplates()

    def Coords(self):
        """ prepare coords template for writing """
        var=StringIO(self.xyz)
        var=np.genfromtxt(var, dtype=str)
        for i in range(len(var)):
            for j in range(len(var[i])-1):
                var[i][j+1]='%%'+'%03d' % (i*3+j+1)
        s=StringIO()
        np.savetxt(s,var,fmt='%s')
        s.seek(0)
        var=s.read()
        return var

    def Grads(self):
        """ prepare gradients template for writing """
        var=StringIO(self.xyz)
        var=np.genfromtxt(var, dtype=str)
        var=var[:,1:]
        fmt='@%03d\n' % (len(var)+1)
        # these come in rows with x,y,z elements 
        for i in range(len(var)):
            fmt+='&%08(f22.10)' + '%03d01'   % (i*3+1)
            fmt+='%08(f22.10)'  + '%03d23'   % (i*3+2)
            fmt+='%08(f22.10)'  + '%03d45\n' % (i*3+3)
        return fmt

    def WriteTemplates(self):
        """ write templates """

        # read XYZ input file
        try:
            # read atomic coordinates
            xyz=open(self.data,'r').readlines()
            xyz=xyz[2:int(xyz[0])+2]
            self.xyz=''.join(xyz)
        except ValueError:
            print("Problem with *.xyz file?")
            sys.exit(1)

        # Set filename
        filename = self.data.replace('.xyz','')

        # Write templates
        finput = self.tmpl_write.substitute(natoms=len(xyz), data=self.Coords())
        open('template.write','w').write(finput)
        open('template.writeg','w').write(finput)
        #
        open('template.read','w').write(self.tmpl_read)
        #
        finput = self.tmpl_control.substitute(natoms=len(xyz), data=self.xyz)
        open('Control.dat','w').write(finput)
        #
        self.tmpl_readg=self.tmpl_readg+self.Grads()
        open('template.readg','w').write(self.tmpl_readg)
        #
        self.tmpl_readg2=self.tmpl_readg2+self.Grads()
        open('template.readg2','w').write(self.tmpl_readg2)

        # Copy adcmp2.sh to $CWD
        try:
                subprocess.check_output('cp ~/bin/adcmp2.sh .', shell=True)
        except:
                print('There is no ~/bin/adcmp2.sh?')


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

    def __init__(self, logfile):
        # initialize
        self.logfile = logfile
        self.E = {}
        self.P = {}
        self.Sep = " "
        self.au2ev = 27.211396132
        # parse files
        self.ParseFile()

    def ParseFile(self):
        """Read properties from this file."""
        pass

    def PrintEnergies(self):
        """Print properties from this file."""
        pass

class EOMCC_LOGS(LOGS):
    """A eom-ccsd parser routine"""

    def ParseFile(self):
        """Read properties from this file."""
        self.log = open(self.logfile)
        # read energies
        line = FindLine(self.log,'1PROGRAM * CCSD')
        line = FindLine(self.log,'Final Results:')
        if line !=-1 :
            SkipLines(self.log,3)
            LOGS.Energies[self.logfile] = self.ReadEnergies()
            LOGS.lLen = max(len(self.logfile), LOGS.lLen)
        else:
            print("Couldn't find EOM-CCSD excitation energies")
        # read oscillator strengths
        line = FindLine(self.log,'EOM-CCSD properties program')
        if line !=-1 :
            line = FindLine(self.log,'Final Results for EOM-CCSD')
            LOGS.Energies[self.logfile] = self.ReadTransitions()
        else:
            print("Couldn't find EOM-CCSD properties")

    def ReadEnergies(self):
        '''Read excitation energies from current log'''
        self.E = {}
        #
        while 1:
            line = self.log.readline()
            if line == '':
                break
            if line.find('*****') !=-1 :
                break
            data = line.split()
            if len(data) >= 7:
                label = data[0]
                exc = float(data[1])
                en  = float(data[2])
                if len(data) == 7:
                    orb = ''
                else:
                    orb = ' '.join(data[-3:])
                LOGS.oLen = max(len(orb), LOGS.oLen)

                self.E[label]=[exc,en,orb]

                if not label in LOGS.eLabels:
                    LOGS.eLabels.append(label)
        # Return energies
        return self.E

    def ReadTransitions(self):
        '''Read excitation energies from current log'''
        #
        while 1:
            line = self.log.readline()
            if line == '':
                break
            if line.find('Final Results for EOM-CCSD') !=-1 :
                break
            if line.find('State    Exc.') !=-1 :
                l = SkipLines(self.log,1).split()[0]
                f = FindLine(self.log,'Oscillator strength').split()[-1]
                self.E[l].append(float(f))
        # Return energies
        return self.E

    def PrintEnergies(self):
        """Print properties from this batch."""
        # Prepare data for printout
        lLen = LOGS.lLen
        eLen = LOGS.eLen
        oLen = LOGS.oLen

        # print total energies header
        log = '# Total energies from %s calculations:\n' % 'EOM-CCSD'
        # print labels
        line = '%'+str(lLen)+'s'+self.Sep
        log += line % '#'.ljust(lLen)
        for l in LOGS.eLabels:
            log += line % l.center(eLen)
        log += '\n'
        # print energies
        for f, p in sorted(LOGS.Energies.items()):
            # print filename ...
            line = '%'+str(lLen)+'s'+self.Sep
            log += line % f.rjust(lLen)
            # ... and energies for all states
            line = '%'+str(eLen)+'.6'+'f'+self.Sep
            for s, e in sorted(p.items()):
                log += line % e[1]
            log += '\n'

        # print excitation energies header
        log += '\n'
        log += '# Ecitation energies from %s calculations:\n' % 'EOM-CCSD'
        # print labels
        line = '%'+str(lLen)+'s'+self.Sep
        log += line % '#'.ljust(lLen)
        for l in LOGS.eLabels:
            log += line % l.center(eLen)
        log += '\n'
        # print energies
        for f, p in sorted(LOGS.Energies.items()):
            # print filename ...
            line = '%'+str(lLen)+'s'+self.Sep
            log += line % f.rjust(lLen)
            # ... and energies for all states
            line = '%'+str(eLen)+'.3'+'f'+self.Sep
            for s, e in sorted(p.items()):
                log += line % (e[0]*self.au2ev)
            log += '\n'

        # print transition types
        log += '\n'
        log += '# Major contributions to transitions from %s calculations:\n' % 'EOM-CCSD'
        # print labels
        line = '%'+str(lLen)+'s'+self.Sep
        log += line % '#'.ljust(lLen)
        for l in LOGS.eLabels:
            log += line % l.center(eLen)
        log += '\n'
        # print transitions
        for f, p in sorted(LOGS.Energies.items()):
            # print filename ...
            line = '%'+str(lLen)+'s'+self.Sep
            log += line % f.rjust(lLen)
            # ... and transitions for all states
            line = '%'+str(eLen)+'s'+self.Sep
            for s, e in sorted(p.items()):
                log += line % e[2]
            log += '\n'

        # print transition strengths header
        log += '\n'
        log += '# Oscillator strengths from %s calculations:\n' % 'EOM-CCSD'
        # print labels
        line = '%'+str(lLen)+'s'+self.Sep
        log += line % '#'.ljust(lLen)
        for l in LOGS.eLabels:
            log += line % l.center(eLen)
        log += '\n'
        # print energies
        for f, p in sorted(LOGS.Energies.items()):
            # print filename ...
            line = '%'+str(lLen)+'s'+self.Sep
            log += line % f.rjust(lLen)
            # ... and f for all states
            for s, e in sorted(p.items()):
                if len(e) == 3:
                    line = '%'+str(eLen)+'s'+self.Sep
                    log += line % '-'
                else:
                    line = '%'+str(eLen)+'.3'+'f'+self.Sep
                    log += line % e[3]
            log += '\n'

        print(log)

class RS2C_LOGS(LOGS):
    """A rs2c parser routine"""

    def ParseFile(self):
        """Read properties from this file."""
        self.log = open(self.logfile)
        # read MCSCF energies
        line = FindLine(self.log,'1PROGRAM * MULTI')
        if line !=-1 :
            self.ReadMultiEnergies()
            LOGS.Energies[self.logfile] = self.E
            LOGS.Properties[self.logfile] = self.P
            LOGS.lLen = max(len(self.logfile), LOGS.lLen)
        else:
            print("Couldn't find MCSCF excitation energies")

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
        while 1:
            line = self.log.readline()
            if line == '':
                break
            if line.find('*******') !=-1 :
                break

            if re.search('!MCSCF STATE .+ Energy',line) !=None :
                data = line.split()
                label = data[2]
                en  = float(data[-1])
                #
                self.E[label]=[en]
                #
                if not label in LOGS.eLabels:
                    LOGS.eLabels.append(label)

            if line.find('!MCSCF trans ') !=-1 :
                data = line.split()
                label = data[2]
                mu  = float(data[3])
                #
                self.P[label]=mu

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
        for l in LOGS.eLabels:
            log += line % l.center(eLen)
        log += '\n'
        # print energies
        for f, p in sorted(LOGS.Energies.items()):
            # print filename ...
            line = '%'+str(lLen)+'s'+self.Sep
            log += line % f.rjust(lLen)
            # ... and energies for all states
            line = '%'+str(eLen)+'.6'+'f'+self.Sep
            for s, e in sorted(p.items()):
                log += line % e[0]
            log += '\n'

        # print total energies header
        log += '\n'
        log += '# Total energies from %s calculations:\n' % 'RS2C'
        # print labels
        line = '%'+str(lLen)+'s'+self.Sep
        log += line % '#'.ljust(lLen)
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
            line = '%'+str(eLen)+'.6'+'f'+self.Sep
            for s, e in sorted(p.items()):
                if len(e)==2:
                    log += line % e[1]
                else:
                    log += line0 % '-'

            log += '\n'

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

def GetXYZ(File):
    '''Get trajectory from CIOpt.log'''

    # open log file
    f=open(File,'r')

    XYZ=[]
    E=[]

    # loop through file
    while 1:
        l = f.readline()
        if l == '' :
            break

        if l.find('Iteration     0') !=-1 :
            label = l.split()[3]
            SkipLines(f,1)
            xyz = []
            while 1:
                l = f.readline()
                if l == '' :
                    break
                if re.search(r'\*\*\*',l) !=None :
                    xyz=np.array(xyz,dtype=str)
                    Natoms=len(xyz)
                    XYZ.append(xyz)
                    E.append(label)
                    break
                xyz.append(l.split())

        if l.find('Iteration') !=-1 :
            label = l.split()[3]

        if l.find('Current Geometry') !=-1 :
            xyz = []
            for i in range(Natoms):
                l=SkipLines(f,1)
                xyz.append(l.split())
            xyz=np.array(xyz,dtype=str)
            XYZ.append(xyz)
            E.append(label)

    # prepare trajectory
    out = ''
    for i in range(len(XYZ)):
        e=E[i]
        xyz=XYZ[i]
        out += '%d\n' % len(xyz)
        line =' %18.9f ' % float(e)
        out += line+'\n'
        for i in range(len(xyz)):
            line = '%5s %14.5f %14.5f %14.5f' % (xyz[i,0], float(xyz[i,1]), float(xyz[i,2]), float(xyz[i,3]))
            out += line+'\n'

    # save trajectory
    open(File.replace('.log','')+'.xyz','w').write(out)
    
#----------------------------------------------------------------------------
# Main routine
#----------------------------------------------------------------------------
if __name__ == "__main__":
    Main(sys.argv[1:])

