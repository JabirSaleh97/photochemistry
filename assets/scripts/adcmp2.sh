#!/bin/bash
#
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

Define=0
Run=0
Help="0"
# Turbomole script defaults
Nnodes="1"
Ncores="4"
Version="7.6"
Arch="smp"
Memory="32"

if [ -z $1 ]
then
    Help="1"
else
    Help="0"
fi

# Read options
while getopts ":drhn:p:v:a:m:" Option
do
    case $Option in
        d ) Define=1 ;;
        r ) Run=1 ;;
        h ) Help="1" ;;
        n ) Nnodes=$OPTARG ;;
        p ) Ncores=$OPTARG ;;
        v ) Version=$OPTARG ;;
        a ) Arch=$OPTARG ;;
        m ) Memory=$OPTARG ;;
  esac
done

shift $(($OPTIND - 1))

export Turbomole="Turbomole -s -n $Nnodes -p $Ncores -v $Version -a $Arch -m $Memory"

# Print help
help() {

# Header
cat <<EOF
This script will prepare directories 'lower' and 'upper' for Turbomole
calculations and run MP2 and ADC(2) calculations, respectively. Please note
that the upper_def() and lower_def() functions for define might need to
be modified.

usage: `basename $1` inputfile.xyz

Where options are: 
  -d run define and prepare directories
  -r run ricc2/mp2 calculations
  -h print this help
  -m memory [$Memory]
  -n # nodes [$Nnodes]
  -p # cores [$Ncores]
  -v version [$Version]
EOF

exit 0
}

upper_def() {
cat <<EOF | $Turbomole define


a coord
sy c1
*
no
b all SV(P)
*
eht



scf
iter
500

cc
freeze
*
cbas
*
memory
4000
ricc2
model adc(2)
maxiter 100
geoopt adc(2) (a 1)
*
exci
irrep=a nexc=1
*
*
*
EOF
}

lower_def() {
cat <<EOF | $Turbomole define


a coord
sy c1
*
no
b all SV(P)
*
eht



scf
iter
500

cc
freeze
*
cbas
*
memory
2000
ricc2
model mp2
maxiter 100
geoopt mp2
*
*
*
*
EOF
}

[ $Help -eq 1 ] && help $0

# Run define for lower and upper state
run_define() {
    # Prepare coords
    $Turbomole x2t $1 > coord
    # Prepare lower dir
    rm -rf lower
    mkdir lower
    cp coord lower
    cd lower
    lower_def
    cd ..
    # Prepare upper dir
    rm -rf upper
    mkdir upper
    cp coord upper
    cd upper
    upper_def
    cd ..
}

# Run ricc2
run_ricc2() {
    # Redirect output
    run_define $1 > define.log
    # Check if prevous mos exist
    if [ -f mos ]
    then
        cat mos > lower/mos
        cat mos > upper/mos
    fi
    # run lower
    cd lower
    $Turbomole dscf | tee -a ../full.log > ../dscf.log
    if [ -f dscf_problem ]
    then
        rm -f dscf_problem
        $Turbomole dscf | tee -a ../full.log >> ../dscf.log
    fi
    if [ -f mos ]
    then
        cat mos > ../mos
        cat mos > ../upper/mos
    fi
    $Turbomole ricc2 | tee -a ../full.log > ../ricc2.log
    cat gradient > ../gradients
    cd ..
    # run upper
    cd upper
    $Turbomole dscf | tee -a ../full.log > ../dscf.log
    if [ -f dscf_problem ]
    then
        rm -f dscf_problem
        $Turbomole dscf | tee -a ../full.log >> ../dscf.log
    fi
    if [ -f mos ]
    then
        cat mos > ../mos
    fi
    $Turbomole ricc2 | tee -a ../full.log >> ../ricc2.log
    cat gradient >> ../gradients
    cd ..
}

if [ $Define -eq 1 ]
then
    run_define $1 > define.log
fi

if [ $Run -eq 1 ]
then
    run_ricc2 $1
fi

