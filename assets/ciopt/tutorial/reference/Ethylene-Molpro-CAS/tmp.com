***ethylene SA-2-CAS(2/2)
! One step of S0/S1 conical intersection search in Cartesian coordinates
memory,50,m
FILE,2,qwerty.wfn
punch,all.pun
geometry={
nosymm,noorient,angstrom;
C,,         0.01270752,      0.08316128,     -0.44951760
C,,        -0.31609960,     -0.10419905,      0.95662401
H,,         1.07575678,     -0.16119799,     -0.55815890
H,,        -0.11651152,      1.10580391,     -0.79699203
H,,        -0.53581374,     -0.57577321,     -1.11894023
H,,        -0.12003944,     -0.34679494,      1.96798438
}
basis=6-31g**
 
INT
HF
MULTI
config,csf;
OCC,9;
CLOSED,7;
WF,16,1,0;
STATE,2;
WEIGHT,1,1;
MAXITER,100;
CPMCSCF,GRAD,  2.1,record=5100.1
CPMCSCF,GRAD,  1.1,record=5101.1
FORCES
SAMC,5100.1
FORCES
SAMC,5101.1
---
