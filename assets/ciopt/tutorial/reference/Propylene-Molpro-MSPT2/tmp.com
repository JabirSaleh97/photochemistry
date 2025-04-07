***ethylene SA-2-CAS(2/2)
! One step of S0/S1 conical intersection search in Cartesian coordinates
memory,50,m
FILE,2,uio.wfn
punch,all.pun
geometry={
nosymm,noorient,angstrom;
C,,        -0.00151839,      0.06643047,     -0.42582970
C,,        -0.44141403,      0.15996170,      0.94674777
C,,        -0.20459871,     -0.21946011,      2.32020699
H,,         0.90729111,     -0.55611433,     -0.42789595
H,,         0.27335644,      1.02284127,     -0.88884198
H,,        -0.71550447,     -0.42192727,     -1.10162612
H,,        -1.01411128,     -0.80550678,      2.77375182
H,,         0.69174033,     -0.86040700,      2.32572333
H,,         0.00475885,      0.61518185,      3.00175845
}
basis=6-31g**
 
INT
 
HF
 
MULTI
config,csf;
OCC,13;
CLOSED,11;
WF,24,1,0;
STATE,2;
WEIGHT,1,1;
MAXITER,100;
 
RS2,shift=.3,mix=2,root=  2
OCC,13;
CLOSED,11;
WF,24,1,0;
STATE,2;
 
forces
 
RS2,shift=.3,mix=2,root=  1
OCC,13;
CLOSED,11;
WF,24,1,0;
STATE,2;
 
forces
---
