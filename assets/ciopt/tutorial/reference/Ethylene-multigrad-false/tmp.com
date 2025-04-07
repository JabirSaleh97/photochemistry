***ethylene SA-2-CAS(2/2)
! One step of S0/S1 conical intersection search in Cartesian coordinates
memory,50,m
FILE,2,uio.wfn
punch,all.pun
geometry={
nosymm,noorient,angstrom;
C,,        -0.00420008,      0.08511219,     -0.43764419
C,,        -0.33680539,     -0.09715873,      0.94681648
H,,         1.06913893,     -0.16610413,     -0.52387791
H,,        -0.11059883,      1.11072186,     -0.80564349
H,,        -0.53205718,     -0.57947706,     -1.12921222
H,,        -0.08547743,     -0.35209415,      1.95056123
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
 
RS2,shift=.3,mix=2
OCC,9;
CLOSED,7;
WF,16,1,0;
STATE,2;
 
---
